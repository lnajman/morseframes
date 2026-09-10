#!/usr/bin/env python3
"""Frozen old/new/toggle study with independent process-pair replication.

Both candidate modes use the SAME resident process and exact executable.
No profiling is interleaved with the timings. All raw runs are retained.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile

import benchmark_simplicial_gradients as common
import benchmark_rk_binary_isolation as isolation
from validate_simplicial_gradient_ab import equivalent
from render_rk_binary_isolation import recompute, check_equal

MODES = ('baseline', 'candidate_off', 'candidate_on')
PAIRS = (('candidate_on', 'baseline'), ('candidate_off', 'baseline'),
         ('candidate_on', 'candidate_off'))
SOURCES = tuple(dict.fromkeys(isolation.SOURCES + (
    'tools/benchmark_rk_parallel_closure.py', 'tools/render_rk_binary_isolation.py',
    'tools/profile_reduction_kernel_levels.py')))


def orders(reverse=False):
    rows = [list(p) for p in itertools.permutations(MODES)]
    return [p[::-1] for p in rows] if reverse else rows


def summarize(blocks):
    return {n + '/' + d: {m: {p: common.summarize([
        dict(candidate=[r[m][p] for r in b['runs'][n]],
             baseline=[r[m][p] for r in b['runs'][d]]) for b in blocks])
        for p in common.PHASES} for m in common.ALGORITHMS} for n, d in PAIRS}


def validate_case(case, data):
    assert len(case['blocks']) == 6
    assert case['launch_order'] == (['candidate', 'baseline'] if case['reverse'] else ['baseline', 'candidate'])
    identity = case['metadata']['baseline']['identity']
    assert identity == case['metadata']['candidate']['identity']
    assert sum(identity['lower_star_sizes']) == identity['simplices']
    for a in identity['algorithms'].values():
        assert sum((-1)**i * n for i, n in enumerate(a['critical_counts'])) == identity['euler']
        assert 2*a['steps'] - sum(a['critical_counts']) == identity['simplices']
    assert all(common.digest(Path(p)) == case['reference_sha256'] for p in case['dumps'].values())
    for block, order in zip(case['blocks'], orders(case['reverse'])):
        assert block['order'] == order and set(block['runs']) == set(MODES)
        for rows in block['runs'].values():
            common.check_runs(rows, 6)
            assert len({tuple(r) for r in rows}) == 6
    equivalent(case['summary'], summarize(case['blocks']))
    # Independent implementation: not the statistics function used by runners.
    for n, d in PAIRS:
        for m in common.ALGORITHMS:
            for p in common.PHASES:
                check_equal(case['summary'][n+'/'+d][m][p], recompute(case['blocks'], n, d, m, p))


def audit(path):
    data = json.loads(path.read_text())
    assert data['completed'] and data['schema'] == 'rk-parallel-closure-v1' and not data['source_status']
    assert set(data['source_sha256']) == set(SOURCES)
    for source, expected in data['source_sha256'].items():
        raw = subprocess.check_output(['git', 'show', data['source_revision'] + ':' + source], cwd=common.ROOT)
        assert hashlib.sha256(raw).hexdigest() == expected
    for v, build in data['builds'].items():
        assert common.digest(Path(build['binary'])) == build['binary_sha256']
        assert common.digest(Path(build['symbols'])) == build['symbols_sha256']
        assert isolation.trace_symbols(Path(build['symbols']).read_text()) == 0
        with tempfile.TemporaryDirectory(prefix='rk-closure-audit-') as tmp:
            root = Path(tmp) / 'snapshot'
            common.snapshot_headers(build['revision'], root)
            assert common.header_digest(root / 'include') == build['headers_sha256']
    settings = data['settings']
    assert settings['workers'] == ([8, 1] if settings['reverse'] else [1, 8])
    assert settings['process_pairs'] == 2 and settings['blocks'] == 6
    assert settings['repetitions'] == 6 and settings['warmups'] == 2
    expected = [(p, w, i) for p in settings['inputs'] for w in settings['workers'] for i in range(2)]
    assert [(c['input'], c['workers'], c['replicate']) for c in data['cases']] == expected
    references = {}
    for c in data['cases']:
        assert common.digest(Path(c['input'])) == c['input_sha256']
        assert c['reverse'] == bool(settings['reverse'] ^ (c['replicate'] % 2))
        validate_case(c, data)
        assert c['reference_sha256'] == references.setdefault(c['input'], c['reference_sha256'])
    print('VALIDATED', path.name, len(data['cases']), 'fresh process pairs;',
          len(data['cases'])*27, 'independently recomputed phase summaries;', common.digest(path), flush=True)
    return data


def run_mode(processes, mode, workers, repeats):
    process = processes['baseline' if mode == 'baseline' else 'candidate']
    command = 'run_closure_serial' if mode == 'candidate_off' else 'run'
    process.process.stdin.write(f'{command} {workers} {repeats}\n')
    process.process.stdin.flush()
    return process.read()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline')
    p.add_argument('--candidate')
    p.add_argument('--inputs', type=Path, nargs='+')
    p.add_argument('--reverse', action='store_true')
    p.add_argument('--output', type=Path)
    p.add_argument('--audit', type=Path, nargs='+')
    args = p.parse_args()
    if args.audit:
        studies = [audit(p) for p in args.audit]
        if len(studies) == 2:
            a, b = studies
            assert a['source_revision'] == b['source_revision'] and a['source_sha256'] == b['source_sha256']
            assert a['settings']['inputs'] == b['settings']['inputs'][::-1]
            assert a['settings']['workers'] == b['settings']['workers'][::-1]
            assert not a['settings']['reverse'] and b['settings']['reverse']
            for v in a['builds']:
                assert a['builds'][v]['revision'] == b['builds'][v]['revision']
                assert a['builds'][v]['headers_sha256'] == b['builds'][v]['headers_sha256']
            refs = lambda d: {c['input']: c['reference_sha256'] for c in d['cases']}
            assert refs(a) == refs(b)
            print('Frozen sources, reversed schedules and exact references agree')
        return
    if not args.baseline or not args.candidate or not args.inputs or not args.output or args.output.exists():
        p.error('Frozen baseline/candidate, inputs and fresh output required')
    status = common.command_output('git', 'status', '--porcelain')
    if status: p.error('Commit sources before timing')
    inputs = [str(p.resolve(strict=True)) for p in args.inputs]
    if len(inputs) != len(set(inputs)): p.error('Duplicate inputs')
    if args.reverse: inputs.reverse()
    workers = [8, 1] if args.reverse else [1, 8]
    artifacts = args.output.resolve().with_suffix('.artifacts'); artifacts.mkdir()
    compiler = shutil.which('clang++')
    flags = ['-std=c++17', '-O3', '-DNDEBUG', '-pthread', '-DMORSEFRAMES_BENCHMARK_ORDINARY_ONLY']
    data = dict(schema='rk-parallel-closure-v1', completed=False, source_status=status,
                source_revision=common.command_output('git','rev-parse','HEAD'),
                source_sha256={s: common.digest(common.ROOT/s) for s in SOURCES},
                compiler=common.command_output(compiler,'--version'), flags=flags,
                platform=platform.platform(), architecture=platform.machine(),
                settings=dict(inputs=inputs, workers=workers, reverse=args.reverse,
                              process_pairs=2, blocks=6, repetitions=6, warmups=2), builds={}, cases=[])
    for v, ref in [('baseline',args.baseline),('candidate',args.candidate)]:
        snapshot = artifacts / v
        revision = common.snapshot_headers(ref, snapshot)
        binary = artifacts / (v + '-ordinary')
        subprocess.run([compiler,*flags,'-I',str(snapshot/'include'),
                        str(common.ROOT/'benchmarks/benchmark_simplicial_gradients.cpp'),'-o',str(binary)], check=True)
        symbols = artifacts / (v + '.symbols'); symbols.write_text(isolation.read_symbols(binary))
        assert isolation.trace_symbols(symbols.read_text()) == 0
        data['builds'][v] = dict(revision=revision, binary=str(binary), binary_sha256=common.digest(binary),
            headers_sha256=common.header_digest(snapshot/'include'), symbols=str(symbols), symbols_sha256=common.digest(symbols))
        print('Built', v, 'ordinary-only executable', flush=True)
    data['started_utc'] = datetime.now(timezone.utc).isoformat()
    with args.output.open('x') as output:
        try:
            for ci, path in enumerate(inputs):
                for count in workers:
                    for replicate in range(2):
                        reverse = bool(args.reverse ^ (replicate % 2))
                        launch = ['candidate','baseline'] if reverse else ['baseline','candidate']
                        with ExitStack() as stack:
                            processes, dumps = {}, {}
                            for v in launch:
                                dump = artifacts / f'c{ci}-w{count}-r{replicate}-{v}.dump'
                                dumps[v] = str(dump)
                                processes[v] = common.Worker(Path(data['builds'][v]['binary']),Path(path),dump)
                                stack.callback(processes[v].close)
                            assert common.digest(Path(dumps['baseline'])) == common.digest(Path(dumps['candidate']))
                            for mode in orders(reverse)[0]: common.check_runs(run_mode(processes,mode,count,2),2)
                            blocks = []
                            for order in orders(reverse):
                                blocks.append(dict(order=order,runs={m:run_mode(processes,m,count,6) for m in order}))
                            data['cases'].append(dict(input=path,input_sha256=common.digest(Path(path)),workers=count,
                                replicate=replicate,reverse=reverse,launch_order=launch,dumps=dumps,
                                reference_sha256=common.digest(Path(dumps['baseline'])),
                                metadata={v:w.metadata for v,w in processes.items()},blocks=blocks))
                        print('Measured',Path(path).name,count,'workers, pair',replicate,flush=True)
            data['finished_utc'] = datetime.now(timezone.utc).isoformat()
            for case in data['cases']: case['summary'] = summarize(case['blocks'])
            data['completed'] = True
        finally:
            json.dump(data,output,indent=2); output.write('\n')
    audit(args.output)


if __name__ == '__main__': main()
