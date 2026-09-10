#!/usr/bin/env python3
"""Old/new headers crossed with ordinary-only/diagnostic-capable native drivers.

Every measured call is unprofiled. Fresh processes for each input/worker-count
configuration remove prior profiled calls and worker-count history as confounders.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile

import benchmark_simplicial_gradients as common
from validate_simplicial_gradient_ab import equivalent

ROOT = common.ROOT
VARIANTS = ('baseline_ordinary', 'candidate_ordinary', 'baseline_full', 'candidate_full')
PAIRS = (('candidate_ordinary', 'baseline_ordinary'),
         ('candidate_full', 'baseline_full'),
         ('candidate_full', 'candidate_ordinary'),
         ('baseline_full', 'baseline_ordinary'))
SOURCES = ('benchmarks/benchmark_simplicial_gradients.cpp', 'benchmarks/pls_profile.hpp',
           'tools/benchmark_rk_binary_isolation.py', 'tools/benchmark_simplicial_gradients.py',
           'tools/benchmark_reduction_kernel_ab.py', 'tools/validate_simplicial_gradient_ab.py',
           'tools/profile_reduction_kernel_simplicial.py', 'tools/pls_phase_profile.py')


def read_symbols(binary):
    # Apple's nm does not implement GNU nm's -C; retain raw names instead.
    return subprocess.check_output(['nm', '-j', str(binary)], text=True)


def trace_symbols(symbols):
    return sum('build_reduction_kernel_implILb1E' in line or 'run_levels' in line
               for line in symbols.splitlines())


def schedule(blocks, reverse=False):
    if blocks < 8 or blocks % 8:
        raise ValueError('Use a positive multiple of eight blocks')
    rows = []
    for i in range(blocks):
        offset = i % 4
        row = list(VARIANTS[offset:] + VARIANTS[:offset])
        if (i // 4) % 2:
            row.reverse()
        rows.append(row[::-1] if reverse else row)
    return rows


def summaries(blocks):
    result = {}
    for numerator, denominator in PAIRS:
        key = numerator + '/' + denominator
        result[key] = {}
        for method in common.ALGORITHMS:
            result[key][method] = {}
            for phase in common.PHASES:
                samples = [{'baseline': [r[method][phase] for r in b['runs'][denominator]],
                            'candidate': [r[method][phase] for r in b['runs'][numerator]]} for b in blocks]
                result[key][method][phase] = common.summarize(samples)
    return result


def validate_case(case, settings):
    orders = schedule(settings['blocks'], settings['reverse'])
    if len(case['blocks']) != len(orders):
        raise ValueError('Incomplete blocks')
    expected_launch = list(VARIANTS[::-1] if settings['reverse'] else VARIANTS)
    if case['launch_order'] != expected_launch or set(case['metadata']) != set(VARIANTS):
        raise ValueError('Incorrect process coverage/order')
    identity = case['metadata'][VARIANTS[0]]['identity']
    if sum(identity['lower_star_sizes']) != identity['simplices']:
        raise ValueError('Invalid star-size accounting')
    for row in identity['algorithms'].values():
        counts = row['critical_counts']
        if (sum((-1)**i * n for i, n in enumerate(counts)) != identity['euler']
                or 2 * row['steps'] - sum(counts) != identity['simplices']):
            raise ValueError('Invalid critical/simplex accounting')
    for v in VARIANTS:
        if case['metadata'][v]['identity'] != identity:
            raise ValueError('Complex/gradient identity changed')
        if common.digest(Path(case['reference_dumps'][v])) != case['reference_sha256']:
            raise ValueError('Ordered reference dump changed')
    permutations = {p: settings['repeats'] // 6 for p in itertools.permutations(common.ALGORITHMS)}
    for block, order in zip(case['blocks'], orders):
        if block['order'] != order or set(block['runs']) != set(VARIANTS):
            raise ValueError('Invalid crossover coverage/order')
        for rows in block['runs'].values():
            common.check_runs(rows, settings['repeats'])
            if Counter(tuple(r) for r in rows) != permutations:
                raise ValueError('Unbalanced algorithm order')
    equivalent(case['summary'], summaries(case['blocks']))


def audit(path):
    data = json.loads(path.read_text())
    if data['schema'] != 'rk-binary-isolation-v1' or not data['completed'] or data['source_status']:
        raise ValueError('Incomplete/unfrozen study')
    if set(data['source_sha256']) != set(SOURCES):
        raise ValueError('Missing source provenance')
    for source, sha in data['source_sha256'].items():
        raw = subprocess.check_output(['git', 'show', data['source_revision'] + ':' + source], cwd=ROOT)
        if hashlib.sha256(raw).hexdigest() != sha:
            raise ValueError('Source provenance mismatch')
    for version in common.VERSIONS:
        with tempfile.TemporaryDirectory(prefix='rk-isolation-audit-') as temp:
            snapshot = Path(temp) / 'source'
            revision = common.snapshot_headers(data['revisions'][version], snapshot)
            sha = common.header_digest(snapshot / 'include')
            for mode in ('ordinary', 'full'):
                row = data['builds'][version + '_' + mode]
                if row['revision'] != revision or row['headers_sha256'] != sha:
                    raise ValueError('Header identity mismatch')
                if common.digest(Path(row['binary'])) != row['binary_sha256']:
                    raise ValueError('Binary identity mismatch')
                if common.digest(Path(row['symbols'])) != row['symbols_sha256']:
                    raise ValueError('Symbol evidence mismatch')
                if trace_symbols(Path(row['symbols']).read_text()) != row['trace_symbol_count']:
                    raise ValueError('Incorrect symbol count')
                if mode == 'ordinary' and row['trace_symbol_count']:
                    raise ValueError('Ordinary binary contains traced instantiation')
    settings = data['settings']
    expected = [(p, w) for p in settings['inputs'] for w in settings['workers']]
    if [(c['input'], c['workers']) for c in data['cases']] != expected:
        raise ValueError('Incorrect case order/coverage')
    for case in data['cases']:
        if common.digest(Path(case['input'])) != case['input_sha256']:
            raise ValueError('Input changed')
        validate_case(case, settings)
    print('VALIDATED', path.name, len(data['cases']), 'configurations', common.digest(path), flush=True)
    return data


def audit_pair(paths):
    studies = [audit(p) for p in paths]
    if len(studies) == 2:
        a, b = studies
        if (a['revisions'] != b['revisions'] or a['source_sha256'] != b['source_sha256']
                or a['source_revision'] != b['source_revision']
                or a['settings']['reverse'] or not b['settings']['reverse']
                or a['settings']['inputs'] != b['settings']['inputs'][::-1]
                or a['settings']['workers'] != b['settings']['workers'][::-1]):
            raise ValueError('Not a frozen reversed confirmation')
        for key in ('blocks', 'repeats', 'warmups'):
            if a['settings'][key] != b['settings'][key]:
                raise ValueError('Settings mismatch')
        identities = lambda d: {(c['input'], c['workers']): (c['reference_sha256'], c['metadata'][VARIANTS[0]]['identity']) for c in d['cases']}
        if identities(a) != identities(b):
            raise ValueError('Session reference identity mismatch')
        print('Reversed schedules and exact ordered identities agree')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline')
    p.add_argument('--candidate')
    p.add_argument('--inputs', type=Path, nargs='+')
    p.add_argument('--workers', type=int, nargs='+', default=[1, 8])
    p.add_argument('--blocks', type=int, default=8)
    p.add_argument('--repeats', type=int, default=6)
    p.add_argument('--warmups', type=int, default=2)
    p.add_argument('--reverse', action='store_true')
    p.add_argument('--output', type=Path)
    p.add_argument('--audit', type=Path, nargs='+')
    args = p.parse_args()
    if args.audit:
        audit_pair(args.audit); return
    schedule(args.blocks, args.reverse)
    if (not args.inputs or not args.baseline or not args.candidate or not args.output
            or args.output.exists() or args.repeats < 6 or args.repeats % 6 or args.warmups < 1
            or min(args.workers) < 1 or len(set(args.workers)) != len(args.workers)):
        p.error('Unique workers, positive counts, repeats divisible by six and fresh output required')
    inputs = [str(path.resolve(strict=True)) for path in args.inputs]
    if len(set(inputs)) != len(inputs):
        p.error('Duplicate inputs')
    workers = args.workers
    if args.reverse:
        inputs.reverse(); workers = workers[::-1]
    status = common.command_output('git', 'status', '--porcelain')
    if status:
        p.error('Freeze source changes before measurement')
    artifacts = args.output.resolve().with_suffix('.artifacts')
    artifacts.mkdir()
    compiler = shutil.which('clang++')
    flags = ['-std=c++17', '-O3', '-DNDEBUG', '-pthread']
    result = dict(schema='rk-binary-isolation-v1', completed=False, source_status=status,
                  source_revision=common.command_output('git', 'rev-parse', 'HEAD'),
                  source_sha256={p: common.digest(ROOT / p) for p in SOURCES},
                  compiler=common.command_output(compiler, '--version'), flags=flags,
                  platform=platform.platform(), architecture=platform.machine(), cpu_count=os.cpu_count(),
                  settings=dict(inputs=inputs, workers=workers, blocks=args.blocks,
                                repeats=args.repeats, warmups=args.warmups, reverse=args.reverse),
                  revisions={}, builds={}, cases=[])
    # All four builds finish before any timing starts. Preserve binaries/symbols.
    for version, ref in zip(common.VERSIONS, (args.baseline, args.candidate)):
        snapshot = artifacts / version
        revision = common.snapshot_headers(ref, snapshot)
        result['revisions'][version] = revision
        for mode in ('ordinary', 'full'):
            variant = version + '_' + mode
            binary = artifacts / variant
            defines = ['-DMORSEFRAMES_BENCHMARK_ORDINARY_ONLY'] if mode == 'ordinary' else []
            subprocess.run([compiler, *flags, *defines, '-I', str(snapshot / 'include'),
                            str(ROOT / SOURCES[0]), '-o', str(binary)], check=True)
            symbols = read_symbols(binary)
            symbol_file = artifacts / (variant + '.symbols')
            symbol_file.write_text(symbols)
            traces = trace_symbols(symbols)
            if mode == 'ordinary' and traces:
                raise ValueError('Unexpected trace instantiation in ordinary binary')
            result['builds'][variant] = dict(revision=revision, binary=str(binary), binary_sha256=common.digest(binary),
                headers_sha256=common.header_digest(snapshot / 'include'), defines=defines,
                format=common.command_output('file', str(binary)), bytes=binary.stat().st_size,
                symbols=str(symbol_file), symbols_sha256=common.digest(symbol_file), trace_symbol_count=traces)
            print('Built', variant, 'trace symbols', traces, flush=True)
    result['started_utc'] = datetime.now(timezone.utc).isoformat()
    with args.output.open('x') as output:
        try:
            for case_index, path in enumerate(inputs):
                for count in workers:
                    with ExitStack() as stack:
                        launch = list(VARIANTS[::-1] if args.reverse else VARIANTS)
                        processes, dumps = {}, {}
                        for variant in launch:
                            dump = artifacts / f'case-{case_index}-w{count}-{variant}.dump'
                            dumps[variant] = str(dump)
                            w = common.Worker(Path(result['builds'][variant]['binary']), Path(path), dump)
                            stack.callback(w.close); processes[variant] = w
                        fingerprints = {common.digest(Path(path)) for path in dumps.values()}
                        if len(fingerprints) != 1:
                            raise ValueError('Complete ordered reference dumps differ')
                        for variant in launch:
                            common.check_runs(processes[variant].run(count, args.warmups), args.warmups)
                        blocks = []
                        for order in schedule(args.blocks, args.reverse):
                            runs = {v: processes[v].run(count, args.repeats) for v in order}
                            blocks.append(dict(order=order, runs=runs))
                        result['cases'].append(dict(input=path, input_sha256=common.digest(Path(path)), workers=count,
                            launch_order=launch, metadata={v: w.metadata for v, w in processes.items()},
                            reference_dumps=dumps, reference_sha256=fingerprints.pop(), blocks=blocks))
                        print('Measured', Path(path).name, count, 'workers', flush=True)
            result['finished_utc'] = datetime.now(timezone.utc).isoformat()
            # Summaries/audits run only after all timing and process teardown.
            for case in result['cases']:
                case['summary'] = summaries(case['blocks'])
                validate_case(case, result['settings'])
            result['completed'] = True
        finally:
            json.dump(result, output, indent=2); output.write('\n')
    audit(args.output)


if __name__ == '__main__':
    main()
