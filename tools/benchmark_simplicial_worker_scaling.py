#!/usr/bin/env python3
"""Frozen F-Max/PLS/RK worker crossover study; construction stays separate.

Each block contains every worker count. Each count runs all six method orders.
The second half reverses the cyclic worker schedules from the first half.
Use --reverse for a separately saved confirmation. Audit without new timings
with --audit MAIN CONFIRMATION. No algorithm or native-driver changes needed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import itertools
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import statistics
import subprocess
import tempfile

import benchmark_simplicial_gradients as common
from validate_simplicial_gradient_ab import equivalent

ROOT = common.ROOT
SOURCES = (
    'tools/benchmark_simplicial_worker_scaling.py',
    'tools/benchmark_simplicial_gradients.py',
    'tools/benchmark_reduction_kernel_ab.py',
    'tools/validate_simplicial_gradient_ab.py',
    'tools/profile_reduction_kernel_simplicial.py',
    'tools/pls_phase_profile.py',
    'benchmarks/benchmark_simplicial_gradients.cpp',
    'benchmarks/pls_profile.hpp',
)
STATIC = ('rounds', 'facet_kernels', 'facet_cell_visits', 'local_candidate_visits',
          'local_coboundary_visits', 'local_membership_tests',
          'local_membership_comparisons', 'local_large_membership_tests',
          'local_large_membership_comparisons', 'local_sparse_scan_passes',
          'local_sparse_candidate_visits', 'local_removed_candidate_visits',
          'local_protected_candidate_visits', 'closure_sparse_cells',
          'closure_sparse_entries', 'closure_boundary_visits', 'closure_duplicate_faces',
          'closure_boundary_index_visits', 'closure_boundary_index_entries')


def schedule(workers, blocks, reverse=False):
    n = len(workers)
    if len(set(workers)) != n or 1 not in workers or min(workers) < 1:
        raise ValueError('distinct positive worker counts including one required')
    if blocks < 2 * n or blocks % (2 * n):
        raise ValueError('blocks must be a positive multiple of twice the worker count')
    rows = []
    for i in range(blocks):
        offset = i % n
        row = workers[offset:] + workers[:offset]
        if (i // n) % 2:
            row = row[::-1]
        rows.append(row[::-1] if reverse else row)
    return rows


def summarize_case(blocks, workers):
    result = {}
    for count in workers:
        w = str(count)
        result[w] = {'methods': {}, 'relative_to_one': {}, 'comparison': {}}
        for method in common.ALGORITHMS:
            result[w]['methods'][method] = {
                phase: statistics.median(r[method][phase] for b in blocks for r in b['runs'][w])
                for phase in common.PHASES}
            paired = [{'baseline': [r[method]['algorithm_seconds'] for r in b['runs']['1']],
                       'candidate': [r[method]['algorithm_seconds'] for r in b['runs'][w]]}
                      for b in blocks]
            result[w]['relative_to_one'][method] = common.summarize(paired)
        for numerator, denominator in [('reduction_kernel', 'process_lower_stars'),
                                       ('reduction_kernel', 'f_max'),
                                       ('process_lower_stars', 'f_max')]:
            paired = [{'baseline': [r[denominator]['algorithm_seconds'] for r in b['runs'][w]],
                       'candidate': [r[numerator]['algorithm_seconds'] for r in b['runs'][w]]}
                      for b in blocks]
            result[w]['comparison'][f'{numerator}/{denominator}'] = common.summarize(paired)
    return result


def validate_profiles(profiles, identity, workers, repeats):
    counts = set()
    assert set(profiles) == set(map(str, workers))
    for w in workers:
        modes = profiles[str(w)]
        assert set(modes) == {'profile', 'rk_coarse', 'rk_detailed'}
        for mode, rows in modes.items():
            assert len(rows) == repeats
            for row in rows:
                if mode == 'profile':
                    assert row['stars'] == identity['vertices']
                    equivalent(row['validated_fine_seconds'], common.validate_pls_profile(
                        row['pls_profile_seconds'], row['setup_seconds'], row['local_wall_seconds'],
                        row['replay_seconds'], row['algorithm_seconds'] - row['builder_seconds']))
                    continue
                equivalent(row['unattributed_seconds'], common.validate_rk_profile(row))
                if mode == 'rk_detailed':
                    counts.add(tuple(row[k] for k in STATIC))
                if w > 1 and identity['vertices'] > 1:
                    n = min(w, identity['vertices'])
                    lo, hi, total = (row[k] for k in ('min_level_task_seconds',
                        'max_level_task_seconds', 'cumulative_level_task_seconds'))
                    assert lo <= hi <= row['level_wall_seconds'] + 1e-8
                    assert n * lo - 1e-8 <= total <= n * hi + 1e-8
                    assert row['min_worker_levels'] * n <= identity['vertices'] <= row['max_worker_levels'] * n
                    assert row['min_worker_simplices'] * n <= identity['simplices'] <= row['max_worker_simplices'] * n
                    assert row['level_chunks'] == math.ceil(identity['vertices'] / row['level_chunk_size'])
    assert len(counts) == 1, 'RK work changed across worker counts or repetitions'


def validate_case(case, settings):
    workers = settings['workers']
    expected = schedule(workers, settings['blocks'], settings['reverse'])
    assert len(case['blocks']) == len(expected)
    identity = case['metadata']['identity']
    assert sum(identity['lower_star_sizes']) == identity['simplices']
    for row in identity['algorithms'].values():
        assert sum((-1)**d * n for d, n in enumerate(row['critical_counts'])) == identity['euler']
    for block, order in zip(case['blocks'], expected):
        assert block['order'] == order and set(block['runs']) == set(map(str, workers))
        for runs in block['runs'].values():
            common.check_runs(runs, settings['repeats'])
            assert Counter(tuple(r) for r in runs) == {
                p: settings['repeats'] // 6 for p in itertools.permutations(common.ALGORITHMS)}
    equivalent(case['summary'], summarize_case(case['blocks'], workers))
    validate_profiles(case['profiles'], identity, workers, settings['profiles'])


def audit(path):
    data = json.loads(path.read_text())
    assert data['schema'] == 'simplicial-worker-scaling-v1' and data['completed']
    assert not data['source_status'] and not data['headers_patch']
    assert set(data['source_sha256']) == set(SOURCES)
    for source, sha in data['source_sha256'].items():
        # Validate against the recorded commit, not a potentially newer checkout.
        contents = subprocess.check_output(['git', 'show', f"{data['source_revision']}:{source}"], cwd=ROOT)
        import hashlib
        assert hashlib.sha256(contents).hexdigest() == sha, source
    with tempfile.TemporaryDirectory(prefix='scaling-audit-') as directory:
        target = Path(directory) / 'headers'
        assert common.snapshot_headers(data['build']['revision'], target) == data['build']['revision']
        assert common.header_digest(target / 'include') == data['build']['headers_sha256']
    assert len(data['cases']) == len(data['settings']['inputs'])
    assert [c['path'] for c in data['cases']] == data['settings']['inputs']
    assert len(set(data['settings']['inputs'])) == len(data['cases'])
    for case in data['cases']:
        assert common.digest(Path(case['path'])) == case['input_sha256']
        assert len(case['exact_reference_dump_sha256']) == 64
        validate_case(case, data['settings'])
    print('VALIDATED', path.name, len(data['cases']), 'inputs', common.digest(path), flush=True)
    return data


def audit_pair(paths):
    data = [audit(p) for p in paths]
    if len(data) == 2:
        a, b = data
        assert not a['settings']['reverse'] and b['settings']['reverse']
        assert a['source_revision'] == b['source_revision']
        assert a['source_sha256'] == b['source_sha256']
        assert a['build']['revision'] == b['build']['revision']
        assert a['build']['headers_sha256'] == b['build']['headers_sha256']
        for key in ('workers', 'blocks', 'repeats', 'warmups', 'profiles'):
            assert a['settings'][key] == b['settings'][key]
        assert a['settings']['inputs'] == b['settings']['inputs'][::-1]
        for x, y in zip(a['cases'], b['cases'][::-1]):
            for key in ('input_sha256', 'exact_reference_dump_sha256'):
                assert x[key] == y[key]
            assert x['metadata']['identity'] == y['metadata']['identity']
        print('VALIDATED reversed confirmation and exact ordered reference dumps')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit', nargs='+', type=Path)
    p.add_argument('--revision', default='HEAD')
    p.add_argument('--inputs', nargs='+', type=Path)
    p.add_argument('--workers', nargs='+', type=int, default=[1, 2, 4, 8])
    p.add_argument('--blocks', type=int, default=8)
    p.add_argument('--repeats', type=int, default=6)
    p.add_argument('--warmups', type=int, default=2)
    p.add_argument('--profiles', type=int, default=3)
    p.add_argument('--reverse', action='store_true')
    p.add_argument('--compiler', default='clang++')
    p.add_argument('--cxx-flags', default='-std=c++17 -O3 -DNDEBUG -pthread')
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    if a.audit:
        if len(a.audit) > 2:
            p.error('audit one study or a main/confirmation pair')
        audit_pair(a.audit)
        return
    if not a.inputs or not a.output or a.output.exists():
        p.error('inputs and a fresh output path are required')
    if a.repeats < 6 or a.repeats % 6 or min(a.warmups, a.profiles) < 1:
        p.error('repeats must be a positive multiple of six; warmups/profiles positive')
    orders = schedule(a.workers, a.blocks, a.reverse)
    inputs = [str(x.resolve(strict=True)) for x in a.inputs]
    if len(set(inputs)) != len(inputs):
        p.error('duplicate inputs')
    if a.reverse:
        inputs.reverse()
    status = common.command_output('git', 'status', '--porcelain')
    if status:
        p.error('freeze and commit runner/source changes before measuring')
    compiler = shutil.which(a.compiler)
    if not compiler:
        p.error('compiler unavailable')
    data = dict(schema='simplicial-worker-scaling-v1', completed=False,
        source_revision=common.command_output('git', 'rev-parse', 'HEAD'),
        source_status=status, headers_patch=common.command_output('git', 'diff', 'HEAD', '--', 'include'),
        started_utc=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
        architecture=platform.machine(), cpu_count=os.cpu_count(),
        compiler=common.command_output(compiler, '--version'), flags=shlex.split(a.cxx_flags),
        source_sha256={s: common.digest(ROOT / s) for s in SOURCES},
        settings=dict(inputs=inputs, workers=a.workers, blocks=a.blocks, repeats=a.repeats,
                      warmups=a.warmups, profiles=a.profiles, reverse=a.reverse),
        timing_scope='Same resident finalized native complex. Fresh builder, all gradient preparation, '
                     'pool startup, local work, replay, natural in-method cleanup. No persistence; '
                     'returned gradient/builder alive at stop. Loading/construction separate.',
        interval_scope='Paired-block bootstrap, within one session/input; not independent-machine '
                       'uncertainty. Worker ratios are t(workers)/t(1); F-Max stays sequential.',
        cases=[])
    with tempfile.TemporaryDirectory(prefix='simplicial-scaling-') as directory, a.output.open('x') as output:
        scratch = Path(directory)
        revision = common.snapshot_headers(a.revision, scratch / 'build')
        binary = scratch / 'worker'
        subprocess.run([compiler, *data['flags'], '-I', str(scratch / 'build/include'),
                        str(ROOT / SOURCES[-2]), '-o', str(binary)], check=True)
        data['build'] = dict(revision=revision,
            headers_sha256=common.header_digest(scratch / 'build/include'),
            binary_sha256=common.digest(binary), format=common.command_output('file', str(binary)))
        # All compilation/input generation is finished before measurement.
        for path in inputs:
            dump = scratch / 'reference.dump'
            worker = common.Worker(binary, Path(path), dump)
            try:
                case = dict(path=path, input_sha256=common.digest(Path(path)),
                    metadata=worker.metadata, exact_reference_dump_sha256=common.digest(dump),
                    blocks=[], profiles={str(w): {m: [] for m in ('profile', 'rk_coarse', 'rk_detailed')}
                                         for w in a.workers})
                for w in orders[0]:
                    common.check_runs(worker.run(w, a.warmups), a.warmups)
                for order in orders:
                    block = dict(order=order, runs={})
                    for w in order:
                        runs = worker.run(w, a.repeats)
                        common.check_runs(runs, a.repeats)
                        block['runs'][str(w)] = runs
                    case['blocks'].append(block)
                # Instrumented calls are separate, after all headline timings.
                for repeat in range(a.profiles):
                    for w in orders[repeat % len(orders)]:
                        for mode in list(case['profiles'][str(w)])[::(-1 if repeat % 2 else 1)]:
                            worker.process.stdin.write(f'{mode} {w}\n')
                            worker.process.stdin.flush()
                            row = worker.read()
                            if mode == 'profile':
                                row['validated_fine_seconds'] = common.validate_pls_profile(
                                    row['pls_profile_seconds'], row['setup_seconds'], row['local_wall_seconds'],
                                    row['replay_seconds'], row['algorithm_seconds'] - row['builder_seconds'])
                            else:
                                row['unattributed_seconds'] = common.validate_rk_profile(row)
                            case['profiles'][str(w)][mode].append(row)
            finally:
                worker.close()
            case['summary'] = summarize_case(case['blocks'], a.workers)
            validate_case(case, data['settings'])
            data['cases'].append(case)
            output.seek(0)
            json.dump(data, output, indent=2)
            output.write('\n'); output.truncate(); output.flush()
            print(Path(path).stem, 'RK ms', {w: round(1000 * row['methods']['reduction_kernel']['algorithm_seconds'], 3)
                  for w, row in case['summary'].items()}, flush=True)
        data['completed'] = True
        data['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.seek(0)
        json.dump(data, output, indent=2)
        output.write('\n'); output.truncate(); output.flush()


if __name__ == '__main__':
    main()
