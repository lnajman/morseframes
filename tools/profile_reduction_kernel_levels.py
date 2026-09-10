#!/usr/bin/env python3
"""Optional RK level timelines and paired instrumentation-overhead diagnostics.

Construction is separate; fresh builder, pool, preparation and replay are timed.
Task IDs are persistent scratch-slot IDs, not physical CPU/thread identities.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import tempfile

from benchmark_reduction_kernel_ab import ROOT, Worker, command_output, header_digest, snapshot_headers, summarize
from profile_reduction_kernel_simplicial import validate as validate_global
from profile_reduction_kernel_simplicial import PARALLEL_CLOSURE_FIELDS, validate_parallel_closure

MODES = ('rk_plain', 'rk_coarse', 'rk_detailed', 'rk_levels_coarse', 'rk_levels_detailed')
TIMES = ('closure', 'facet', 'essential', 'facet_execution', 'aggregation', 'merge',
         'core', 'local_reduction', 'closure_initial', 'closure_packed',
         'closure_boundary_index', 'closure_traversal', 'closure_sort', 'closure_materialize')
COUNTS = ('kernel_rounds', 'facet_kernels', 'reductions', 'perforations', 'parallel_batches',
          'closure_sparse_cells', 'closure_sparse_entries', 'closure_boundary_visits',
          'closure_duplicate_faces', 'closure_boundary_index_visits', 'closure_boundary_index_entries',
          'local_candidate_visits', 'local_coboundary_visits', 'local_membership_tests',
          'local_membership_comparisons')
PHASES = tuple(name + '_seconds' for name in TIMES[:6])
TOL = 1e-10


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mode_orders(reverse=False):
    forward = [list(MODES[i:] + MODES[:i]) for i in range(len(MODES))]
    orders = forward + [row[::-1] for row in forward]
    return orders[::-1] if reverse else orders


def level_owners(values):
    if not values or any(not math.isfinite(v) for v in values) or len(set(values)) != len(values):
        raise ValueError('This focused diagnostic requires injective finite vertex values')
    return sorted(range(len(values)), key=values.__getitem__)


def validate(row, mode, identity, values, workers):
    flat = {k: v for k, v in row.items() if k != 'level_trace'}
    if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in flat.values()):
        raise ValueError('Invalid global values')
    if flat['algorithm_seconds'] <= 0 or not math.isclose(
            flat['algorithm_seconds'], flat['builder_seconds'] + flat['kernel_seconds'], rel_tol=1e-10):
        raise ValueError('Invalid outer timing partition')
    if mode != 'rk_plain':
        validate_global(flat)
    if not mode.startswith('rk_levels_'):
        if 'level_trace' in row:
            raise ValueError('Unexpected level trace')
        return
    trace = row['level_trace']
    owners = level_owners(values)
    detailed = mode == 'rk_levels_detailed'
    if (trace['completed'] is not True or trace['detailed'] is not detailed
            or trace['executor_workers'] != workers
            or trace['level_tasks'] != min(workers, len(owners))
            or len(trace['levels']) != len(owners)
            or not math.isfinite(trace['level_wall_seconds'])
            or not 0 <= trace['level_wall_seconds'] <= flat['level_wall_seconds'] + TOL):
        raise ValueError('Invalid trace header or coverage')
    timelines = [[] for _ in range(trace['level_tasks'])]
    for level, item in enumerate(trace['levels']):
        numeric = {k: v for k, v in item.items() if k != 'completed'}
        if any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in numeric.values()):
            raise ValueError('Invalid level value')
        if any(type(item[k]) is not int for k in ('level', 'task', 'simplices', 'events') + COUNTS):
            raise ValueError('Noninteger level count')
        if (item['completed'] is not True or item['level'] != level
                or not 0 <= item['task'] < len(timelines)
                or item['simplices'] != identity['lower_star_sizes'][owners[level]]
                or not 0 < item['events'] <= item['simplices']):
            raise ValueError('Invalid level identity, coverage or events')
        start, duration = item['start_seconds'], item['duration_seconds']
        if start + duration > trace['level_wall_seconds'] + TOL:
            raise ValueError('Level exceeds phase wall time')
        timelines[item['task']].append((start, start + duration))
        validate_parallel_closure(item, workers)
        if detailed:
            if (item['reductions'] + item['perforations'] != item['events']
                    or 2 * item['reductions'] + item['perforations'] != item['simplices']
                    or sum(item[k] for k in PHASES) > duration + TOL
                    or item['closure_packed_seconds'] + item['closure_boundary_index_seconds']
                       > item['closure_initial_seconds'] + TOL
                    or sum(item['closure_' + k + '_seconds'] for k in ('initial', 'traversal', 'sort', 'materialize'))
                       > item['closure_seconds'] + TOL):
                raise ValueError('Invalid per-level phase or event accounting')
            # This runner requires multiple injective levels: no nested facet workers.
            if len(owners) > 1 and item['core_seconds'] + item['local_reduction_seconds'] > item['facet_execution_seconds'] + TOL:
                raise ValueError('Facet child phases exceed execution time')
        elif any(item[k + '_seconds'] for k in TIMES) or any(item[k] for k in COUNTS):
            raise ValueError('Coarse trace unexpectedly ran detailed counters')
    for timeline in timelines:
        timeline.sort()
        if any(a[1] > b[0] + TOL for a, b in zip(timeline, timeline[1:])):
            raise ValueError('Overlapping work within one persistent task')
    if sum(r['simplices'] for r in trace['levels']) != identity['simplices']:
        raise ValueError('Missing simplices')
    if sum(r['events'] for r in trace['levels']) != identity['algorithms']['reduction_kernel']['steps']:
        raise ValueError('Missing events')
    parallel_fields = tuple(k for k in PARALLEL_CLOSURE_FIELDS if k in flat)
    if parallel_fields and any(not all(k in item for k in parallel_fields) for item in trace['levels']):
        raise ValueError('Missing per-level parallel closure profile')
    for key in tuple(k + '_seconds' for k in TIMES) + COUNTS + parallel_fields:
        global_key = 'rounds' if key == 'kernel_rounds' else key
        if global_key in flat and not math.isclose(sum(r[key] for r in trace['levels']), flat[global_key], rel_tol=1e-10, abs_tol=TOL):
            raise ValueError('Per-level/global mismatch: ' + key)


def overheads(blocks):
    return {f'{a}/{b}': summarize([{'baseline': [block['rows'][b]['algorithm_seconds']],
                                  'candidate': [block['rows'][a]['algorithm_seconds']]} for block in blocks])
            for a, b in [(m, 'rk_plain') for m in MODES[1:]] +
                        [('rk_levels_coarse', 'rk_coarse'), ('rk_levels_detailed', 'rk_detailed')]}


def audit(data):
    if data['schema'] != 'rk-level-profile-v1' or not data['completed'] or data['source_status']:
        raise ValueError('Incomplete/unfrozen profile')
    for path, expected in data['source_hashes'].items():
        content = subprocess.check_output(['git', 'show', data['source_revision'] + ':' + path], cwd=ROOT)
        if hashlib.sha256(content).hexdigest() != expected:
            raise ValueError('Source provenance mismatch')
    if digest(data['binary']) != data['binary_sha256']:
        raise ValueError('Binary hash mismatch')
    with tempfile.TemporaryDirectory(prefix='rk-level-audit-') as temp:
        snapshot = Path(temp) / 'source'
        snapshot_headers(data['source_revision'], snapshot)
        if header_digest(snapshot / 'include') != data['headers_sha256']:
            raise ValueError('Header provenance mismatch')
    seen = set()
    input_invariants = {}
    for case in data['cases']:
        key = (case['input'], case['workers'])
        if key in seen:
            raise ValueError('Duplicate case')
        seen.add(key)
        if digest(case['input']) != case['input_sha256'] or digest(case['reference_dump']) != case['reference_sha256']:
            raise ValueError('Input/reference dump hash mismatch')
        tokens = Path(case['input']).read_text().split()
        values = list(map(float, tokens[4:4 + int(tokens[2])]))
        if values != case['vertex_values']:
            raise ValueError('Vertex values do not match input')
        identity = case['metadata']['identity']
        if sum(identity['lower_star_sizes']) != identity['simplices']:
            raise ValueError('Invalid simplex identity')
        for algorithm in identity['algorithms'].values():
            criticals = algorithm['critical_counts']
            if (sum((-1)**i * n for i, n in enumerate(criticals)) != identity['euler']
                    or 2 * algorithm['steps'] - sum(criticals) != identity['simplices']):
                raise ValueError('Invalid critical counts')
        if len(case['blocks']) != len(mode_orders()):
            raise ValueError('Incomplete crossover')
        for block, order in zip(case['blocks'], mode_orders(data['reverse'])):
            if block['order'] != order or set(block['rows']) != set(MODES):
                raise ValueError('Invalid mode order')
            for mode, row in block['rows'].items():
                validate(row, mode, identity, values, case['workers'])
            detailed = block['rows']['rk_levels_detailed']['level_trace']['levels']
            invariant = (case['reference_sha256'], identity,
                         [[r[k] for k in ('simplices', 'events') + COUNTS] for r in detailed])
            if input_invariants.setdefault(case['input'], invariant) != invariant:
                raise ValueError('Per-level work changed between repetitions/worker counts')
        if case['overheads'] != overheads(case['blocks']):
            raise ValueError('Incorrect overhead summaries')
    if len(seen) != len(data['inputs']) * len(data['workers']):
        raise ValueError('Missing cases')
    expected = [(p, w) for p in data['inputs'] for w in data['workers']]
    if [(c['input'], c['workers']) for c in data['cases']] != expected:
        raise ValueError('Unexpected case order')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, nargs='+')
    p.add_argument('--workers', type=int, nargs='+', default=[1, 4, 8])
    p.add_argument('--warmups', type=int, default=2)
    p.add_argument('--reverse', action='store_true')
    p.add_argument('--output', type=Path)
    p.add_argument('--audit', type=Path, nargs='+')
    a = p.parse_args()
    if a.audit:
        studies = [json.loads(path.read_text()) for path in a.audit]
        for data in studies:
            audit(data)
        if len(studies) == 2:
            x, y = studies
            if x['source_revision'] != y['source_revision'] or x['inputs'] != y['inputs'][::-1] or x['workers'] != y['workers'][::-1] or x['reverse'] == y['reverse']:
                raise ValueError('Not a reversed frozen confirmation')
            def identities(d):
                return {(c['input'], c['workers']): (c['reference_sha256'], c['metadata']['identity']) for c in d['cases']}
            if identities(x) != identities(y):
                raise ValueError('Reference identity changed')
        print('Level trace audit passed:', ', '.join(map(str, a.audit)))
        return
    if (not a.inputs or not a.output or a.output.exists() or a.warmups < 1
            or min(a.workers) < 1 or len(set(a.workers)) != len(a.workers)
            or len(set(map(Path.resolve, a.inputs))) != len(a.inputs)):
        p.error('Unique inputs/workers, positive counts and a fresh output path required')
    status = command_output('git', 'status', '--porcelain')
    if status:
        p.error('Freeze all tracked changes before profiling')
    inputs = [str(path.resolve()) for path in a.inputs]
    workers = a.workers
    if a.reverse:
        inputs.reverse(); workers = workers[::-1]
    artifacts = a.output.resolve().with_suffix('.artifacts')
    artifacts.mkdir()
    binary = artifacts / 'rk-level-worker'
    source_paths = ['benchmarks/benchmark_simplicial_gradients.cpp', 'benchmarks/pls_profile.hpp',
                    'tools/profile_reduction_kernel_levels.py', 'tools/profile_reduction_kernel_simplicial.py',
                    'tools/benchmark_reduction_kernel_ab.py']
    compiler = shutil.which('clang++')
    flags = ['-std=c++17', '-O3', '-DNDEBUG', '-pthread']
    data = {'schema': 'rk-level-profile-v1', 'completed': False,
            'source_revision': command_output('git', 'rev-parse', 'HEAD'), 'source_status': status,
            'source_hashes': {path: digest(ROOT / path) for path in source_paths},
            'headers_sha256': header_digest(ROOT / 'include'),
            'compiler': command_output(compiler, '--version'), 'flags': flags,
            'platform': platform.platform(), 'architecture': platform.machine(), 'cpu_count': os.cpu_count(),
            'inputs': inputs, 'workers': workers, 'warmups': a.warmups, 'reverse': a.reverse, 'cases': []}
    with tempfile.TemporaryDirectory(prefix='rk-level-build-') as temp:
        snapshot = Path(temp) / 'source'
        snapshot_headers(data['source_revision'], snapshot)
        subprocess.run([compiler, *flags, '-I', str(snapshot / 'include'),
                        str(ROOT / source_paths[0]), '-o', str(binary)], check=True)
    data.update(binary=str(binary), binary_sha256=digest(binary), started_utc=datetime.now(timezone.utc).isoformat())
    with a.output.open('x') as output:
        for index, path in enumerate(inputs):
            path = Path(path)
            tokens = path.read_text().split()
            values = list(map(float, tokens[4:4 + int(tokens[2])]))
            level_owners(values)
            dump = artifacts / f'reference-{index}.dump'
            w = Worker(binary, path, dump)
            try:
                for count in workers:
                    def run(mode):
                        w.process.stdin.write(f'{mode} {count}\n'); w.process.stdin.flush()
                        return w.read()
                    for _ in range(a.warmups):
                        for mode in MODES:
                            run(mode)
                    blocks = []
                    for order in mode_orders(a.reverse):
                        rows = {mode: run(mode) for mode in order}
                        blocks.append({'order': order, 'rows': rows})
                    # Validation/statistical work is outside the timing block.
                    for block in blocks:
                        for mode, row in block['rows'].items():
                            validate(row, mode, w.metadata['identity'], values, count)
                    case = {'input': str(path), 'input_sha256': digest(path), 'vertex_values': values,
                            'workers': count, 'metadata': w.metadata, 'reference_dump': str(dump),
                            'reference_sha256': digest(dump), 'blocks': blocks, 'overheads': overheads(blocks)}
                    data['cases'].append(case)
                    print(path.name, count, 'RK plain ms', round(1000 * statistics.median(b['rows']['rk_plain']['algorithm_seconds'] for b in blocks), 3), flush=True)
            finally:
                w.close()
        data.update(completed=True, completed_utc=datetime.now(timezone.utc).isoformat())
        audit(data)
        json.dump(data, output, indent=2); output.write('\n')


if __name__ == '__main__':
    main()
