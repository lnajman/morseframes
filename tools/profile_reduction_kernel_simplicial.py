#!/usr/bin/env python3
"""Coarse versus detailed RK diagnostics on resident simplicial complexes."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import statistics
import tempfile

from benchmark_reduction_kernel_ab import Worker, ROOT, command_output, header_digest

SEARCH_FIELDS = (
    'local_membership_comparisons', 'local_large_membership_tests',
    'local_large_membership_comparisons', 'local_sparse_scan_passes',
    'local_sparse_candidate_visits', 'local_removed_candidate_visits',
    'local_protected_candidate_visits',
)
CLOSURE_TIMES = ('closure_initial_seconds', 'closure_packed_seconds',
                 'closure_traversal_seconds', 'closure_sort_seconds',
                 'closure_materialize_seconds')
CLOSURE_COUNTS = ('closure_sparse_cells', 'closure_sparse_entries',
                  'closure_boundary_visits', 'closure_duplicate_faces',
                  'closure_index_growths', 'closure_entry_growths')
CLOSURE_FIELDS = CLOSURE_TIMES + CLOSURE_COUNTS
BOUNDARY_INDEX_FIELDS = ('closure_boundary_index_seconds', 'closure_boundary_index_visits',
                        'closure_boundary_index_entries')


def validate(row):
    if any(not math.isfinite(v) or v < 0 for v in row.values()):
        raise ValueError("Invalid RK diagnostic value")
    if not math.isclose(row['algorithm_seconds'], row['builder_seconds'] + row['kernel_seconds'], rel_tol=1e-10):
        raise ValueError("Outer RK phases do not partition algorithm time")
    remainder = row['kernel_seconds'] - sum(row[k] for k in ['setup_seconds','level_wall_seconds','replay_seconds'])
    if remainder < -1e-10:
        raise ValueError("Coarse RK phases exceed kernel time")
    if any(k in row for k in SEARCH_FIELDS):
        if not all(k in row for k in SEARCH_FIELDS):
            raise ValueError('Incomplete RK local-search counters')
        if (row['local_large_membership_tests'] > row['local_membership_tests']
                or row['local_large_membership_comparisons'] > row['local_membership_comparisons']
                or row['local_membership_comparisons'] < row['local_membership_tests']
                or row['local_large_membership_comparisons'] < row['local_large_membership_tests']
                or row['local_sparse_candidate_visits'] > row['local_candidate_visits']
                or row['local_removed_candidate_visits'] + row['local_protected_candidate_visits']
                   > row['local_sparse_candidate_visits']):
            raise ValueError('Inconsistent RK local-search counters')
    if any(k in row for k in CLOSURE_FIELDS):
        if not all(k in row for k in CLOSURE_FIELDS):
            raise ValueError('Incomplete RK closure profile')
        # Packed preparation is nested inside initial setup, not additive.
        children = sum(row[k] for k in CLOSURE_TIMES if k != 'closure_packed_seconds')
        if (row['closure_packed_seconds'] > row['closure_initial_seconds'] + 1e-10
                or children > row['closure_seconds'] + 1e-10
                or row['closure_sparse_cells'] > row['closure_sparse_entries']
                or row['closure_duplicate_faces'] > row['closure_boundary_visits']
                or row['closure_index_growths'] > row['closure_sparse_entries']
                or row['closure_entry_growths'] > row['closure_sparse_entries']):
            raise ValueError('Inconsistent RK closure phases or counters')
        same_level_visits = (row['closure_sparse_entries'] - row['closure_sparse_cells']
                             + row['closure_duplicate_faces'])
        if row['closure_boundary_visits'] < same_level_visits:
            raise ValueError('Missing RK closure traversal visits')
    if any(k in row for k in BOUNDARY_INDEX_FIELDS):
        if not all(k in row for k in BOUNDARY_INDEX_FIELDS + CLOSURE_FIELDS):
            raise ValueError('Incomplete RK boundary-index profile')
        if (row['closure_boundary_index_seconds'] + row['closure_packed_seconds']
                > row['closure_initial_seconds'] + 1e-10
                or row['closure_boundary_index_entries'] > row['closure_boundary_index_visits']):
            raise ValueError('Inconsistent RK boundary-index preparation')
        if row['closure_boundary_index_visits'] and row['closure_boundary_visits'] != same_level_visits:
            raise ValueError('Indexed RK traversal must visit only same-level boundaries')
    # Detailed kernel fields sum over levels/tasks, not global elapsed time.
    # Core/local are children of facet execution and must not be added to it.
    return max(0., remainder)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--inputs',type=Path,nargs='+',required=True)
    p.add_argument('--workers',type=int,nargs='+',default=[1,8])
    p.add_argument('--repeats',type=int,default=5)
    p.add_argument('--warmups',type=int,default=2)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if min(a.workers)<1 or a.repeats<1 or a.warmups<1 or a.output.exists():
        p.error('positive counts and a fresh output path are required')
    digest=lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    data={'schema':'rk-simplicial-profile-v1','started_utc':datetime.now(timezone.utc).isoformat(),
          'source_revision':command_output('git','rev-parse','HEAD'),
          'source_status':command_output('git','status','--porcelain'),
          'headers_sha256':header_digest(ROOT/'include'),'binary_sha256':digest(a.binary),
          'driver_sha256':digest(ROOT/'benchmarks/benchmark_simplicial_gradients.cpp'),
          'repeats':a.repeats,'warmups':a.warmups,'workers':a.workers,'cases':[]}
    with a.output.open('x') as output, tempfile.TemporaryDirectory(prefix='rk-profile-') as directory:
        for path in a.inputs:
            w=Worker(a.binary.resolve(),path.resolve(),Path(directory)/'reference.dump')
            try:
                for count in a.workers:
                    rows={mode:[] for mode in ['rk_coarse','rk_detailed']}
                    for i in range(a.warmups+a.repeats):
                        for mode in list(rows)[::(-1 if i%2 else 1)]:
                            w.process.stdin.write(f'{mode} {count}\n'); w.process.stdin.flush()
                            row=w.read(); row['unattributed_seconds']=validate(row)
                            if i>=a.warmups: rows[mode].append(row)
                    data['cases'].append({'input':str(path.resolve()),'input_sha256':digest(path),
                        'workers':count,'metadata':w.metadata,'profiles':rows})
                    for mode,values in rows.items():
                        print(path.name,count,mode,{k:round(1000*statistics.median(r[k] for r in values),3)
                            for k in values[0] if k.endswith('_seconds')},flush=True)
                        print('COUNTS',{k:statistics.median(r[k] for r in values) for k in values[0] if not k.endswith('_seconds')},flush=True)
            finally: w.close()
        data['completed_utc']=datetime.now(timezone.utc).isoformat()
        json.dump(data,output,indent=2); output.write('\n')


if __name__=='__main__': main()
