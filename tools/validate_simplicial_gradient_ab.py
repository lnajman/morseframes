#!/usr/bin/env python3
"""Audit saved simplicial A/B evidence without running new timings."""
import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import statistics
import tempfile

import benchmark_simplicial_gradients as benchmark
from pls_phase_profile import validate as validate_pls
from profile_reduction_kernel_simplicial import validate as validate_rk
from profile_reduction_kernel_simplicial import SEARCH_FIELDS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equivalent(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            equivalent(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            equivalent(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-15), (a, b)
    else:
        assert a == b, (a, b)


def validate_protected_scan_elision(baseline, candidate):
    """Allow only skipped protected visits, not changed non-protected work.

    Pair this with the ordinary work invariants and complete ordered dumps.
    Preparation is included in elapsed time, but is not a repeated scan visit.
    """
    protected = 'local_protected_candidate_visits'
    assert candidate[protected] <= baseline[protected]
    saved = baseline[protected] - candidate[protected]
    for key in ('local_candidate_visits', 'local_sparse_candidate_visits'):
        assert baseline[key] - candidate[key] == saved, key
    for key in ('local_removed_candidate_visits', 'local_sparse_scan_passes',
                'local_membership_tests', 'local_large_membership_tests',
                'local_membership_comparisons', 'local_large_membership_comparisons',
                'local_coboundary_visits'):
        assert baseline[key] == candidate[key], key


def audit(path, allow_protected_scan_elision=False):
    data = json.loads(path.read_text())
    assert data['schema'] == 'simplicial-gradient-ab-v1' and data['completed']
    assert not data['source_status'] and not data['headers_patch']
    for field, source in {
        'driver_sha256': 'benchmarks/benchmark_simplicial_gradients.cpp',
        'runner_sha256': 'tools/benchmark_simplicial_gradients.py',
        'helper_sha256': 'tools/benchmark_reduction_kernel_ab.py',
        'profile_helper_sha256': 'benchmarks/pls_profile.hpp',
        'rk_profile_validator_sha256': 'tools/profile_reduction_kernel_simplicial.py',
    }.items():
        assert data[field] == digest(benchmark.ROOT / source), field
    with tempfile.TemporaryDirectory(prefix='simplicial-audit-') as directory:
        for version in benchmark.VERSIONS:
            build = data['builds'][version]
            target = Path(directory) / version
            revision = benchmark.snapshot_headers(build['revision'], target)
            assert revision == build['revision']
            assert benchmark.header_digest(target / 'include') == build['headers_sha256']
    settings = data['settings']
    coverage = {(i, w) for i in range(len(data['inputs'])) for w in settings['workers']}
    assert len(data['cases']) == len(coverage)
    assert {(c['input_index'], c['workers']) for c in data['cases']} == coverage
    for item in data['inputs']:
        assert item['sha256'] == digest(Path(item['path'])) and item['exact_dumps_match']
        identity = item['metadata']['candidate']['identity']
        assert identity == item['metadata']['baseline']['identity']
        assert sum(identity['lower_star_sizes']) == identity['simplices']
        for algorithm in identity['algorithms'].values():
            counts = algorithm['critical_counts']
            assert sum((-1)**d * n for d, n in enumerate(counts)) == identity['euler']
    static_counts = ('rounds', 'facet_kernels', 'facet_cell_visits',
                     'local_candidate_visits', 'local_coboundary_visits',
                     'local_membership_tests', 'inline_cell_overflows', 'inline_event_overflows')
    if settings['rk_profiles'] and all('local_sparse_scan_passes' in
            data['cases'][0]['rk_profiles'][v]['rk_detailed'][0] for v in benchmark.VERSIONS):
        static_counts += ('local_large_membership_tests', 'local_sparse_scan_passes',
                          'local_sparse_candidate_visits', 'local_removed_candidate_visits',
                          'local_protected_candidate_visits')
    if allow_protected_scan_elision:
        assert settings['rk_profiles'] > 0
        assert 'local_protected_candidate_visits' in static_counts
        static_counts = tuple(k for k in static_counts if k not in (
            'local_candidate_visits', 'local_sparse_candidate_visits',
            'local_protected_candidate_visits'))
    if settings['rk_profiles'] and all('closure_sparse_cells' in
            data['cases'][0]['rk_profiles'][v]['rk_detailed'][0] for v in benchmark.VERSIONS):
        static_counts += ('closure_sparse_cells', 'closure_sparse_entries', 'closure_duplicate_faces')
    for case in data['cases']:
        assert len(case['samples']) == settings['blocks']
        assert Counter(tuple(s['order']) for s in case['samples']) == {
            ('baseline', 'candidate'): settings['blocks'] // 2,
            ('candidate', 'baseline'): settings['blocks'] // 2}
        counts = set()
        for version in benchmark.VERSIONS:
            orders = Counter()
            for sample in case['samples']:
                benchmark.check_runs(sample[version], settings['repeats'])
                orders.update(tuple(r) for r in sample[version])
            assert orders == {p: settings['blocks'] * settings['repeats'] // 6
                              for p in itertools.permutations(benchmark.ALGORITHMS)}
            assert len(case['profiles'][version]) == settings['profiles']
            for p in case['profiles'][version]:
                equivalent(p['validated_fine_seconds'], validate_pls(
                    p['pls_profile_seconds'], p['setup_seconds'], p['local_wall_seconds'],
                    p['replay_seconds'], p['algorithm_seconds'] - p['builder_seconds']))
            for mode in ('rk_coarse', 'rk_detailed'):
                rows = case['rk_profiles'][version][mode]
                assert len(rows) == settings['rk_profiles']
                for p in rows:
                    equivalent(p['unattributed_seconds'], validate_rk(p))
                    if mode == 'rk_detailed':
                        counts.add(tuple(p[k] for k in static_counts))
            for rows in case['memory'][version].values():
                assert len(rows) == settings['memory_repeats']
        assert len(counts) <= 1, 'RK work changed despite expected exact sequence parity'
        if settings['rk_profiles'] and all('closure_boundary_index_visits' in
                case['rk_profiles'][v]['rk_detailed'][0] for v in benchmark.VERSIONS):
            reference = case['rk_profiles']['baseline']['rk_detailed'][0]
            for version in benchmark.VERSIONS:
                for row in case['rk_profiles'][version]['rk_detailed']:
                    # Boundary indexing may skip cross-level records, but
                    # cannot change any local reduction or same-level closure.
                    assert row['closure_boundary_visits'] <= reference['closure_boundary_visits']
                    for key in SEARCH_FIELDS:
                        assert row[key] == reference[key], key
        if allow_protected_scan_elision:
            reference = case['rk_profiles']['baseline']['rk_detailed'][0]
            for version in benchmark.VERSIONS:
                rows = case['rk_profiles'][version]['rk_detailed']
                # Check repeat stability, including counters deliberately
                # excluded from the cross-version equality check above.
                for row in rows:
                    for key in ('local_candidate_visits', 'local_sparse_candidate_visits',
                                'local_protected_candidate_visits'):
                        assert row[key] == rows[0][key]
                    validate_protected_scan_elision(reference, row)
        equivalent(case['summary'], benchmark.summaries(case['samples']))
        equivalent(case['comparison'], benchmark.comparisons(case['samples']))
    print('VALIDATED', path.name, len(data['cases']), 'configurations', digest(path))
    for case in data['cases']:
        name = Path(data['inputs'][case['input_index']]['path']).stem
        times = case['summary']['reduction_kernel']['algorithm_seconds']
        methods = {a: round(1000 * case['summary'][a]['algorithm_seconds']['median_seconds']['candidate'], 3)
                   for a in benchmark.ALGORITHMS}
        closures = {v: round(1000 * statistics.median(
            p['closure_seconds'] for p in case['rk_profiles'][v]['rk_detailed']), 3)
            for v in benchmark.VERSIONS} if settings['rk_profiles'] else {}
        print(name, case['workers'], 'RK ms',
              {v: round(1000*t, 3) for v, t in times['median_seconds'].items()},
              'ratio', round(times['median_paired_ratio'], 3),
              'CI', [round(t, 3) for t in times['paired_ratio_bootstrap_95_interval']],
              'methods', methods, 'closure', closures)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', type=Path, nargs='+')
    parser.add_argument('--reversed-confirmation', action='store_true')
    parser.add_argument('--allow-protected-scan-elision', action='store_true',
                        help='allow only protected-visit removal; require all other local work unchanged')
    args = parser.parse_args()
    datasets = [audit(path, args.allow_protected_scan_elision) for path in args.inputs]
    if args.reversed_confirmation:
        assert len(datasets) == 2
        a, b = datasets
        assert [i['path'] for i in a['inputs']] == list(reversed([i['path'] for i in b['inputs']]))
        assert a['settings']['workers'] == list(reversed(b['settings']['workers']))
        for version in benchmark.VERSIONS:
            assert a['builds'][version]['headers_sha256'] == b['builds'][version]['headers_sha256']
        print('Reversed confirmation schedules and frozen header identities verified')


if __name__ == '__main__':
    main()
