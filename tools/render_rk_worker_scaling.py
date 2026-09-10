#!/usr/bin/env python3
"""Render exact lookup tables for the frozen worker-scaling study.

Audit raw evidence first. This renderer preserves seed/session granularity;
across-seed medians/ranges are descriptive, not confidence intervals.
"""
import argparse
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
median = statistics.median
METHODS = ('reduction_kernel', 'process_lower_stars', 'f_max')
FILES = ('rk-worker-scaling-main.json', 'rk-worker-scaling-confirmation.json',
         'rk-worker-scaling-controls-main.json', 'rk-worker-scaling-controls-confirmation.json')


def ratio(row):
    value = row['median_paired_ratio']
    lo, hi = row['paired_ratio_bootstrap_95_interval']
    return f'{value:.3f} [{lo:.3f}, {hi:.3f}]'


def speed(row):
    value = row['median_paired_ratio']
    lo, hi = row['paired_ratio_bootstrap_95_interval']
    return f'{1/value:.2f} [{1/hi:.2f}, {1/lo:.2f}]'


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |',
                      *['| ' + ' | '.join(map(str, row)) + ' |' for row in rows]])


def name(case):
    return Path(case['path']).stem


def short(case):
    return name(case).replace('grid-d', '').replace('-n2-', 'D/').replace('-n3-', 'D/').replace('-n4-', 'D/').replace('seed', 's')


def render(data):
    out = ['# Frozen worker scaling: numerical appendix',
           'Generated from the four raw studies linked in the [main report](rk_worker_scaling_benchmark.md). '
           'Times are milliseconds. Intervals are paired-block bootstrap 95% intervals within a session; '
           'cross-seed ranges are not intervals. All outliers and exceptions remain included. '
           'F-Max is sequential at every worker setting.']
    out += ['## Per-seed total-gradient times and comparisons',
            'Each row uses 48 unprofiled runs per method. Times are marginal medians; ratios are medians of '
            'paired-block ratios, not quotients of the displayed medians. RK speedup compares its one-worker '
            'time to its time at the listed count. Values below one for RK/PLS favor RK.']
    for label, study in zip(('Main', 'Confirmation'), data[:2]):
        rows = []
        for case in sorted(study['cases'], key=name):
            for w in ('1', '2', '4', '8'):
                s = case['summary'][w]
                rows.append([short(case), w, *[f"{1000*s['methods'][m]['algorithm_seconds']:.3f}" for m in METHODS],
                    speed(s['relative_to_one']['reduction_kernel']),
                    ratio(s['comparison']['reduction_kernel/process_lower_stars']),
                    ratio(s['comparison']['reduction_kernel/f_max'])])
        out += [f'### {label}', table(['Input', 'Workers', 'RK', 'PLS', 'F-Max',
                'RK speedup [interval]', 'RK/PLS [interval]', 'RK/F-Max [interval]'], rows)]

    out += ['## Critical simplices by dimension',
            'These ordered vectors are identical across sessions and worker counts within each method. '
            'The original 4D/seed0 method difference is retained. No cross-method equality is imposed.']
    rows = []
    for c in sorted(data[0]['cases'], key=name):
        i = c['metadata']['identity']
        rows.append([short(c), i['simplices'], *[str(i['algorithms'][m]['critical_counts']) for m in METHODS]])
    out.append(table(['Input', 'Simplices', 'RK', 'PLS', 'F-Max'], rows))

    out += ['## Loading and construction, kept separate',
            'One startup observation per input/study, not a dedicated construction experiment. '
            'These values are excluded from all gradient comparisons above.']
    other = {name(c): c for c in data[1]['cases']}
    rows = [[short(c), *[f"{1000*x['metadata'][key]:.3f}" for x in (c, other[name(c)])
                       for key in ('loading_seconds', 'construction_seconds')]]
            for c in sorted(data[0]['cases'], key=name)]
    out.append(table(['Input', 'Main load', 'Main construct', 'Confirm load', 'Confirm construct'], rows))

    out += ['## RK worker-task diagnostics',
            'Three separate coarse profiles per count. Durations are elapsed persistent-task lifetimes, '
            'not CPU time or per-level timings. The balance indicator is the median of '
            '`sum(task lifetimes)/(workers × longest task lifetime)` in each call. Level/simplex extrema '
            'need not belong to the same task. For one worker these parallel-only fields are not populated.']
    for label, study in zip(('Main', 'Confirmation'), data[:2]):
        rows = []
        for c in sorted(study['cases'], key=name):
            if not any(f'grid-d{d}-' in c['path'] for d in (5, 6, 7)):
                continue
            for w in ('2', '4', '8'):
                p = c['profiles'][w]['rk_coarse']
                mid = lambda k: median(r[k] for r in p)
                balance = median(r['cumulative_level_task_seconds']/(int(w)*r['max_level_task_seconds']) for r in p)
                rows.append([short(c), w, int(mid('level_chunk_size')),
                    f"{1000*mid('level_wall_seconds'):.3f}",
                    f"{1000*mid('min_level_task_seconds'):.3f}–{1000*mid('max_level_task_seconds'):.3f}",
                    f'{balance:.3f}',
                    f"{mid('min_worker_levels'):g}–{mid('max_worker_levels'):g}",
                    f"{mid('min_worker_simplices'):g}–{mid('max_worker_simplices'):g}"])
        out += [f'### {label} task diagnostics', table(['Input', 'Workers', 'Chunk', 'Level wall',
                    'Task min–max', 'Balance', 'Levels min–max', 'Simplices min–max'], rows)]

    out += ['## Coarse phases for every method',
            'RK and PLS phase rows below use separate instrumented calls (three per count), '
            'not the unprofiled headline samples. Each raw row passes phase accounting; marginal phase '
            'medians need not sum to the median total. F-Max exposes builder/kernel/total only. '
            'Full nested RK and PLS detail remains in the raw studies.']
    for label, study in zip(('Main', 'Confirmation'), data[:2]):
        for method in ('RK', 'PLS', 'F-Max'):
            rows = []
            for c in sorted(study['cases'], key=name):
                for w in ('1', '8'):
                    if method == 'RK':
                        p = c['profiles'][w]['rk_coarse']
                        keys = ['builder_seconds', 'setup_seconds', 'level_wall_seconds', 'replay_seconds',
                                'unattributed_seconds', 'algorithm_seconds']
                        vals = [median(r[k] for r in p) for k in keys]
                    elif method == 'PLS':
                        p = c['profiles'][w]['profile']
                        vals = [median(r[k] for r in p) for k in
                                ('builder_seconds', 'setup_seconds', 'local_wall_seconds', 'replay_seconds')]
                        vals += [median(r['algorithm_seconds'] - sum(r[k] for k in
                            ('builder_seconds', 'setup_seconds', 'local_wall_seconds', 'replay_seconds')) for r in p),
                            median(r['algorithm_seconds'] for r in p)]
                    else:
                        vals = [c['summary'][w]['methods']['f_max'][k] for k in common_phases()]
                    rows.append([short(c), w, *[f'{1000*v:.3f}' for v in vals]])
            headers = (['Builder', 'Kernel', 'Total'] if method == 'F-Max' else
                       ['Builder', 'Setup', 'Level wall' if method == 'RK' else 'Local wall', 'Replay',
                        'Remainder/cleanup', 'Total'])
            out += [f'### {label}: {method}', table(['Input', 'Workers', *headers], rows)]

    out += ['## Repeated boundary-index regression controls',
            'Current headers versus the pre-index baseline. Six paired blocks × two repetitions, '
            'separate from the scaling sweep. Ratios above one mean slower current RK; all controls are retained. '
            'PLS and F-Max code are unchanged and serve as context, not a correction factor.']
    for label, study in zip(('Main', 'Confirmation'), data[2:]):
        rows = []
        for c in study['cases']:
            nm = Path(study['inputs'][c['input_index']]['path']).stem
            rows.append([nm, c['workers'], *[ratio(c['summary'][m]['algorithm_seconds']) for m in METHODS]])
        out += [f'### {label} controls', table(['Input', 'Workers', 'RK after/before', 'PLS after/before',
                                                'F-Max after/before'], rows)]
    return '\n\n'.join(out) + '\n'


def common_phases():
    return ('builder_seconds', 'kernel_seconds', 'algorithm_seconds')


def figure(data, output, png=None):
    # Chart contract: scientific file figure for the existing Sphinx report.
    # Question: how do RK and PLS scale, and how much does the seed matter?
    # Faceted dot-and-range: 3 dimensions x 4 counts x 2 methods x 2 sessions,
    # each mark backed by six seed-specific paired speedups (48 runs/method).
    # Show medians and full seed ranges, not pooled confidence intervals.
    # Hard two-root palette: blue RK / orange PLS, circles / triangles;
    # main filled / confirmation open. No lines interpolating worker counts.
    # 8 x 9.5 inch vertically faceted SVG: labels stay legible in Sphinx's
    # narrow article column. QA in Sphinx and a PNG at output size.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams['svg.hashsalt'] = 'rk-worker-scaling-frozen'
    fig, axes = plt.subplots(3, 1, figsize=(8, 9.5), sharex=True, sharey=True)
    for axis, dimension in zip(axes, (5, 6, 7)):
        for study_index, study in enumerate(data[:2]):
            cases = [c for c in study['cases'] if f'grid-d{dimension}-' in c['path']]
            assert len(cases) == 6
            for method_index, (method, color, marker) in enumerate([
                    ('reduction_kernel', '#174A7E', 'o'),
                    ('process_lower_stars', '#C76822', '^')]):
                for x, w in enumerate(('1', '2', '4', '8')):
                    values = [1 / c['summary'][w]['relative_to_one'][method]['median_paired_ratio'] for c in cases]
                    mid = median(values)
                    axis.errorbar(x + (method_index - .5)*.22 + (study_index - .5)*.09, mid,
                        yerr=[[mid-min(values)], [max(values)-mid]], fmt=marker, color=color,
                        markerfacecolor=color if study_index == 0 else 'white',
                        markersize=4.5, capsize=2, elinewidth=1, markeredgewidth=1)
        axis.set_title(f'{dimension}D', loc='left', fontsize=12)
        axis.set_xticks(range(4), ['1', '2', '4', '8'])
        axis.set_xlim(-.4, 3.4)
        axis.set_ylim(0, 6)
        axis.axhline(1, color='#666666', linestyle='--', linewidth=.8)
        axis.grid(axis='y', color='#dddddd', linewidth=.6)
        axis.spines[['top', 'right']].set_visible(False)
    axes[1].set_ylabel('Speedup relative to one worker')
    axes[2].set_xlabel('Workers')
    fig.suptitle('Gradient speedup by worker count', x=.065, y=.98, ha='left', fontsize=15)
    fig.text(.065, .94, 'Six seeds per dimension/session: median and full range, not confidence intervals', fontsize=9)
    handles = [Line2D([], [], marker='o', color='#174A7E', linestyle='', label='RK'),
               Line2D([], [], marker='^', color='#C76822', linestyle='', label='PLS'),
               Line2D([], [], marker='o', color='#444444', linestyle='', label='Main: filled'),
               Line2D([], [], marker='o', color='#444444', markerfacecolor='white', linestyle='', label='Confirmation: open')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.06,.925), ncol=4, frameon=False, fontsize=9)
    fig.text(.065, .018, 'Frozen algorithms · 10 September 2026 · F-Max remains sequential\nGradient preparation included; loading and construction separate', fontsize=8)
    fig.subplots_adjust(left=.11, right=.985, bottom=.10, top=.86, hspace=.24)
    fig.savefig(output, metadata={'Date': None, 'Title': 'Gradient speedup by worker count'})
    if png:
        fig.savefig(png, dpi=140)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir', type=Path, default=ROOT.parent)
    p.add_argument('--output', type=Path, default=ROOT / 'docs/rk_worker_scaling_tables.md')
    p.add_argument('--check', action='store_true')
    p.add_argument('--figure', type=Path)
    p.add_argument('--png', type=Path)
    a = p.parse_args()
    data = [json.loads((a.data_dir / f).read_text()) for f in FILES]
    assert all(d['completed'] for d in data)
    result = render(data)
    if a.check:
        assert a.output.read_text() == result, 'Stale scaling tables'
        print('VALIDATED exact generated scaling tables')
    else:
        a.output.write_text(result)
    if a.figure:
        figure(data, a.figure, a.png)


if __name__ == '__main__':
    main()
