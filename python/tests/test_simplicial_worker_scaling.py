"""Crossover schedules and saved-evidence validation, without performance claims."""
from collections import Counter
import copy
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import benchmark_simplicial_worker_scaling as scaling


class WorkerScalingTests(unittest.TestCase):
    def test_schedule_balances_positions_and_pair_order(self):
        for workers in ([1, 2], [1, 2, 4, 8]):
            for blocks in (2 * len(workers), 4 * len(workers)):
                rows = scaling.schedule(workers, blocks)
                for position in range(len(workers)):
                    self.assertEqual(Counter(row[position] for row in rows),
                                     {w: blocks // len(workers) for w in workers})
                for a, b in itertools.combinations(workers, 2):
                    self.assertEqual(sum(row.index(a) < row.index(b) for row in rows), blocks // 2)
                self.assertEqual(scaling.schedule(workers, blocks, True), [r[::-1] for r in rows])

    def test_invalid_schedules(self):
        for workers, blocks in [([2, 4], 4), ([1, 1], 4), ([0, 1], 4),
                                ([1, 2, 4, 8], 6), ([1, 2], 0), ([1, 2], 6)]:
            with self.assertRaises(ValueError):
                scaling.schedule(workers, blocks)

    @staticmethod
    def blocks():
        blocks = []
        for i, order in enumerate(scaling.schedule([1, 2], 4)):
            runs = {}
            for w in order:
                runs[str(w)] = []
                for permutation in itertools.permutations(scaling.common.ALGORITHMS):
                    run = {}
                    for method in permutation:
                        factor = 1 if method == 'f_max' else 1 / w
                        total = (i + 1) * factor
                        run[method] = dict(builder_seconds=.1 * total,
                                           kernel_seconds=.9 * total, algorithm_seconds=total)
                    runs[str(w)].append(run)
            blocks.append(dict(order=order, runs=runs))
        return blocks

    def test_ratio_direction_and_sequential_control(self):
        summary = scaling.summarize_case(self.blocks(), [1, 2])
        rk = summary['2']['relative_to_one']['reduction_kernel']
        self.assertEqual(rk['median_paired_ratio'], .5)
        self.assertEqual(rk['paired_ratio_bootstrap_95_interval'], [.5, .5])
        self.assertEqual(summary['2']['relative_to_one']['f_max']['median_paired_ratio'], 1)
        self.assertEqual(summary['2']['comparison']['reduction_kernel/process_lower_stars']['median_paired_ratio'], 1)

    def test_block_ratio_not_ratio_of_global_medians(self):
        blocks = self.blocks()
        for b, base, candidate in zip(blocks, [1, 10, 100, 1000], [2, 5, 25, 125]):
            for w, total in [('1', base), ('2', candidate)]:
                for r in b['runs'][w]:
                    r['reduction_kernel']['algorithm_seconds'] = total
        row = scaling.summarize_case(blocks, [1, 2])['2']['relative_to_one']['reduction_kernel']
        self.assertEqual(row['median_paired_ratio'], .375)
        self.assertNotEqual(row['median_paired_ratio'], row['median_seconds']['candidate'] / row['median_seconds']['baseline'])

    def test_bad_runs_rejected(self):
        blocks = self.blocks()
        for b in blocks:
            for rows in b['runs'].values():
                scaling.common.check_runs(rows, 6)
        bad = copy.deepcopy(blocks[0]['runs']['1'])
        bad[0]['f_max']['algorithm_seconds'] = float('nan')
        with self.assertRaises(AssertionError):
            scaling.common.check_runs(bad, 6)


if __name__ == '__main__':
    unittest.main()
