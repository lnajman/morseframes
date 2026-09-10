import csv
import sys
import unittest
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_process_lower_stars as bench  # noqa: E402
import render_process_lower_stars_scaling as render  # noqa: E402


class ProcessLowerStarsBenchmarkTest(unittest.TestCase):
    def test_balanced_and_skewed_fan_workloads(self):
        balanced = bench.make_balanced_case(anchor_count=2, fan_size=2)
        skewed = bench.make_skewed_case(
            anchor_count=2,
            light_fan_size=1,
            heavy_fan_size=3,
        )

        self.assertEqual(sorted(balanced.lower_star_sizes), [1, 1, 1, 1, 2, 2, 2, 2, 7, 7])
        self.assertEqual(sum(balanced.lower_star_sizes), balanced.complex.size)
        self.assertEqual(max(skewed.lower_star_sizes), 10)
        self.assertGreater(
            max(skewed.lower_star_sizes) / min(skewed.lower_star_sizes),
            max(balanced.lower_star_sizes) / min(balanced.lower_star_sizes),
        )

    def test_lpt_load_estimate_is_deterministic(self):
        self.assertEqual(bench.estimated_task_loads((3, 2, 1), 2), (3, 3))
        self.assertEqual(bench.estimated_task_loads((3, 2, 1), 4), (3, 2, 1))

    def test_worker_sweep_preserves_sequence_and_critical_counts(self):
        case = bench.make_balanced_case(anchor_count=2, fan_size=1)
        rows = bench.benchmark_case(case, workers=(1, 2), repeats=1, warmups=0)

        self.assertEqual(len(rows), 3)
        self.assertEqual(
            {row.algorithm for row in rows},
            {"process-lower-stars", "process-lower-stars-parallel"},
        )
        self.assertTrue(all(row.exact_sequence_matches_sequential for row in rows))
        self.assertEqual(len({row.num_critical_simplices for row in rows}), 1)
        self.assertEqual(len({row.critical_simplices_by_dimension for row in rows}), 1)
        self.assertEqual(rows[0].sequence_speedup_vs_sequential, 1.0)
        self.assertEqual(rows[0].total_speedup_vs_sequential, 1.0)
        for row in rows:
            self.assertEqual(
                row.num_critical_simplices + 2 * row.num_regular_pairs,
                row.num_simplices,
            )
            self.assertGreaterEqual(row.sequence_seconds, 0.0)
            self.assertGreaterEqual(row.persistence_seconds, 0.0)
            self.assertAlmostEqual(
                row.total_seconds,
                row.sequence_seconds + row.persistence_seconds,
            )
            self.assertGreater(row.sequence_speedup_vs_sequential, 0.0)
            self.assertGreater(row.total_speedup_vs_sequential, 0.0)
            self.assertGreater(row.sequence_parallel_efficiency, 0.0)

        output = StringIO()
        bench.write_rows(rows, output, "csv")
        serialized = list(csv.DictReader(StringIO(output.getvalue())))
        self.assertEqual(len(serialized), len(rows))
        self.assertIn("critical_simplices_by_dimension", serialized[0])
        self.assertIn("sequence_seconds_per_eliminated_simplex", serialized[0])
        self.assertIn("sequence_speedup_vs_sequential", serialized[0])

        table_path = ROOT / "build" / "test_process_lower_stars_scaling_table.tex"
        table_path.parent.mkdir(parents=True, exist_ok=True)
        render.render_table(serialized, table_path)
        table = table_path.read_text()
        self.assertIn("Critical (by dim.)", table)
        self.assertIn("Balanced", table)

    def test_multi_scale_cli_computes_matched_heavy_fans(self):
        args = bench.parse_args(
            [
                "--anchors",
                "4",
                "--balanced-fans",
                "2",
                "5",
                "--light-fan",
                "1",
            ]
        )
        self.assertEqual(args.balanced_fans, [2, 5])
        self.assertIsNone(args.heavy_fan)


if __name__ == "__main__":
    unittest.main()
