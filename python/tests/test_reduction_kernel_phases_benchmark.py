import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
import benchmark_reduction_kernel_phases as bench
import render_reduction_kernel_phases as render


class ReductionKernelPhasesBenchmarkTest(unittest.TestCase):
    def samples(self):
        row = dict(total_seconds=10.0, builder_seconds=1.0, setup_seconds=2.0,
                   level_wall_seconds=4.0, replay_seconds=1.0,
                   cumulative_level_task_seconds=12.0, max_parallel_levels=4,
                   local_reduction_seconds=9.0)
        return {"run": [[5.0, 8.0, 20.0]], "coarse": [[row]], "detailed": [[row]]}

    def test_medians_shares_and_nested_times(self):
        samples = self.samples()
        original = copy.deepcopy(samples)
        summary = bench.summarize(samples)
        self.assertEqual(summary["performance_seconds"]["median"], 8)
        coarse = summary["coarse"]
        self.assertEqual(coarse["shares"]["level_wall_seconds"]["median"], 0.4)
        self.assertAlmostEqual(coarse["other_share"]["median"], 0.2)
        self.assertEqual(coarse["effective_level_task_parallelism"]["median"], 3)
        self.assertEqual(coarse["overhead_ratio"], 1.25)
        self.assertEqual(samples, original)

    def test_invalid_timings_and_double_counting_rejected(self):
        for values in ([], [-1], [float("nan")], [float("inf")]):
            with self.assertRaises(ValueError):
                bench.distribution(values)
        samples = self.samples()
        samples["coarse"][0][0]["replay_seconds"] = 6.0
        with self.assertRaisesRegex(ValueError, "exceed"):
            bench.summarize(samples)

    def test_balanced_order_and_absent_level_parallelism(self):
        orders = [bench.worker_order([1, 2, 4, 8], b) for b in range(8)]
        for position in range(4):
            self.assertEqual(sorted(order[position] for order in orders),
                             [1, 1, 2, 2, 4, 4, 8, 8])
        samples = self.samples()
        samples["coarse"][0][0]["max_parallel_levels"] = 0
        self.assertIsNone(bench.summarize(samples)["coarse"]["effective_level_task_parallelism"])

    def test_renderer_requires_complete_results_and_baselines(self):
        with self.assertRaisesRegex(ValueError, "incomplete"):
            render.largest_case_rows({"cases": []})
        summary = bench.summarize(self.samples())
        summary["speedup_vs_one"] = 2
        case = dict(all_sequences_match=True, mode="lower-star", family="volume",
                    size=8, num_simplices=100,
                    measurements=[{"workers": 1, "summary": summary},
                                  {"workers": 8, "summary": summary}])
        data = {"completed_utc": "test", "cases": [case]}
        self.assertIn("Lower-star 3D", render.render_table(data))
        self.assertEqual(render.largest_case_rows(data)[0]["speedup"], 2)
        case["measurements"].pop()
        with self.assertRaisesRegex(ValueError, "one- and eight"):
            render.largest_case_rows(data)


if __name__ == "__main__":
    unittest.main()
