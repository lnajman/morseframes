import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_ttk_reduction_kernel as bench


class TtkReductionKernelBenchmarkTest(unittest.TestCase):
    def test_ratios_use_paired_samples_and_keep_raw_variability(self):
        result = {
            "f_max_samples_seconds": [1, 10, 100],
            "reduction_kernel_samples_seconds": [50, 100, 2],
            "ttk_samples_seconds": [25, 50, 1],
            # Legacy best times must not silently replace the raw medians.
            "f_max_seconds": 1,
            "reduction_kernel_seconds": 2,
            "gradient_seconds": 1,
        }
        summary = bench.summarize_timings(result, 3)
        self.assertEqual(summary["f_max_seconds"], 10)
        self.assertEqual(summary["reduction_kernel_seconds"], 50)
        self.assertEqual(summary["reduction_kernel_ratio_vs_f_max"], 10)
        self.assertEqual(summary["reduction_kernel_ratio_vs_ttk"], 2)
        self.assertEqual(summary["f_max_iqr_seconds"], (5.5, 55))
        self.assertEqual(summary["ttk_samples_seconds"], (25, 50, 1))

    def test_incomplete_or_nonfinite_measurements_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "rebuild"):
            bench.summarize_timings({"f_max_seconds": 1}, 3)
        for values in ([], [1], [1, 2, float("nan")], [1, 2, 0]):
            with self.subTest(values=values):
                with self.assertRaises(ValueError):
                    bench.summarize_timings({"f_max_samples_seconds": values}, 3)


if __name__ == "__main__":
    unittest.main()
