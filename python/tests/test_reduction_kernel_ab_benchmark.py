import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
import benchmark_reduction_kernel_ab as bench


class ReductionKernelAbBenchmarkTest(unittest.TestCase):
    def test_legacy_defaults_and_explicit_plateau(self):
        with tempfile.TemporaryDirectory() as directory:
            common = ["--baseline", "HEAD", "--output", str(Path(directory) / "raw.json")]
            default = bench.parse_args(common)
            self.assertEqual((default.family, default.filtration), ("volume", "lower-star"))
            plateau = bench.parse_args(common + ["--family", "terrain", "--filtration", "plateau"])
            self.assertEqual((plateau.family, plateau.filtration), ("terrain", "plateau"))

    def test_paired_block_summary_retains_variability(self):
        samples = [{"baseline": [1.0, 2.0, 9.0], "candidate": [0.5, 1.0, 8.0]}] * 4
        summary = bench.summarize(samples)
        self.assertEqual(summary["median_paired_ratio"], 0.5)
        self.assertEqual(summary["median_seconds"], {"baseline": 2.0, "candidate": 1.0})
        self.assertEqual(summary["paired_ratio_bootstrap_95_interval"], [0.5, 0.5])
        self.assertEqual(summary["faster_blocks"], 4)


if __name__ == "__main__":
    unittest.main()
