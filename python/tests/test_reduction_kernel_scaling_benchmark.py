import csv
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_reduction_kernel_scaling as bench  # noqa: E402
import render_reduction_kernel_scaling as render  # noqa: E402


class ReductionKernelScalingBenchmarkTest(unittest.TestCase):
    def test_worker_scaling_preserves_sequence_and_barcode(self):
        rows = bench.benchmark_scaling(
            seed=0,
            grid_size=4,
            workers=(1, 2),
            repeats=1,
            warmups=0,
        )

        self.assertEqual([row.max_workers for row in rows], [1, 2])
        self.assertTrue(all(row.sequence_matches_sequential for row in rows))
        self.assertTrue(all(row.barcode_matches_standard for row in rows))
        self.assertEqual(rows[0].sequence_speedup_vs_one_worker, 1.0)
        self.assertEqual(rows[0].total_speedup_vs_one_worker, 1.0)
        self.assertGreater(rows[1].sequence_speedup_vs_one_worker, 0.0)
        self.assertGreater(rows[1].total_speedup_vs_one_worker, 0.0)
        self.assertEqual(
            rows[0].critical_simplices_by_dimension,
            rows[1].critical_simplices_by_dimension,
        )

        output = StringIO()
        bench.write_rows(rows, output, "csv")
        serialized = list(csv.DictReader(StringIO(output.getvalue())))
        with tempfile.TemporaryDirectory() as directory:
            table_path = Path(directory) / "scaling.tex"
            render.render_table(serialized, table_path)
            table = table_path.read_text()
        self.assertIn("Sequence speedup", table)
        self.assertIn("Level batches", table)


if __name__ == "__main__":
    unittest.main()
