import csv
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_tetrahedral_worker_scaling as bench  # noqa: E402
import render_tetrahedral_worker_scaling as render  # noqa: E402


class TetrahedralWorkerScalingBenchmarkTest(unittest.TestCase):
    def test_worker_scaling_contract(self):
        rows = bench.benchmark_scaling(
            seed=0,
            grid_size=3,
            workers=(1, 2),
            repeats=1,
            warmups=0,
        )

        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row.family == "injective-volume" for row in rows))
        self.assertTrue(all(row.sequence_matches_sequential for row in rows))
        self.assertTrue(all(row.barcode_matches_standard for row in rows))
        for strategy in bench.STRATEGIES:
            strategy_rows = [row for row in rows if row.strategy == strategy]
            self.assertEqual([row.max_workers for row in strategy_rows], [1, 2])
            self.assertEqual(strategy_rows[0].sequence_speedup_vs_one_worker, 1.0)
            self.assertEqual(strategy_rows[0].total_speedup_vs_one_worker, 1.0)
            self.assertEqual(
                strategy_rows[0].critical_simplices_by_dimension,
                strategy_rows[1].critical_simplices_by_dimension,
            )

        output = StringIO()
        bench.write_rows(rows, output, "csv")
        serialized = list(csv.DictReader(StringIO(output.getvalue())))
        with tempfile.TemporaryDirectory() as directory:
            table_path = Path(directory) / "scaling.tex"
            render.render_table(serialized, table_path)
            table = table_path.read_text()
        self.assertIn("ProcessLowerStars", table)
        self.assertIn("Reduction kernel", table)
        self.assertIn("Sequence speedup", table)


if __name__ == "__main__":
    unittest.main()
