import csv
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_gradient_strategies as bench  # noqa: E402
import render_gradient_strategy_benchmark as render  # noqa: E402


class GradientStrategyBenchmarkTest(unittest.TestCase):
    def test_2d_and_3d_gradient_contract(self):
        options = {"workers": (1, 2), "repeats": 1, "warmups": 0}
        rows = bench.benchmark_terrain(0, 3, **options)
        rows += bench.benchmark_volume(0, 3, **options)

        self.assertEqual(len(rows), 14)
        self.assertEqual({row.dimension for row in rows}, {2, 3})
        self.assertTrue(all(row.construction_seconds > 0.0 for row in rows))
        self.assertTrue(all(row.simplices_per_second > 0.0 for row in rows))
        self.assertTrue(all(row.matches_sequential for row in rows))
        self.assertTrue(
            all(
                sum(row.critical_simplices_by_dimension)
                == row.num_critical_simplices
                for row in rows
            )
        )
        f_max_rows = [row for row in rows if row.algorithm == bench.mp.F_MAX_SEQUENCE]
        self.assertTrue(
            all(row.construction_time_ratio_vs_f_max == 1.0 for row in f_max_rows)
        )
        for family in {row.family for row in rows}:
            family_rows = [row for row in rows if row.family == family]
            f_max_count = next(
                row.num_critical_simplices
                for row in family_rows
                if row.algorithm == bench.mp.F_MAX_SEQUENCE
            )
            self.assertTrue(
                all(row.num_critical_simplices == f_max_count for row in family_rows)
            )

        output = StringIO()
        bench.write_rows(rows, output, "csv")
        serialized = list(csv.DictReader(StringIO(output.getvalue())))
        with tempfile.TemporaryDirectory() as directory:
            table_path = Path(directory) / "gradient.tex"
            render.render_table(serialized, table_path)
            table = table_path.read_text()
        self.assertIn("2D triangulated terrain", table)
        self.assertIn("3D tetrahedral volume", table)
        self.assertIn("Critical/F-Max", table)


if __name__ == "__main__":
    unittest.main()
