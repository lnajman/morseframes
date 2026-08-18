import csv
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_simplicial_strategies as bench  # noqa: E402
import render_simplicial_strategy_benchmark as render  # noqa: E402


class SimplicialStrategyBenchmarkTest(unittest.TestCase):
    def test_injective_terrain_contract(self):
        complex_ = bench.make_injective_terrain(seed=0, grid_size=4)
        vertex_values = [
            complex_.filtration(simplex)
            for simplex in range(complex_.size)
            if complex_.dimension(simplex) == 0
        ]

        self.assertEqual(len(vertex_values), 16)
        self.assertEqual(len(set(vertex_values)), len(vertex_values))
        self.assertEqual(complex_.num_levels, len(vertex_values))
        for simplex in range(complex_.size):
            maximum = max(
                complex_.filtration(complex_.find_simplex((vertex,)))
                for vertex in complex_.vertices(simplex)
            )
            self.assertEqual(complex_.filtration(simplex), maximum)

    def test_strategy_comparison_tracks_critical_quality(self):
        algorithms = (
            "process-lower-stars",
            "process-lower-stars-parallel",
            "f-max",
            "f-min",
        )
        rows = bench.benchmark_terrain(
            seed=0,
            grid_size=4,
            algorithms=algorithms,
            parallel_workers=2,
            repeats=1,
            warmups=0,
        )

        self.assertEqual([row.algorithm for row in rows], list(algorithms))
        self.assertTrue(all(row.barcode_matches_standard for row in rows))
        process = rows[0]
        parallel = rows[1]
        self.assertTrue(process.sequence_matches_process_lower_stars)
        self.assertTrue(parallel.sequence_matches_process_lower_stars)
        self.assertEqual(
            process.critical_simplices_by_dimension,
            parallel.critical_simplices_by_dimension,
        )
        for row in rows:
            self.assertEqual(
                sum(row.critical_simplices_by_dimension),
                row.num_critical_simplices,
            )
            self.assertEqual(
                row.num_critical_simplices + 2 * row.num_regular_pairs,
                row.num_simplices,
            )
            self.assertGreater(row.sequence_speedup_vs_f_max, 0.0)
            self.assertGreater(row.total_speedup_vs_f_max, 0.0)

        output = StringIO()
        bench.write_rows(rows, output, "csv")
        serialized = list(csv.DictReader(StringIO(output.getvalue())))
        self.assertIn("critical_count_delta_vs_f_max", serialized[0])
        self.assertIn("sequence_speedup_vs_f_max", serialized[0])

        with tempfile.TemporaryDirectory() as directory:
            table_path = Path(directory) / "comparison.tex"
            render.render_table(serialized, table_path)
            table = table_path.read_text()
        self.assertIn("Critical/F-Max", table)
        self.assertIn("ProcessLowerStars", table)


if __name__ == "__main__":
    unittest.main()
