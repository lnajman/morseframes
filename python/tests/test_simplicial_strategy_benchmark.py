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

    def test_injective_volume_contract(self):
        grid_size = 3
        complex_ = bench.make_injective_volume(seed=0, grid_size=grid_size)
        vertex_values = [
            complex_.filtration(simplex)
            for simplex in range(complex_.size)
            if complex_.dimension(simplex) == 0
        ]
        tetrahedra = [
            simplex
            for simplex in range(complex_.size)
            if complex_.dimension(simplex) == 3
        ]

        self.assertEqual(len(vertex_values), grid_size**3)
        self.assertEqual(len(set(vertex_values)), len(vertex_values))
        self.assertEqual(complex_.num_levels, len(vertex_values))
        self.assertEqual(len(tetrahedra), 6 * (grid_size - 1) ** 3)
        for simplex in range(complex_.size):
            maximum = max(
                complex_.filtration(complex_.find_simplex((vertex,)))
                for vertex in complex_.vertices(simplex)
            )
            self.assertEqual(complex_.filtration(simplex), maximum)

    def test_tetrahedral_strategy_comparison(self):
        algorithms = (
            "process-lower-stars",
            "process-lower-stars-parallel",
            "flooding-reduction-kernel",
            "flooding-reduction-kernel-parallel",
            "f-max",
        )
        rows = bench.benchmark_volume(
            seed=0,
            grid_size=3,
            algorithms=algorithms,
            parallel_workers=2,
            repeats=1,
            warmups=0,
        )

        self.assertTrue(all(row.family == "injective-volume" for row in rows))
        self.assertTrue(all(row.max_dimension == 3 for row in rows))
        self.assertTrue(all(row.num_vertices == 27 for row in rows))
        self.assertTrue(all(row.barcode_matches_standard for row in rows))
        self.assertEqual(
            rows[0].critical_simplices_by_dimension,
            rows[1].critical_simplices_by_dimension,
        )
        self.assertEqual(
            rows[2].critical_simplices_by_dimension,
            rows[3].critical_simplices_by_dimension,
        )

    def test_strategy_comparison_tracks_critical_quality(self):
        algorithms = (
            "process-lower-stars",
            "process-lower-stars-parallel",
            "flooding-reduction-kernel",
            "flooding-reduction-kernel-parallel",
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
        kernel = rows[2]
        parallel_kernel = rows[3]
        self.assertTrue(process.sequence_matches_process_lower_stars)
        self.assertTrue(parallel.sequence_matches_process_lower_stars)
        self.assertEqual(
            process.critical_simplices_by_dimension,
            parallel.critical_simplices_by_dimension,
        )
        self.assertEqual(parallel.max_workers, 2)
        self.assertEqual(parallel_kernel.max_workers, 2)
        self.assertEqual(
            kernel.critical_simplices_by_dimension,
            parallel_kernel.critical_simplices_by_dimension,
        )
        if process.cpp_backend:
            self.assertTrue(kernel.reduction_kernel_metrics_available)
            self.assertTrue(parallel_kernel.reduction_kernel_metrics_available)
            self.assertGreater(kernel.reduction_kernel_rounds, 0)
            self.assertLessEqual(parallel_kernel.reduction_kernel_executor_workers, 2)
        else:
            self.assertFalse(kernel.reduction_kernel_metrics_available)
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
        self.assertIn("Reduction kernel", table)

        with tempfile.TemporaryDirectory() as directory:
            kernel_table_path = Path(directory) / "kernel.tex"
            render.render_kernel_table(serialized, kernel_table_path)
            kernel_table = kernel_table_path.read_text()
        self.assertIn("Facet batches", kernel_table)
        self.assertIn("Reduction kernel", kernel_table)


if __name__ == "__main__":
    unittest.main()
