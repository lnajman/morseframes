import csv
import math
import sys
import tempfile
import unittest
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_tetrahedral_phase_profile as bench  # noqa: E402
import render_tetrahedral_phase_profile as render  # noqa: E402


class TetrahedralPhaseProfileBenchmarkTest(unittest.TestCase):
    def test_zero_load_ratio_is_undefined(self):
        self.assertTrue(math.isnan(bench._ratio(10.0, 0.0)))
        self.assertTrue(math.isnan(bench._ratio(0.0, 0.0)))
        self.assertEqual(bench._ratio(10.0, 5.0), 2.0)

    @unittest.skipUnless(bench.mp.cpp_backend_available(), "requires native timings")
    def test_phase_profile_contract(self):
        rows = bench.benchmark_profile(
            seed=0,
            grid_size=4,
            workers=(1, 2),
            repeats=1,
        )

        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row.construction_seconds > 0.0 for row in rows))
        self.assertTrue(all(row.simplices_per_second > 0.0 for row in rows))
        self.assertTrue(all(row.matches_sequential for row in rows))
        self.assertTrue(all(row.critical_simplices_by_dimension for row in rows))
        process_rows = [row for row in rows if row.strategy == "process-lower-stars"]
        self.assertTrue(
            all(row.process_lower_stars_setup_seconds > 0.0 for row in process_rows)
        )
        self.assertTrue(
            all(
                row.process_lower_stars_local_wall_seconds > 0.0 for row in process_rows
            )
        )
        self.assertTrue(
            all(row.process_lower_stars_replay_seconds > 0.0 for row in process_rows)
        )
        self.assertTrue(
            all(
                row.process_lower_stars_task_time_imbalance >= 1.0
                for row in process_rows
            )
        )
        self.assertEqual(process_rows[0].process_lower_stars_setup_parallel_tasks, 0)
        self.assertGreater(process_rows[1].process_lower_stars_setup_parallel_tasks, 0)
        kernel_rows = [row for row in rows if row.strategy == "reduction-kernel"]
        self.assertTrue(
            all(math.isnan(row.process_lower_stars_setup_share) for row in kernel_rows)
        )
        self.assertTrue(
            all(row.reduction_kernel_level_wall_seconds > 0.0 for row in kernel_rows)
        )
        self.assertTrue(
            all(row.reduction_kernel_level_wall_share > 0.0 for row in kernel_rows)
        )
        self.assertTrue(all(row.reduction_kernel_rounds > 0 for row in kernel_rows))
        self.assertTrue(
            all(row.reduction_kernel_local_candidate_visits > 0 for row in kernel_rows)
        )
        self.assertTrue(
            all(row.reduction_kernel_local_coboundary_visits > 0 for row in kernel_rows)
        )
        self.assertTrue(
            all(row.reduction_kernel_local_coboundary_mask_tests > 0 for row in kernel_rows)
        )
        self.assertTrue(
            all(row.reduction_kernel_facet_discovery_mask_tests > 0 for row in kernel_rows)
        )
        self.assertTrue(
            all(row.reduction_kernel_inline_cell_overflows == 0 for row in kernel_rows)
        )
        self.assertTrue(
            all(row.reduction_kernel_inline_event_overflows == 0 for row in kernel_rows)
        )
        parallel_kernel = next(row for row in kernel_rows if row.max_workers == 2)
        self.assertGreater(parallel_kernel.reduction_kernel_level_chunks, 0)
        self.assertGreater(parallel_kernel.reduction_kernel_level_chunk_size, 0)
        self.assertGreater(parallel_kernel.reduction_kernel_task_parallelism, 0.0)
        self.assertGreaterEqual(
            parallel_kernel.reduction_kernel_task_time_imbalance, 1.0
        )
        # Dynamic claiming need not give every worker a chunk on this tiny
        # input. An idle worker has zero chunks/levels/simplices, so all three
        # max/min ratios are undefined; requiring >= 1 depends on scheduling.
        load_ratios = (
            parallel_kernel.reduction_kernel_chunk_imbalance,
            parallel_kernel.reduction_kernel_level_imbalance,
            parallel_kernel.reduction_kernel_simplex_imbalance,
        )
        if any(math.isnan(ratio) for ratio in load_ratios):
            self.assertTrue(all(math.isnan(ratio) for ratio in load_ratios))
        else:
            self.assertTrue(all(math.isfinite(ratio) and ratio >= 1.0
                                for ratio in load_ratios))

        output = StringIO()
        bench.write_rows(rows, output, "csv")
        serialized = list(csv.DictReader(StringIO(output.getvalue())))
        with tempfile.TemporaryDirectory() as directory:
            table_path = Path(directory) / "profile.tex"
            render.render_table(serialized, table_path)
            table = table_path.read_text()
        self.assertIn("Task parallelism", table)
        self.assertIn("ProcessLowerStars", table)
        self.assertIn("Reduction kernel", table)


if __name__ == "__main__":
    unittest.main()
