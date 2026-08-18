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
    @unittest.skipUnless(bench.mp.cpp_backend_available(), "requires native timings")
    def test_phase_profile_contract(self):
        rows = bench.benchmark_profile(
            seed=0,
            grid_size=4,
            workers=(1, 2),
            repeats=1,
        )

        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row.sequence_total_seconds > 0.0 for row in rows))
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
