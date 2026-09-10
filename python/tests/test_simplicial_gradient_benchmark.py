"""Protocol, measurement-boundary and higher-dimensional generator checks."""
import importlib.util
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import benchmark_simplicial_gradients as benchmark
import pls_phase_profile as profile
import profile_reduction_kernel_simplicial as rk_profile
import validate_simplicial_gradient_ab as ab_validation


class SimplicialGradientBenchmarkTests(unittest.TestCase):
    def test_closure_profile_accounting(self):
        row = dict(builder_seconds=.1, kernel_seconds=1., algorithm_seconds=1.1,
                   setup_seconds=.1, level_wall_seconds=.6, replay_seconds=.2,
                   closure_seconds=.4, closure_initial_seconds=.1, closure_packed_seconds=.08,
                   closure_traversal_seconds=.1, closure_sort_seconds=.1,
                   closure_materialize_seconds=.05, closure_sparse_cells=2,
                   closure_sparse_entries=10, closure_boundary_visits=20,
                   closure_duplicate_faces=5, closure_index_growths=1, closure_entry_growths=1)
        self.assertAlmostEqual(rk_profile.validate(row), .1)
        indexed = dict(row, closure_boundary_visits=13, closure_boundary_index_seconds=.01,
                       closure_boundary_index_visits=10, closure_boundary_index_entries=5)
        self.assertAlmostEqual(rk_profile.validate(indexed), .1)
        for bad in [dict(indexed, closure_boundary_index_seconds=.03),
                    dict(indexed, closure_boundary_index_entries=11),
                    dict(indexed, closure_boundary_visits=14),
                    {k:v for k,v in indexed.items() if k != 'closure_boundary_index_visits'}]:
            with self.assertRaises(ValueError):
                rk_profile.validate(bad)
        for bad in [dict(row, closure_packed_seconds=.11), dict(row, closure_sort_seconds=.2),
                    dict(row, closure_boundary_visits=12),
                    dict(row, closure_sparse_cells=11), dict(row, closure_duplicate_faces=21),
                    dict(row, closure_index_growths=11), dict(row, closure_entry_growths=11),
                    {k:v for k,v in row.items() if k != 'closure_traversal_seconds'}]:
            with self.assertRaises(ValueError):
                rk_profile.validate(bad)

    def test_protected_scan_elision_validation(self):
        before = dict(local_candidate_visits=100, local_sparse_candidate_visits=80,
                      local_protected_candidate_visits=50, local_removed_candidate_visits=10,
                      local_sparse_scan_passes=5, local_membership_tests=20,
                      local_large_membership_tests=15, local_membership_comparisons=70,
                      local_large_membership_comparisons=60, local_coboundary_visits=30)
        after = dict(before, local_candidate_visits=60, local_sparse_candidate_visits=40,
                     local_protected_candidate_visits=10)
        ab_validation.validate_protected_scan_elision(before, before)
        ab_validation.validate_protected_scan_elision(before, after)
        # Reject changed queries/removals as well as a false total-visit saving.
        for key in after:
            with self.assertRaises(AssertionError, msg=key):
                ab_validation.validate_protected_scan_elision(before, dict(after, **{key: after[key] + 1}))
        with self.assertRaises(AssertionError):
            ab_validation.validate_protected_scan_elision(after, before)

    def test_rk_profile_accounting(self):
        row = dict(builder_seconds=.1, kernel_seconds=1., algorithm_seconds=1.1,
                   setup_seconds=.1, level_wall_seconds=.6, replay_seconds=.2)
        self.assertAlmostEqual(rk_profile.validate(row), .1)
        for bad in [dict(row, algorithm_seconds=2.), dict(row, setup_seconds=.9),
                    dict(row, kernel_seconds=float('nan')), dict(row, replay_seconds=-1.)]:
            with self.assertRaises(ValueError):
                rk_profile.validate(bad)

    def test_fine_profile_accounting(self):
        values = {key: .1 for key in profile.FIELDS}
        values["cleanup"] = .5
        parsed = profile.validate(values, 1., 1., .2, 3.)
        self.assertAlmostEqual(parsed["gradient_unattributed"], .3)
        self.assertAlmostEqual(parsed["setup_unattributed"], .4)
        self.assertAlmostEqual(parsed["local_unattributed"], .8)
        for bad in [dict(values, keys_cleanup=2.), dict(values, cleanup=.6),
                    dict(values, execution=float('nan')), dict(values, owner_keys=-1.),
                    {k:v for k,v in values.items() if k != 'schedule'}]:
            with self.assertRaises(ValueError):
                profile.validate(bad, 1., 1., .2, 3.)
        with self.assertRaises(ValueError):
            profile.validate(values, .1, 1., .2, 3.)
        with self.assertRaises(ValueError):
            profile.validate(values, 1., 1., .2, 2.)

    def test_rk_search_counters(self):
        row = dict(builder_seconds=.1, kernel_seconds=1., algorithm_seconds=1.1,
                   setup_seconds=.1, level_wall_seconds=.6, replay_seconds=.2,
                   local_membership_tests=10, local_candidate_visits=100,
                   local_membership_comparisons=40, local_large_membership_tests=5,
                   local_large_membership_comparisons=25, local_sparse_scan_passes=4,
                   local_sparse_candidate_visits=80, local_removed_candidate_visits=20,
                   local_protected_candidate_visits=30)
        self.assertAlmostEqual(rk_profile.validate(row), .1)
        for bad in [dict(row, local_large_membership_tests=11),
                    dict(row, local_membership_comparisons=9),
                    dict(row, local_large_membership_comparisons=41),
                    dict(row, local_large_membership_comparisons=4),
                    dict(row, local_sparse_candidate_visits=101),
                    dict(row, local_removed_candidate_visits=51),
                    {k:v for k,v in row.items() if k != 'local_sparse_scan_passes'}]:
            with self.assertRaises(ValueError):
                rk_profile.validate(bad)

    @unittest.skipUnless(os.environ.get("MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK"),
                         "native benchmark executable not supplied")
    def test_native_rk_search_profile(self):
        binary = Path(os.environ["MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK"]).resolve()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            path = directory / "input.txt"
            path.write_text(benchmark.grid_input(6, 2, 2))
            worker = benchmark.Worker(binary, path, directory / "dump")
            try:
                identities = []
                for count in [1, 4]:
                    worker.process.stdin.write(f"rk_detailed {count}\n")
                    worker.process.stdin.flush()
                    row = worker.read()
                    rk_profile.validate(row)
                    self.assertGreater(row['local_large_membership_tests'], 0)
                    self.assertGreater(row['local_removed_candidate_visits'], 0)
                    self.assertGreater(row['closure_sparse_cells'], 0)
                    self.assertGreater(row['closure_duplicate_faces'], 0)
                    self.assertGreater(row['closure_boundary_index_entries'], 0)
                    self.assertEqual(row['closure_boundary_visits'],
                                     row['closure_sparse_entries'] - row['closure_sparse_cells']
                                     + row['closure_duplicate_faces'])
                    identities.append(tuple(row[k] for k in rk_profile.SEARCH_FIELDS + rk_profile.CLOSURE_COUNTS[:4]))
                self.assertEqual(*identities)
            finally:
                worker.close()

    def test_grid_shape_and_injective_values(self):
        for dimension in range(1, 8):
            text = benchmark.grid_input(dimension, 2, 0)
            lines = text.splitlines()
            _, d, vertices, cells = lines[0].split()
            self.assertEqual(int(d), dimension)
            self.assertEqual(int(vertices), 2 ** dimension)
            self.assertEqual(int(cells), len(list(itertools.permutations(range(dimension)))))
            values = list(map(int, lines[1].split()))
            self.assertEqual(sorted(values), list(range(int(vertices))))
            self.assertEqual(len(lines) - 2, int(cells))
            self.assertEqual(len(set(lines[2:])), int(cells))
            used = set()
            for line in lines[2:]:
                cell = list(map(int, line.split()))
                self.assertEqual(len(set(cell)), dimension + 1)
                used.update(cell)
            self.assertEqual(used, set(range(int(vertices))))
            self.assertEqual(text, benchmark.grid_input(dimension, 2, 0))
        self.assertNotEqual(benchmark.grid_input(4, 2, 0), benchmark.grid_input(4, 2, 2))

    def test_grid_budget(self):
        for d, side in [(0, 2), (8, 2), (4, 1), (7, 10)]:
            with self.assertRaises(ValueError):
                benchmark.grid_input(d, side, 0)

    def test_phase_validation(self):
        run = {a: dict(builder_seconds=1., kernel_seconds=2., algorithm_seconds=3.)
               for a in benchmark.ALGORITHMS}
        benchmark.check_runs([run], 1)
        run["f_max"]["algorithm_seconds"] = 4.
        with self.assertRaises(AssertionError):
            benchmark.check_runs([run], 1)
        run["f_max"]["algorithm_seconds"] = float("nan")
        with self.assertRaises(AssertionError):
            benchmark.check_runs([run], 1)
        with self.assertRaises(AssertionError):
            benchmark.check_runs([], 1)

    def test_summary_direction_and_raw_preservation(self):
        runs = [{a: {p: float(i + 1) for p in benchmark.PHASES}
                 for i, a in enumerate(benchmark.ALGORITHMS)}] * 3
        samples = [{"baseline": runs, "candidate": runs}] * 4
        original = json.dumps(samples)
        summary = benchmark.summaries(samples)
        self.assertEqual(summary["process_lower_stars"]["algorithm_seconds"]["median_paired_ratio"], 1.)
        comparison = benchmark.comparisons(samples)
        self.assertEqual(comparison["reduction_kernel/process_lower_stars"]["median_paired_ratio"], 1.5)
        self.assertEqual(original, json.dumps(samples))

    @unittest.skipUnless(os.environ.get("MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK"),
                         "native benchmark executable not supplied")
    def test_native_worker_protocol(self):
        binary = Path(os.environ["MORSEFRAMES_SIMPLICIAL_GRADIENT_BENCHMARK"]).resolve()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            path = directory / "input.txt"
            path.write_text(benchmark.grid_input(4, 2, 0))
            worker = benchmark.Worker(binary, path, directory / "dump")
            try:
                identity = worker.metadata["identity"]
                self.assertEqual(identity["dimension"], 4)
                self.assertEqual(sum(identity["lower_star_sizes"]), identity["simplices"])
                self.assertEqual(set(identity["algorithms"]), set(benchmark.ALGORITHMS))
                for count in [1, 4]:
                    runs = worker.run(count, 6)
                    benchmark.check_runs(runs, 6)
                    self.assertEqual(len({tuple(r) for r in runs}), 6)
                    worker.process.stdin.write(f"profile {count}\n"); worker.process.stdin.flush()
                    row = worker.read()
                    fine = profile.validate(row["pls_profile_seconds"], row["setup_seconds"],
                        row["local_wall_seconds"], row["replay_seconds"],
                        row["algorithm_seconds"] - row["builder_seconds"])
                    self.assertGreater(fine["cleanup"], 0)
                    for mode in ["rk_coarse", "rk_detailed"]:
                        worker.process.stdin.write(f"{mode} {count}\n")
                        worker.process.stdin.flush()
                        row = worker.read()
                        self.assertGreaterEqual(rk_profile.validate(row), 0)
                        self.assertGreater(row["level_wall_seconds"], 0)
                        if mode == "rk_detailed":
                            self.assertGreater(row["rounds"], 0)
                            self.assertGreater(row["closure_seconds"], 0)
                        else:
                            self.assertEqual(row["rounds"], 0)
                            self.assertEqual(row["closure_seconds"], 0)
                            self.assertTrue(all(row[k] == 0 for k in rk_profile.CLOSURE_FIELDS))
                            self.assertTrue(all(row[k] == 0 for k in rk_profile.BOUNDARY_INDEX_FIELDS))
            finally:
                worker.close()
            self.assertEqual(worker.process.returncode, 0)
            for algorithm in range(3):
                memory = json.loads(subprocess.check_output(
                    [str(binary), str(path), "--memory", str(algorithm), "2"], text=True))
                self.assertGreaterEqual(memory["gradient_peak_bytes"], memory["complex_peak_bytes"])
                self.assertGreater(memory["steps"], 0)
            for bad in ["morseframes-ttk-v1 1 2 1\n0 0\n0 1\n",
                        "morseframes-ttk-v1 1 2 1\n0 1\n0\n",
                        "morseframes-ttk-v1 1 2 1\n0 1\n0 0\n"]:
                path.write_text(bad)
                process = subprocess.run([str(binary), str(path), str(directory / "bad")],
                                         capture_output=True, text=True)
                self.assertNotEqual(process.returncode, 0)


if __name__ == "__main__":
    unittest.main()
