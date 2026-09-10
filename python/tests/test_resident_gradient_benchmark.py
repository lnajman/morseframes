import copy
from collections import Counter
from contextlib import redirect_stderr
from io import StringIO
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import benchmark_resident_gradients as bench
import render_resident_gradients as render


class ResidentGradientBenchmarkTest(unittest.TestCase):
    def raw(self):
        algorithms = {}
        for algorithm in bench.ALGORITHMS:
            outer = {k: 1.0 for k in bench.OUTER_PHASES[algorithm]}
            outer["gradient"] = 10.0
            inner = {k: 0.5 for k in bench.DETAIL_PHASES[algorithm]}
            algorithms[algorithm] = dict(
                critical_counts=[1, 0, 0], performance_seconds=[1.0, 10.0, 100.0],
                diagnostics=[dict(total_seconds=sum(outer.values()),
                                  phases_seconds=outer, gradient_details_seconds=inner)])
        algorithms["reduction_kernel"]["performance_seconds"] = [50.0, 100.0, 2.0]
        algorithms["ttk"]["performance_seconds"] = [25.0, 50.0, 1.0]
        return dict(schema="resident-gradient-v1", ttk_revision=bench.TTK_REVISION,
                    dimension=2, vertices=3, simplices=7, workers=1,
                    exact_reference_checks=True, critical_counts_match=True,
                    algorithms=algorithms)

    def test_totals_use_only_uninstrumented_paired_samples(self):
        raw = self.raw()
        original = copy.deepcopy(raw)
        s = bench.summarize(raw, 3, 1)
        self.assertEqual(s["f_max"]["total_seconds"]["median"], 10)
        self.assertEqual(s["reduction_kernel"]["total_seconds"]["median"], 50)
        self.assertEqual(s["paired_ratios"]["reduction_kernel/f_max"]["median"], 10)
        self.assertEqual(s["paired_ratios"]["reduction_kernel/ttk"]["median"], 2)
        self.assertEqual(raw, original)

    def test_nested_phases_are_separate_and_residual_is_explicit(self):
        s = bench.summarize(self.raw(), 3, 1)
        rk = s["reduction_kernel"]
        self.assertEqual(rk["diagnostic_total_seconds"]["median"], 12)
        self.assertAlmostEqual(rk["phase_shares"]["gradient"]["median"], 10 / 12)
        self.assertEqual(rk["gradient_details_seconds"]["unattributed"]["median"], 8.5)
        self.assertIsNone(s["ttk"]["gradient_details_seconds"])

    def test_incomplete_nonfinite_or_overlapping_phases_rejected(self):
        for mutate in (
            lambda r: r.update(schema="legacy"),
            lambda r: r.update(ttk_revision="unknown"),
            lambda r: r.update(exact_reference_checks=False),
            lambda r: r["algorithms"]["ttk"]["performance_seconds"].pop(),
            lambda r: r["algorithms"]["ttk"]["performance_seconds"].__setitem__(0, math.nan),
            lambda r: r["algorithms"]["ttk"]["diagnostics"][0]["phases_seconds"].pop("vertex_order"),
            lambda r: r["algorithms"]["ttk"]["diagnostics"][0].update(total_seconds=999.0),
            lambda r: r["algorithms"]["reduction_kernel"]["diagnostics"][0]["gradient_details_seconds"].update(replay=99.0),
        ):
            raw = self.raw()
            mutate(raw)
            with self.assertRaises(ValueError):
                bench.summarize(raw, 3, 1)

    def test_different_critical_counts_are_reported_not_hidden(self):
        raw = self.raw()
        raw["algorithms"]["ttk"]["critical_counts"] = [2, 1, 0]
        raw["critical_counts_match"] = False
        bench.summarize(raw, 3, 1)
        raw["critical_counts_match"] = True
        with self.assertRaisesRegex(ValueError, "agreement"):
            bench.summarize(raw, 3, 1)

    def test_orders_and_arguments(self):
        self.assertEqual(len(set(bench.ORDERS)), 6)
        for position in range(3):
            self.assertEqual(sorted(row[position] for row in bench.ORDERS), [0, 0, 1, 1, 2, 2])
        with tempfile.TemporaryDirectory() as directory:
            args = ["--benchmark", "native", "--input-dir", directory,
                    "--output", str(Path(directory) / "raw.json")]
            parsed = bench.parse_args(args)
            self.assertEqual(parsed.repeats % 6, 0)
            self.assertEqual(parsed.diagnostics % 6, 0)
            with redirect_stderr(StringIO()), self.assertRaises(SystemExit):
                bench.parse_args(args + ["--repeats", "5"])

    def test_renderer_recomputes_and_rejects_incomplete_studies(self):
        data = dict(schema="resident-gradient-study-v1", completed_utc="test",
                    arguments=dict(repeats=3, diagnostics=1, workers=[1],
                                   terrain_sizes=[3], volume_sizes=[], seeds=[0]),
                    cases=[dict(family="terrain", size=3, seed=0,
                                measurements=[dict(workers=1, raw=self.raw(), summary={})])])
        totals, phases = render.render_tables(data)
        self.assertIn("50000.000", totals)  # 50 seconds, displayed in milliseconds
        self.assertIn("2.000", totals)  # paired RK/TTK ratio, recomputed from raw data
        self.assertIn("vertex order", phases)
        self.assertIn("Total (diagnostic)", phases)
        del data["completed_utc"]
        with self.assertRaisesRegex(ValueError, "completed"):
            render.render_tables(data)

    def split_raw(self):
        raw = self.raw()
        raw.update(schema="resident-gradient-v2", input_loading_seconds=0.1)
        for algorithm in bench.ALGORITHMS:
            construction = [10.0, 100.0, 1.0]
            algorithm_times = ([1.0, 10.0, 100.0] if algorithm == "f_max" else
                               [2.0, 2.0, 100.0] if algorithm == "reduction_kernel" else
                               [2.0, 2.0, 51.0])
            rows = []
            for c, a in zip(construction, algorithm_times):
                row = {k: c / len(bench.CONSTRUCTION_PHASES[algorithm])
                       for k in bench.CONSTRUCTION_PHASES[algorithm]}
                row.update({k: a / len(bench.ALGORITHM_PHASES[algorithm])
                            for k in bench.ALGORITHM_PHASES[algorithm]})
                rows.append(row)
            raw["algorithms"][algorithm]["performance_phases_seconds"] = rows
            raw["algorithms"][algorithm]["performance_seconds"] = [sum(r.values()) for r in rows]
        return raw

    def test_split_comparison_uses_raw_samples_not_diagnostic_or_median_subtraction(self):
        raw = self.split_raw()
        original = copy.deepcopy(raw)
        s = bench.summarize(raw, 3, 1)
        self.assertEqual(s["f_max"]["algorithm_seconds"]["median"], 10)
        self.assertEqual(s["f_max"]["total_seconds"]["median"], 101)
        self.assertEqual(s["f_max"]["construction_seconds"]["median"], 10)
        self.assertNotEqual(101 - 10, s["f_max"]["algorithm_seconds"]["median"])
        self.assertEqual(s["algorithm_paired_ratios"]["reduction_kernel/ttk"]["median"], 1)
        self.assertEqual(bench.ALGORITHM_PHASES["ttk"], {"vertex_order", "gradient"})
        self.assertEqual(bench.ALGORITHM_PHASES["reduction_kernel"], {"builder_setup", "gradient"})
        for a in bench.ALGORITHMS:
            components = bench.performance_components(raw, a)
            for i, total in enumerate(raw["algorithms"][a]["performance_seconds"]):
                self.assertAlmostEqual(components["construction_seconds"][i] +
                                       components["algorithm_seconds"][i], total)
        self.assertEqual(raw, original)

    def test_split_rejects_missing_invalid_or_double_counted_performance_phases(self):
        for mutate in (
            lambda r: r.pop("input_loading_seconds"),
            lambda r: r.update(input_loading_seconds=math.nan),
            lambda r: r["algorithms"]["ttk"].pop("performance_phases_seconds"),
            lambda r: r["algorithms"]["f_max"]["performance_phases_seconds"].pop(),
            lambda r: r["algorithms"]["ttk"]["performance_phases_seconds"][0].update(vertex_order=99.0),
            lambda r: r["algorithms"]["ttk"]["performance_phases_seconds"][0].update(gradient=math.inf),
        ):
            raw = self.split_raw()
            mutate(raw)
            with self.assertRaises(ValueError):
                bench.summarize(raw, 3, 1)
        with self.assertRaisesRegex(ValueError, "v2"):
            bench.performance_components(self.raw(), "ttk")

    def test_split_renderer_defaults_to_algorithm_and_separates_loading_construction(self):
        data = dict(schema="resident-gradient-study-v2", completed_utc="test",
                    construction_phases={a: sorted(bench.CONSTRUCTION_PHASES[a]) for a in bench.ALGORITHMS},
                    algorithm_phases={a: sorted(bench.ALGORITHM_PHASES[a]) for a in bench.ALGORITHMS},
                    arguments=dict(repeats=3, diagnostics=1, workers=[1],
                                   terrain_sizes=[3], volume_sizes=[], seeds=[0]),
                    cases=[dict(family="terrain", size=3, seed=0,
                                measurements=[dict(workers=1, raw=self.split_raw(), summary={})])])
        algorithms, phases = render.render_tables(data)
        self.assertIn("10000.000 & 2000.000 & 2000.000 & 1.000", algorithms)
        self.assertIn("Total (boundary clocks)", phases)
        self.assertNotIn("Total (diagnostic)", phases)
        total, _ = render.render_tables(data, "total")
        self.assertIn("101000.000", total)
        construction = render.render_construction(data)
        self.assertIn("Shared loading", construction)
        self.assertIn("100.000", construction)
        data["cases"][0]["measurements"][0]["raw"]["schema"] = "resident-gradient-v1"
        with self.assertRaisesRegex(ValueError, "Mixed"):
            render.render_tables(data)

    def test_v1_evidence_cannot_supply_new_comparison(self):
        data = dict(schema="resident-gradient-study-v1", completed_utc="test",
                    arguments=dict(repeats=3, diagnostics=1, workers=[1],
                                   terrain_sizes=[3], volume_sizes=[], seeds=[0]),
                    cases=[dict(family="terrain", size=3, seed=0,
                                measurements=[dict(workers=1, raw=self.raw(), summary={})])])
        with self.assertRaisesRegex(ValueError, "v2"):
            render.render_tables(data, "algorithm")
        with self.assertRaisesRegex(ValueError, "v2"):
            render.render_construction(data)

    def pls_raw(self):
        raw = self.split_raw()
        raw.update(schema="resident-gradient-v3", euler_characteristic=1,
                   performance_orders=[[bench.PLS_ALGORITHMS[a] for a in bench.PLS_ORDERS[i % 4]]
                                       for i in range(12)])
        raw["algorithms"]["process_lower_stars"] = copy.deepcopy(raw["algorithms"]["f_max"])
        raw["algorithms"]["process_lower_stars"]["diagnostics"][0]["gradient_details_seconds"] = {
            k: 0.5 for k in bench.DETAIL_PHASES["process_lower_stars"]}
        for a in raw["algorithms"].values():
            a["performance_seconds"] *= 4
            a["performance_phases_seconds"] *= 4
        return raw

    def test_pls_partition_ratios_and_count_differences(self):
        raw = self.pls_raw()
        original = copy.deepcopy(raw)
        summary = bench.summarize(raw, 12, 1)
        self.assertEqual(summary["process_lower_stars"]["algorithm_seconds"]["median"], 10)
        self.assertEqual(bench.ALGORITHM_PHASES["process_lower_stars"], {"builder_setup", "gradient"})
        ratio = summary["algorithm_paired_ratios"]["process_lower_stars/ttk"]
        self.assertAlmostEqual(ratio["median"], 100 / 51)
        self.assertIn("paired_repetition_bootstrap_95_interval", ratio)
        self.assertEqual(summary, bench.summarize(raw, 12, 1))
        self.assertEqual(original, raw)
        raw["algorithms"]["process_lower_stars"]["critical_counts"] = [2, 1, 0]
        raw["critical_counts_match"] = False
        bench.summarize(raw, 12, 1)  # Valid difference, same Euler characteristic.
        raw["algorithms"]["process_lower_stars"]["critical_counts"] = [3, 0, 0]
        with self.assertRaisesRegex(ValueError, "Euler"):
            bench.summarize(raw, 12, 1)

    def test_pls_order_balance_and_arguments(self):
        for position in range(4):
            self.assertEqual(sorted(row[position] for row in bench.PLS_ORDERS), list(range(4)))
        adjacent = Counter((a, b) for row in bench.PLS_ORDERS for a, b in zip(row, row[1:]))
        self.assertEqual(adjacent, {(a, b): 1 for a in range(4) for b in range(4) if a != b})
        with tempfile.TemporaryDirectory() as directory:
            args = ["--benchmark", "native", "--input-dir", directory, "--include-pls",
                    "--output", str(Path(directory) / "raw.json")]
            parsed = bench.parse_args(args)
            self.assertEqual((parsed.repeats, parsed.diagnostics), (12, 4))
            for bad in [["--repeats", "6"], ["--diagnostics", "6"], ["--order-offset", "4"]]:
                with redirect_stderr(StringIO()), self.assertRaises(SystemExit):
                    bench.parse_args(args + bad)

    def test_pls_rejects_missing_or_inconsistent_evidence(self):
        for mutate in (
            lambda r: r["algorithms"].pop("process_lower_stars"),
            lambda r: r["performance_orders"].pop(),
            lambda r: r["performance_orders"][0].__setitem__(0, "ttk"),
            lambda r: r.pop("euler_characteristic"),
            lambda r: r["algorithms"]["process_lower_stars"]["diagnostics"][0]["gradient_details_seconds"].update(replay=99.),
            lambda r: r.update(critical_counts_match=False),
        ):
            raw = self.pls_raw()
            mutate(raw)
            with self.assertRaises(ValueError):
                bench.summarize(raw, 12, 1)

    def test_pls_renderer_preserves_all_four_methods(self):
        data = dict(schema="resident-gradient-study-v3", completed_utc="test",
                    construction_phases={a: sorted(bench.CONSTRUCTION_PHASES[a]) for a in bench.PLS_ALGORITHMS},
                    algorithm_phases={a: sorted(bench.ALGORITHM_PHASES[a]) for a in bench.PLS_ALGORITHMS},
                    arguments=dict(repeats=12, diagnostics=1, workers=[1],
                                   terrain_sizes=[3], volume_sizes=[], seeds=[0]),
                    cases=[dict(family="terrain", size=3, seed=0,
                                measurements=[dict(workers=1, raw=self.pls_raw(), summary={})])])
        table, phases = render.render_tables(data)
        self.assertIn("PLS algorithm (ms)", table)
        self.assertIn("PLS/TTK", table)
        self.assertIn("10000.000 & 2000.000 & 2000.000 & 10000.000 & 1.000 & 1.961", table)
        self.assertIn("PLS & builder setup", phases)
        self.assertIn("PLS construction (ms)", render.render_construction(data))


@unittest.skipUnless(os.environ.get("MORSEFRAMES_RESIDENT_BENCHMARK"),
                     "requires the separately built TTK resident benchmark")
class ResidentGradientNativeTest(unittest.TestCase):
    def test_native_pls_reference_phases_and_orders(self):
        executable = Path(os.environ["MORSEFRAMES_RESIDENT_BENCHMARK"])
        examples = (
            "morseframes-ttk-v1 1 3 3\n2 0 1\n0 1\n1 2\n0 2\n",
            "morseframes-ttk-v1 2 3 1\n2 0 1\n0 1 2\n",
            "morseframes-ttk-v1 3 5 2\n3 0 4 2 1\n0 1 2 3\n1 2 3 4\n",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "complex.txt"
            for source in examples:
                path.write_text(source)
                for workers in (1, 4):
                    raw = bench.run_native(executable, path, workers, 4, 4, 0, True, 2)
                    summary = bench.summarize(raw, 4, 4)
                    self.assertEqual(raw["schema"], "resident-gradient-v3")
                    self.assertTrue(raw["exact_reference_checks"])
                    self.assertEqual(raw["performance_orders"][0],
                                     [bench.PLS_ALGORITHMS[a] for a in bench.PLS_ORDERS[2]])
                    self.assertGreater(summary["process_lower_stars"]["gradient_details_seconds"]["lower_star_setup"]["median"], 0)
                    self.assertGreater(summary["process_lower_stars"]["pls_profile_seconds"]["cleanup"]["median"], 0)
                    for algorithm in bench.PLS_ALGORITHMS:
                        c = bench.performance_components(raw, algorithm)
                        for i in range(4):
                            self.assertAlmostEqual(c["construction_seconds"][i] + c["algorithm_seconds"][i],
                                                   raw["algorithms"][algorithm]["performance_seconds"][i])

    def test_native_pls_rejects_ties_without_changing_legacy_mode(self):
        executable = Path(os.environ["MORSEFRAMES_RESIDENT_BENCHMARK"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "complex.txt"
            path.write_text("morseframes-ttk-v1 2 3 1\n0 0 0\n0 1 2\n")
            self.assertEqual(bench.run_native(executable, path, 1, 1, 1, 0)["schema"], "resident-gradient-v2")
            with self.assertRaises(subprocess.CalledProcessError):
                bench.run_native(executable, path, 1, 4, 4, 0, True)

    def test_native_fresh_runs_ties_and_phase_accounting(self):
        executable = Path(os.environ["MORSEFRAMES_RESIDENT_BENCHMARK"])
        examples = (
            "morseframes-ttk-v1 1 3 3\n2 0 1\n0 1\n1 2\n0 2\n",
            "morseframes-ttk-v1 2 3 1\n0 0 0\n0 1 2\n",
            "morseframes-ttk-v1 3 5 2\n3 0 4 2 1\n0 1 2 3\n1 2 3 4\n",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "complex.txt"
            for source in examples:
                path.write_text(source)
                for workers in (1, 4):
                    # Repeated tiny diagnostics exercise emission/callback
                    # boundary accounting as well as exact gradient checks.
                    raw = bench.run_native(executable, path, workers, 1, 6, 0)
                    summary = bench.summarize(raw, 1, 6)
                    self.assertTrue(raw["exact_reference_checks"])
                    self.assertEqual(raw["schema"], "resident-gradient-v2")
                    self.assertGreater(raw["input_loading_seconds"], 0)
                    for algorithm in bench.ALGORITHMS:
                        c = bench.performance_components(raw, algorithm)
                        self.assertAlmostEqual(c["construction_seconds"][0] + c["algorithm_seconds"][0],
                                               raw["algorithms"][algorithm]["performance_seconds"][0])
                    self.assertGreater(summary["ttk"]["phases_seconds"]["vertex_order"]["median"], 0)

    def test_native_rejects_unsupported_or_invalid_inputs(self):
        executable = Path(os.environ["MORSEFRAMES_RESIDENT_BENCHMARK"])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "complex.txt"
            for source in (
                "morseframes-ttk-v1 2 3 1\n0 1 2\n0 0 2\n",  # degenerate
                "morseframes-ttk-v1 2 4 1\n0 1 2 3\n0 1 2\n",  # isolated
                "morseframes-ttk-v1 2 3 1\n0 nan 2\n0 1 2\n",
            ):
                path.write_text(source)
                with self.assertRaises(subprocess.CalledProcessError):
                    bench.run_native(executable, path, 1, 1, 1, 0)


if __name__ == "__main__":
    unittest.main()
