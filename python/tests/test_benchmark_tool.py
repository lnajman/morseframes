import sys
import tempfile
import unittest
from dataclasses import replace
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

import benchmark_persistence as bench  # noqa: E402
import analyze_selector_features as selfeat  # noqa: E402
import calibrate_profile_gate as gatecal  # noqa: E402
import compare_profile_gates as gatecmp  # noqa: E402
import run_fair_profile_validation as fairval  # noqa: E402
import render_synthetic_scale_table as synthscale  # noqa: E402
import summarize_selector_decisions as seld  # noqa: E402


class BenchmarkToolTest(unittest.TestCase):
    def test_run_benchmarks_smoke(self):
        rows = bench.run_benchmarks(
            families=["lower-star", "plateau", "rips"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            plateau_levels=3,
        )

        self.assertEqual(len(rows), 3)
        self.assertEqual({row.family for row in rows}, {"lower-star", "plateau", "rips"})
        for row in rows:
            self.assertEqual(row.size, 5)
            self.assertEqual(row.seed, 0)
            self.assertEqual(row.frame_mode, "fused")
            self.assertIn(row.validation_mode, {"core", "materialized"})
            self.assertEqual(row.barcodes_materialized, row.validation_mode == "materialized")
            self.assertGreater(row.num_simplices, 0)
            if row.family == "plateau":
                self.assertLessEqual(row.num_levels, 3)
            self.assertGreater(row.num_critical_simplices, 0)
            self.assertGreaterEqual(row.reducer_working_set_size, 0)
            self.assertGreaterEqual(row.reducer_initial_nonempty_annotations, 0)
            self.assertGreaterEqual(row.reducer_initial_total_annotation_size, 0)
            self.assertGreaterEqual(row.reducer_initial_max_annotation_size, 0)
            self.assertGreaterEqual(row.reducer_initial_inverse_list_entries, 0)
            self.assertEqual(
                row.reducer_initial_total_annotation_size,
                row.reducer_initial_inverse_list_entries,
            )
            self.assertGreaterEqual(row.reducer_boundary_plan_face_scans, 0)
            self.assertGreaterEqual(row.reducer_boundary_annotation_candidate_criticals, 0)
            self.assertGreaterEqual(row.reducer_boundary_annotation_zero_skipped_criticals, 0)
            self.assertGreaterEqual(row.reducer_boundary_annotation_zero_skipped_faces, 0)
            self.assertGreaterEqual(row.reducer_boundary_annotation_xors, 0)
            self.assertGreaterEqual(row.reducer_boundary_annotation_total_output_size, 0)
            self.assertGreaterEqual(row.reducer_boundary_annotation_max_output_size, 0)
            self.assertGreaterEqual(row.reducer_pivot_eliminations, 0)
            self.assertGreaterEqual(row.reducer_remove_candidate_scans, 0)
            self.assertGreaterEqual(row.reducer_remove_total_annotation_size, 0)
            self.assertGreaterEqual(row.reducer_xor_candidate_scans, 0)
            self.assertGreaterEqual(row.reducer_xor_applied, 0)
            self.assertGreaterEqual(row.reducer_xor_total_input_size, 0)
            self.assertGreaterEqual(row.reducer_xor_total_output_size, 0)
            self.assertGreaterEqual(row.reducer_xor_inserted_labels, 0)
            self.assertGreaterEqual(row.reducer_xor_removed_labels, 0)
            self.assertGreaterEqual(row.sequence_seconds, 0.0)
            self.assertGreaterEqual(row.reference_seconds, 0.0)
            self.assertGreaterEqual(row.morse_reduction_seconds, 0.0)
            self.assertGreaterEqual(row.reducer_setup_seconds, 0.0)
            self.assertGreaterEqual(row.reducer_compute_seconds, 0.0)
            self.assertGreaterEqual(row.morse_seconds, 0.0)
            self.assertGreaterEqual(row.standard_seconds, 0.0)
            self.assertIsNone(row.gudhi_cam_seconds)
            self.assertIsNone(row.gudhi_cam_speedup)
            self.assertIsNone(row.perseus_seconds)
            self.assertIsNone(row.perseus_speedup)
            self.assertIsNone(row.perseus_matches)
            self.assertGreater(row.speedup, 0.0)

    def test_cam_s4_rips_family_smoke(self):
        rows = bench.run_benchmarks(
            families=["cam-s4-rips"],
            sizes=[8],
            seeds=[0],
            repeats=1,
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].family, "cam-s4-rips")
        self.assertEqual(rows[0].size, 8)
        self.assertGreater(rows[0].num_simplices, 0)

    def test_terrain_family_smoke(self):
        rows = bench.run_benchmarks(
            families=["terrain"],
            sizes=[4],
            seeds=[0],
            repeats=1,
            plateau_levels=4,
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].family, "terrain")
        self.assertEqual(rows[0].name, "terrain-g4-L4-seed0")
        self.assertEqual(rows[0].size, 4)
        self.assertEqual(rows[0].num_simplices, 67)
        self.assertLessEqual(rows[0].num_levels, 4)
        self.assertGreater(rows[0].num_critical_simplices, 0)

    def test_image_grid_family_smoke(self):
        rows = bench.run_benchmarks(
            families=["image-grid"],
            sizes=[4],
            seeds=[0],
            repeats=1,
            plateau_levels=4,
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].family, "image-grid")
        self.assertEqual(rows[0].name, "image-grid-g4-L4-seed0")
        self.assertEqual(rows[0].size, 4)
        self.assertEqual(rows[0].num_simplices, 67)
        self.assertLessEqual(rows[0].num_levels, 4)
        self.assertGreater(rows[0].num_critical_simplices, 0)

    def test_run_benchmarks_can_force_separate_frame_mode(self):
        rows = bench.run_benchmarks(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            plateau_levels=3,
            frame_mode="separate",
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].frame_mode, "separate")
        self.assertGreaterEqual(rows[0].reference_seconds, 0.0)

    def test_run_benchmarks_can_expand_sequence_algorithm_portfolio(self):
        rows = bench.run_benchmarks(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            plateau_levels=3,
            sequence_algorithm="portfolio",
        )

        self.assertEqual(len(rows), len(bench.mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO))
        self.assertEqual(
            {row.sequence_algorithm for row in rows},
            set(bench.mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO),
        )

    def test_roadmap_distance_matrix_reader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "distances.txt"
            path.write_text("0 1 2\n1 0 3\n2 3 0\n")

            distances = bench.read_distance_matrix(path)

            self.assertEqual(distances, [[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])

    def test_roadmap_distance_matrix_reader_rejects_non_square_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "distances.txt"
            path.write_text("0 1\n1 0\n2 3\n")

            with self.assertRaises(ValueError):
                bench.read_distance_matrix(path)

    def test_roadmap_family_uses_cached_dataset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache = Path(tmpdir)
            dataset = bench.ROADMAP_DATASETS[50]
            (cache / dataset.filename).write_text(
                "0 1 2 3\n"
                "1 0 4 5\n"
                "2 4 0 6\n"
                "3 5 6 0\n"
            )

            rows = bench.run_benchmarks(
                families=["roadmap-rips"],
                sizes=[50],
                seeds=[0],
                repeats=1,
                roadmap_cache=cache,
            )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].family, "roadmap-rips")
        self.assertEqual(rows[0].name, "roadmap-rips-random50-16d")
        self.assertEqual(rows[0].size, 50)
        self.assertEqual(rows[0].seed, 0)
        self.assertEqual(rows[0].num_simplices, 14)

    def test_roadmap_family_requires_cached_or_downloaded_dataset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(FileNotFoundError):
                bench.run_benchmarks(
                    families=["roadmap-rips"],
                    sizes=[50],
                    seeds=[0],
                    repeats=1,
                    roadmap_cache=Path(tmpdir),
                )

    def test_optional_gudhi_cam_timing(self):
        rows = bench.run_benchmarks(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            plateau_levels=3,
            time_gudhi=True,
        )

        self.assertEqual(len(rows), 1)
        self.assertTrue(rows[0].gudhi_cam_seconds is None or rows[0].gudhi_cam_seconds >= 0.0)
        self.assertTrue(rows[0].gudhi_cam_speedup is None or rows[0].gudhi_cam_speedup >= 0.0)

    def test_perseus_nmfsimtop_writer_uses_positive_integer_births(self):
        complex_ = bench.mp.FilteredComplex.from_simplices([
            ([0], 0.0),
            ([1], 0.0),
            ([0, 1], 2.5),
        ])

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "input_nmfsimtop.txt"
            birth_to_level = bench.write_perseus_nmfsimtop(complex_, path)
            lines = path.read_text().splitlines()

        self.assertEqual(birth_to_level, {1: 0.0, 2: 2.5})
        self.assertEqual(lines[0], "1")
        self.assertIn("0 1 1", lines)
        self.assertIn("0 2 1", lines)
        self.assertIn("1 1 2 2", lines)

    def test_perseus_barcode_reader_maps_integer_births_back_to_levels(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "perseus"
            Path(f"{prefix}_0.txt").write_text("1 2\n1 -1\n")

            finite, essential = bench.read_perseus_barcode(
                prefix,
                birth_to_level={1: 0.0, 2: 3.0},
                max_dimension=0,
            )

        self.assertEqual(finite, ((0, 0.0, 3.0),))
        self.assertEqual(essential, ((0, 0.0),))

    def test_optional_perseus_timing_skips_missing_executable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rows = bench.run_benchmarks(
                families=["plateau"],
                sizes=[5],
                seeds=[0],
                repeats=1,
                plateau_levels=3,
                time_perseus=True,
                perseus_executable=Path(tmpdir) / "missing-perseus",
            )

        self.assertEqual(len(rows), 1)
        self.assertIsNone(rows[0].perseus_seconds)
        self.assertIsNone(rows[0].perseus_speedup)
        self.assertIsNone(rows[0].perseus_matches)

    def test_portfolio_reuses_external_timings_per_complex(self):
        counts = {"gudhi": 0, "perseus": 0}
        original_gudhi = bench.time_gudhi_cam
        original_perseus = bench.time_perseus_persistence

        def fake_gudhi(complex_, *, repeats=3):
            counts["gudhi"] += 1
            return 0.001

        def fake_perseus(complex_, *, repeats=3, executable="perseus"):
            counts["perseus"] += 1
            return bench.PerseusTiming(
                seconds=0.002,
                finite_barcode=(),
                essential_barcode=(),
                matches_standard=True,
            )

        try:
            bench.time_gudhi_cam = fake_gudhi
            bench.time_perseus_persistence = fake_perseus
            rows = bench.run_benchmarks(
                families=["plateau"],
                sizes=[5],
                seeds=[0],
                repeats=1,
                plateau_levels=3,
                sequence_algorithm="portfolio",
                time_gudhi=True,
                time_perseus=True,
            )
        finally:
            bench.time_gudhi_cam = original_gudhi
            bench.time_perseus_persistence = original_perseus

        self.assertEqual(len(rows), len(bench.mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO))
        self.assertEqual(counts, {"gudhi": 1, "perseus": 1})
        for row in rows:
            self.assertEqual(row.gudhi_cam_seconds, 0.001)
            self.assertEqual(row.perseus_seconds, 0.002)
            self.assertTrue(row.perseus_matches)

    def test_profile_vs_measured_smoke(self):
        rows = bench.run_profile_vs_measured(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            profile_repeats=1,
            plateau_levels=3,
            sequence_algorithm=["saturated", "same-level-reduction"],
            profile_selection_metric="profile_total_work",
        )

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.family, "plateau")
        self.assertEqual(row.candidate_count, 2)
        self.assertEqual(row.profile_candidate_count, 2)
        self.assertEqual(row.profile_candidate_gate, "all")
        self.assertIn(
            row.profile_selected_algorithm,
            {"saturated", "same-level-reduction"},
        )
        self.assertIn(
            row.measured_best_algorithm,
            {"saturated", "same-level-reduction"},
        )
        self.assertGreaterEqual(row.profile_selected_measured_rank, 1)
        self.assertGreaterEqual(row.measured_best_profile_rank, 1)
        self.assertGreaterEqual(row.profile_penalty_ratio, 1.0)
        self.assertGreaterEqual(row.profile_penalty_percent, 0.0)
        self.assertGreaterEqual(row.total_profile_seconds, 0.0)
        self.assertGreaterEqual(row.total_measured_morse_seconds, 0.0)
        self.assertGreaterEqual(row.measured_best_sequence_seconds, 0.0)
        self.assertGreaterEqual(row.profile_selected_sequence_seconds, 0.0)
        self.assertGreaterEqual(row.measured_best_reducer_work, 0.0)
        self.assertGreaterEqual(row.profile_selected_reducer_work, 0.0)

    def test_adaptive_profile_selection_metric_resolves_by_family(self):
        self.assertEqual(
            bench.resolve_profile_selection_metric("plateau", "adaptive"),
            "critical_ratio",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("terrain", "adaptive"),
            "critical_ratio",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("image-grid", "adaptive"),
            "critical_ratio",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("cam-s4-rips", "adaptive"),
            "critical_ratio",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("roadmap-rips", "adaptive"),
            "profile_total_work",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("plateau", "adaptive_structured"),
            "critical_ratio",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("terrain", "adaptive_structured"),
            "critical_ratio",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("image-grid", "adaptive_structured"),
            "profile_total_work",
        )
        self.assertEqual(
            bench.resolve_profile_selection_metric("roadmap-rips", "adaptive_structured"),
            "profile_total_work",
        )

    def test_profile_vs_measured_adaptive_metric_smoke(self):
        rows = bench.run_profile_vs_measured(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            profile_repeats=1,
            plateau_levels=3,
            sequence_algorithm=["f-max", "f-min"],
            profile_selection_metric="adaptive",
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].profile_selection_metric, "adaptive")
        self.assertEqual(rows[0].effective_profile_selection_metric, "critical_ratio")
        self.assertIn(rows[0].profile_selected_algorithm, {"f-max", "f-min"})

    def test_profile_vs_measured_can_evaluate_multiple_metrics_once(self):
        rows = bench.run_profile_vs_measured(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            profile_repeats=1,
            plateau_levels=3,
            sequence_algorithm=["f-max", "f-min"],
            profile_selection_metric=["profile_total_work", "adaptive"],
        )

        self.assertEqual(len(rows), 2)
        self.assertEqual(
            {row.profile_selection_metric for row in rows},
            {"profile_total_work", "adaptive"},
        )
        self.assertEqual(
            {
                row.profile_selection_metric: row.effective_profile_selection_metric
                for row in rows
            },
            {
                "profile_total_work": "profile_total_work",
                "adaptive": "critical_ratio",
            },
        )
        self.assertEqual({row.measured_best_morse_seconds for row in rows}, {rows[0].measured_best_morse_seconds})
        self.assertEqual({row.total_measured_morse_seconds for row in rows}, {rows[0].total_measured_morse_seconds})

    def test_family_aware_profile_candidate_gate(self):
        algorithms = bench.mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO

        self.assertEqual(
            bench.profile_candidate_algorithm_options(
                family="plateau",
                size=16,
                algorithms=algorithms,
                gate="family-aware",
            ),
            (
                bench.mp.F_MAX_SEQUENCE,
                bench.mp.F_MIN_SEQUENCE,
                bench.mp.SATURATED_SEQUENCE,
            ),
        )
        self.assertEqual(
            bench.profile_candidate_algorithm_options(
                family="terrain",
                size=8,
                algorithms=algorithms,
                gate="family-aware",
            ),
            (
                bench.mp.F_MAX_SEQUENCE,
                bench.mp.F_MIN_SEQUENCE,
                bench.mp.SATURATED_SEQUENCE,
            ),
        )
        self.assertEqual(
            bench.profile_candidate_algorithm_options(
                family="image-grid",
                size=8,
                algorithms=algorithms,
                gate="family-aware",
            ),
            (
                bench.mp.F_MAX_SEQUENCE,
                bench.mp.F_MIN_SEQUENCE,
                bench.mp.SATURATED_SEQUENCE,
            ),
        )
        self.assertEqual(
            bench.profile_candidate_algorithm_options(
                family="roadmap-rips",
                size=100,
                algorithms=algorithms,
                gate="family-aware",
            ),
            (bench.mp.F_MAX_SEQUENCE, bench.mp.COREDUCTION_SEQUENCE),
        )

    def test_complex_shape_profile_and_shape_aware_gate(self):
        complex_ = bench.mp.FilteredComplex.from_simplices(
            [
                ([0], 0.0),
                ([1], 0.0),
                ([2], 0.0),
                ([0, 1], 1.0),
                ([0, 2], 1.0),
                ([1, 2], 1.0),
                ([0, 1, 2], 1.0),
            ]
        )
        shape = bench.profile_complex_shape(complex_)

        self.assertEqual(shape.max_dimension, 2)
        self.assertEqual(shape.num_vertices, 3)
        self.assertEqual(shape.num_edges, 3)
        self.assertEqual(shape.num_triangles, 1)
        self.assertAlmostEqual(shape.simplex_vertex_ratio, 7.0 / 3.0)
        self.assertAlmostEqual(shape.vertex_sqrt, 3.0 ** 0.5)
        self.assertAlmostEqual(shape.edge_density, 1.0)
        self.assertAlmostEqual(shape.triangle_density, 1.0)
        self.assertEqual(shape.max_coboundary_size, 2)
        self.assertEqual(shape.largest_level_size, 4)
        self.assertAlmostEqual(shape.avg_level_size, 3.5)
        self.assertAlmostEqual(shape.largest_level_ratio, 4.0 / 7.0)
        self.assertAlmostEqual(shape.singleton_level_ratio, 0.0)
        self.assertAlmostEqual(shape.level_count_ratio, 2.0 / 7.0)
        self.assertAlmostEqual(shape.level_concentration, 25.0 / 49.0)
        self.assertGreater(shape.level_entropy_ratio, 0.98)
        self.assertLessEqual(shape.level_entropy_ratio, 1.0)
        self.assertEqual(
            bench.profile_candidate_algorithm_options(
                family="lower-star",
                size=3,
                algorithms=bench.mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
                gate="shape-aware",
                shape=shape,
            ),
            (bench.mp.F_MAX_SEQUENCE, bench.mp.COREDUCTION_SEQUENCE),
        )

    def test_tabulated_profile_candidate_gate_uses_table_keys(self):
        complex_ = bench.mp.FilteredComplex.from_simplices(
            [
                ([0], 0.0),
                ([1], 0.0),
                ([2], 0.0),
                ([0, 1], 1.0),
                ([0, 2], 1.0),
                ([1, 2], 1.0),
                ([0, 1, 2], 1.0),
            ]
        )
        shape = bench.profile_complex_shape(complex_)
        candidate_table = {
            "lower-star:complete-2-skeleton": (bench.mp.F_MIN_SEQUENCE,),
            "lower-star": (bench.mp.F_MAX_SEQUENCE,),
        }

        self.assertEqual(
            bench.profile_candidate_algorithm_options(
                family="lower-star",
                size=3,
                algorithms=bench.mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
                gate="tabulated",
                shape=shape,
                candidate_table=candidate_table,
            ),
            (bench.mp.F_MIN_SEQUENCE,),
        )

    def test_profile_candidate_table_reader(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "gate.csv"
            path.write_text(
                "group,proposed_candidates\n"
                f"rips,{bench.mp.F_MAX_SEQUENCE};{bench.mp.F_MIN_SEQUENCE}\n"
            )

            table = bench.read_profile_candidate_table(path)

        self.assertEqual(
            table,
            {"rips": (bench.mp.F_MAX_SEQUENCE, bench.mp.F_MIN_SEQUENCE)},
        )

    def test_profile_vs_measured_family_aware_gate(self):
        rows = bench.run_profile_vs_measured(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            profile_repeats=1,
            plateau_levels=3,
            sequence_algorithm="portfolio",
            profile_selection_metric="profile_total_work",
            profile_candidate_gate="family-aware",
        )

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.candidate_count, len(bench.mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO))
        self.assertEqual(row.profile_candidate_count, 3)
        self.assertEqual(row.profile_candidate_gate, "family-aware")
        self.assertIn(
            row.profile_selected_algorithm,
            {
                bench.mp.F_MAX_SEQUENCE,
                bench.mp.F_MIN_SEQUENCE,
                bench.mp.SATURATED_SEQUENCE,
            },
        )
        self.assertGreaterEqual(row.measured_best_profile_rank, 0)

    def test_profile_vs_measured_tabulated_gate(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "gate.csv"
            path.write_text(
                "group,proposed_candidates\n"
                f"plateau,{bench.mp.F_MAX_SEQUENCE};{bench.mp.F_MIN_SEQUENCE}\n"
            )
            rows = bench.run_profile_vs_measured(
                families=["plateau"],
                sizes=[5],
                seeds=[0],
                repeats=1,
                profile_repeats=1,
                plateau_levels=3,
                sequence_algorithm="portfolio",
                profile_selection_metric="profile_total_work",
                profile_candidate_gate="tabulated",
                profile_candidate_table=path,
            )

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.profile_candidate_gate, "tabulated")
        self.assertEqual(row.profile_candidate_count, 2)
        self.assertIn(
            row.profile_selected_algorithm,
            {bench.mp.F_MAX_SEQUENCE, bench.mp.F_MIN_SEQUENCE},
        )

    def test_profile_vs_measured_shape_aware_gate(self):
        rows = bench.run_profile_vs_measured(
            families=["rips"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            profile_repeats=1,
            plateau_levels=3,
            sequence_algorithm="portfolio",
            profile_selection_metric="profile_total_work",
            profile_candidate_gate="shape-aware",
        )

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row.profile_candidate_gate, "shape-aware")
        self.assertGreaterEqual(row.shape_max_dimension, 1)
        self.assertGreater(row.shape_num_vertices, 0)
        self.assertGreater(row.shape_simplex_vertex_ratio, 0.0)
        self.assertGreater(row.shape_vertex_sqrt, 0.0)
        self.assertGreaterEqual(row.shape_edge_density, 0.0)
        self.assertGreater(row.shape_largest_level_size, 0)
        self.assertGreater(row.shape_avg_level_size, 0.0)
        self.assertGreaterEqual(row.shape_singleton_level_ratio, 0.0)
        self.assertGreaterEqual(row.shape_level_concentration, 0.0)
        self.assertLessEqual(row.shape_level_entropy_ratio, 1.0)
        self.assertLessEqual(row.profile_candidate_count, row.candidate_count)

    def test_gate_calibration_shape_bins_sparse_graphs(self):
        row = {
            "shape_max_dimension": "1",
            "shape_edge_density": "0.03",
            "shape_triangle_density": "0.0",
            "shape_largest_level_ratio": "0.0",
            "shape_max_coboundary_size": "3",
        }

        self.assertEqual(gatecal.shape_bin(row), "graph-sparse")

    def test_gate_calibration_proposes_winners_by_family_shape(self):
        rows = [
            {
                "family": "roadmap-rips",
                "size": "50",
                "measured_best_algorithm": bench.mp.F_MAX_SEQUENCE,
                "measured_best_morse_seconds": "0.000010",
                "shape_max_dimension": "2",
                "shape_edge_density": "1.0",
                "shape_triangle_density": "1.0",
                "shape_largest_level_ratio": "0.0",
                "shape_max_coboundary_size": "12",
            },
            {
                "family": "roadmap-rips",
                "size": "50",
                "measured_best_algorithm": bench.mp.COREDUCTION_SEQUENCE,
                "measured_best_morse_seconds": "0.000020",
                "shape_max_dimension": "2",
                "shape_edge_density": "1.0",
                "shape_triangle_density": "1.0",
                "shape_largest_level_ratio": "0.0",
                "shape_max_coboundary_size": "12",
            },
        ]

        proposals = gatecal.propose_candidate_gates(rows, group_by="family-shape")

        self.assertEqual(len(proposals), 1)
        proposal = proposals[0]
        self.assertEqual(proposal.group, "roadmap-rips:complete-2-skeleton")
        self.assertEqual(proposal.cases, 2)
        self.assertEqual(proposal.candidate_count, 2)
        self.assertEqual(proposal.coverage_percent, 100.0)
        self.assertEqual(
            set(proposal.proposed_candidates.split(";")),
            {bench.mp.F_MAX_SEQUENCE, bench.mp.COREDUCTION_SEQUENCE},
        )

    def test_gate_calibration_candidate_limit_trades_coverage(self):
        rows = [
            {"family": "rips", "size": "16", "measured_best_algorithm": bench.mp.F_MAX_SEQUENCE},
            {"family": "rips", "size": "16", "measured_best_algorithm": bench.mp.F_MAX_SEQUENCE},
            {
                "family": "rips",
                "size": "16",
                "measured_best_algorithm": bench.mp.COREDUCTION_SEQUENCE,
            },
        ]

        proposals = gatecal.propose_candidate_gates(rows, group_by="family", max_candidates=1)

        self.assertEqual(len(proposals), 1)
        self.assertEqual(proposals[0].proposed_candidates, bench.mp.F_MAX_SEQUENCE)
        self.assertAlmostEqual(proposals[0].coverage_percent, 100.0 * 2.0 / 3.0)

    def test_gate_calibration_guard_candidates_are_appended(self):
        rows = [
            {
                "family": "lower-star",
                "size": "16",
                "measured_best_algorithm": bench.mp.F_MAX_SEQUENCE,
            }
        ]

        proposals = gatecal.propose_candidate_gates(
            rows,
            group_by="family",
            include_candidates={"lower-star": (bench.mp.F_MIN_SEQUENCE,)},
        )

        self.assertEqual(
            proposals[0].proposed_candidates,
            f"{bench.mp.F_MAX_SEQUENCE};{bench.mp.F_MIN_SEQUENCE}",
        )
        self.assertEqual(proposals[0].coverage_percent, 100.0)

    def test_profile_gate_comparison_summarizes_unavailable_winners(self):
        rows = [
            {
                "family": "plateau",
                "profile_candidate_count": "2",
                "profile_matches_measured": "True",
                "measured_best_profile_rank": "1",
                "profile_penalty_percent": "0.0",
                "profile_selected_measured_rank": "1",
                "total_profile_seconds": "0.1",
                "total_measured_morse_seconds": "1.0",
            },
            {
                "family": "plateau",
                "profile_candidate_count": "2",
                "profile_matches_measured": "False",
                "measured_best_profile_rank": "0",
                "profile_penalty_percent": "25.0",
                "profile_selected_measured_rank": "3",
                "total_profile_seconds": "0.2",
                "total_measured_morse_seconds": "1.0",
            },
        ]

        summaries = gatecmp.summarize_gate_rows("guarded", rows)

        plateau = next(row for row in summaries if row.family == "plateau")
        all_row = next(row for row in summaries if row.family == "all")
        self.assertEqual(plateau.label, "guarded")
        self.assertEqual(plateau.cases, 2)
        self.assertEqual(plateau.avg_profile_candidates, 2.0)
        self.assertEqual(plateau.match_percent, 50.0)
        self.assertEqual(plateau.unavailable_count, 1)
        self.assertEqual(plateau.unavailable_percent, 50.0)
        self.assertEqual(plateau.avg_penalty_percent, 12.5)
        self.assertEqual(plateau.max_penalty_percent, 25.0)
        self.assertEqual(plateau.avg_selected_rank, 2.0)
        self.assertAlmostEqual(plateau.profile_overhead_percent, 15.0)
        self.assertEqual(all_row.cases, 2)

    def test_profile_gate_comparison_parses_labeled_input_specs(self):
        label, path = gatecmp.parse_input_spec("expanded=work/table.csv")

        self.assertEqual(label, "expanded")
        self.assertEqual(path, Path("work/table.csv"))

    def test_profile_gate_comparison_can_label_by_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "metrics.csv"
            path.write_text(
                "family,profile_selection_metric,profile_candidate_count,"
                "profile_matches_measured,measured_best_profile_rank,"
                "profile_penalty_percent,profile_selected_measured_rank,"
                "total_profile_seconds,total_measured_morse_seconds\n"
                "plateau,profile_total_work,2,True,1,0.0,1,0.1,1.0\n"
                "plateau,adaptive,2,False,2,10.0,2,0.1,1.0\n"
            )

            summaries = gatecmp.compare_gate_outputs(
                [("metrics", path)],
                include_all=False,
                label_column="profile_selection_metric",
            )

        self.assertEqual(
            {(row.label, row.family, row.avg_penalty_percent) for row in summaries},
            {("adaptive", "plateau", 10.0), ("profile_total_work", "plateau", 0.0)},
        )

    def test_profile_gate_comparison_can_group_by_family_size(self):
        rows = [
            {
                "family": "plateau",
                "size": "16",
                "profile_candidate_count": "3",
                "profile_matches_measured": "True",
                "measured_best_profile_rank": "1",
                "profile_penalty_percent": "0.0",
                "profile_selected_measured_rank": "1",
                "total_profile_seconds": "0.1",
                "total_measured_morse_seconds": "1.0",
            },
            {
                "family": "plateau",
                "size": "24",
                "profile_candidate_count": "4",
                "profile_matches_measured": "False",
                "measured_best_profile_rank": "2",
                "profile_penalty_percent": "10.0",
                "profile_selected_measured_rank": "2",
                "total_profile_seconds": "0.2",
                "total_measured_morse_seconds": "1.0",
            },
        ]

        summaries = gatecmp.summarize_gate_rows(
            "plateau-gate",
            rows,
            include_all=False,
            group_by="family-size",
        )

        self.assertEqual([row.family for row in summaries], ["plateau:n=16", "plateau:n=24"])
        self.assertEqual([row.avg_profile_candidates for row in summaries], [3.0, 4.0])

    def test_profile_vs_measured_summary_writer(self):
        rows = bench.run_profile_vs_measured(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            profile_repeats=1,
            plateau_levels=3,
            sequence_algorithm=["saturated", "same-level-reduction"],
        )
        output = StringIO()

        bench.write_profile_vs_measured_summary(rows, output)

        text = output.getvalue()
        self.assertIn("match%", text)
        self.assertIn("avg_pen%", text)
        self.assertIn("eff_metric", text)
        self.assertIn("prof_cands", text)
        self.assertIn("prof_over%", text)
        self.assertIn("plateau", text)

    def test_fair_profile_validation_builds_repeatable_recipe(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--roadmap-cache",
                str(Path(tmpdir) / "roadmap-data"),
            ])
            commands = fairval.build_all_commands(args)

        self.assertEqual(len(commands), 7)
        benchmark_commands = [
            command for command in commands
            if "benchmark_persistence.py" in command[1]
        ]
        summary_commands = [
            command for command in commands
            if "compare_profile_gates.py" in command[1]
        ]
        self.assertEqual(len(benchmark_commands), 3)
        self.assertEqual(len(summary_commands), 4)
        self.assertIn("--profile-selection-metrics", benchmark_commands[0])
        self.assertIn("profile_total_work", benchmark_commands[0])
        self.assertIn("critical_ratio", benchmark_commands[0])
        self.assertIn("adaptive", benchmark_commands[0])
        self.assertIn("--time-gudhi-cam", benchmark_commands[0])
        self.assertTrue(
            any(
                "effective_profile_selection_metric" in command
                for command in summary_commands
            )
        )

    def test_fair_profile_validation_can_regenerate_summaries_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--summaries-only",
            ])
            commands = fairval.build_all_commands(args)

        self.assertEqual(len(commands), 4)
        self.assertTrue(
            all("compare_profile_gates.py" in command[1] for command in commands)
        )

    def test_fair_profile_validation_extended_holdout_uses_separate_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--validation-preset",
                "extended-holdout",
            ])
            commands = fairval.build_all_commands(args)
            splits = fairval.validation_splits("extended-holdout")

        self.assertEqual(splits[0].label, "base-holdout")
        self.assertEqual(splits[0].sizes, (20, 28, 36, 44))
        self.assertEqual(splits[0].seeds, (20, 21, 22, 23, 24))
        self.assertIn("extended_holdout", splits[0].output_name)
        self.assertIn("extended_holdout", str(args.table_output))
        self.assertIn("extended_holdout", str(args.prose_output))
        self.assertIn("extended_holdout", str(args.manifest_output))
        self.assertTrue(
            any(
                "profile_metric_extended_holdout_comparison.csv"
                in " ".join(command)
                for command in commands
            )
        )

    def test_fair_profile_validation_plateau_holdout_uses_broader_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--validation-preset",
                "plateau-holdout",
            ])
            commands = fairval.build_all_commands(args)
            splits = fairval.validation_splits("plateau-holdout")

        self.assertEqual(len(splits), 1)
        self.assertEqual(splits[0].label, "plateau-holdout")
        self.assertEqual(splits[0].families, ("plateau",))
        self.assertEqual(splits[0].sizes, (12, 16, 20, 24, 28, 32, 36, 40, 44, 48))
        self.assertIn("boundary_annotation_work", args.profile_selection_metrics)
        self.assertIn("initial_annotation_size", args.profile_selection_metrics)
        self.assertIn("adaptive_structured", args.profile_selection_metrics)
        self.assertIn("plateau_holdout", str(args.table_output))
        self.assertTrue(
            any(
                "profile_vs_measured_plateau_holdout_multi_metric.csv"
                in " ".join(command)
                for command in commands
            )
        )

    def test_fair_profile_validation_terrain_holdout_uses_separate_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--validation-preset",
                "terrain-holdout",
            ])
            commands = fairval.build_all_commands(args)
            splits = fairval.validation_splits("terrain-holdout")

        self.assertEqual(len(splits), 1)
        self.assertEqual(splits[0].label, "terrain-holdout")
        self.assertEqual(splits[0].families, ("terrain",))
        self.assertEqual(splits[0].sizes, (6, 8, 10, 12))
        self.assertEqual(splits[0].seeds, (40, 41, 42, 43, 44))
        self.assertIn("boundary_annotation_work", args.profile_selection_metrics)
        self.assertIn("adaptive_structured", args.profile_selection_metrics)
        self.assertIn("terrain_holdout", str(args.table_output))
        self.assertTrue(
            any(
                "profile_vs_measured_terrain_holdout_multi_metric.csv"
                in " ".join(command)
                for command in commands
            )
        )

    def test_fair_profile_validation_image_holdout_uses_separate_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--validation-preset",
                "image-holdout",
            ])
            commands = fairval.build_all_commands(args)
            splits = fairval.validation_splits("image-holdout")

        self.assertEqual(len(splits), 1)
        self.assertEqual(splits[0].label, "image-holdout")
        self.assertEqual(splits[0].families, ("image-grid",))
        self.assertEqual(splits[0].sizes, (6, 8, 10, 12))
        self.assertEqual(splits[0].seeds, (50, 51, 52, 53, 54))
        self.assertIn("boundary_annotation_work", args.profile_selection_metrics)
        self.assertIn("adaptive_structured", args.profile_selection_metrics)
        self.assertIn("image_holdout", str(args.table_output))
        self.assertTrue(
            any(
                "profile_vs_measured_image_holdout_multi_metric.csv"
                in " ".join(command)
                for command in commands
            )
        )

    def test_fair_profile_validation_image_stress_uses_larger_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--validation-preset",
                "image-stress",
            ])
            commands = fairval.build_all_commands(args)
            splits = fairval.validation_splits("image-stress")

        self.assertEqual(len(splits), 1)
        self.assertEqual(splits[0].label, "image-stress")
        self.assertEqual(splits[0].families, ("image-grid",))
        self.assertEqual(splits[0].sizes, (14, 18, 22, 26))
        self.assertEqual(splits[0].seeds, (60, 61, 62, 63, 64))
        self.assertIn("adaptive_structured", args.profile_selection_metrics)
        self.assertIn("image_stress", str(args.table_output))
        self.assertTrue(
            any(
                "profile_vs_measured_image_stress_multi_metric.csv"
                in " ".join(command)
                for command in commands
            )
        )

    def test_fair_profile_validation_renders_latex_table_from_raw_csvs(self):
        header = (
            "profile_selection_metric,effective_profile_selection_metric,"
            "profile_matches_measured,measured_best_profile_rank,"
            "profile_penalty_percent,profile_selected_measured_rank\n"
        )
        rows = (
            "profile_total_work,profile_total_work,True,1,0.0,1\n"
            "profile_total_work,profile_total_work,False,1,10.0,2\n"
            "critical_ratio,critical_ratio,True,1,0.0,1\n"
            "adaptive,critical_ratio,True,1,0.0,1\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)
            for split in fairval.default_validation_splits():
                (output_dir / split.output_name).write_text(header + rows)

            table_rows = fairval.summarize_fair_metric_rows(
                splits=fairval.default_validation_splits(),
                output_dir=output_dir,
                metric_order=fairval.DEFAULT_PROFILE_METRICS,
            )
            latex = fairval.render_fair_metric_table(table_rows)
            prose = fairval.render_fair_metric_prose(table_rows)
            holdout_prose = fairval.render_fair_metric_prose(
                replace(row, split=f"{row.split}-holdout")
                for row in table_rows
            )

        self.assertIn("Base & PTW (PTW) & 2 & 50.0", latex)
        self.assertIn("Base & Adapt (CR) & 1 & 100.0", latex)
        self.assertIn(r"\label{tab:fair-profile-metrics}", latex)
        self.assertIn(r"from $5.00\%$ for \texttt{profile-total-work}", prose)
        self.assertIn(r"to $0.00\%$ for \texttt{critical-ratio}/\texttt{adaptive}", prose)
        self.assertIn(r"from $5.00\%$ for \texttt{profile-total-work}", holdout_prose)

    def test_fair_profile_validation_renders_plateau_prose(self):
        rows = [
            fairval.FairMetricTableRow(
                split="plateau-holdout",
                requested_metric="profile_total_work",
                effective_metric="PTW",
                cases=4,
                match_percent=50.0,
                unavailable_count=0,
                avg_penalty_percent=2.0,
                max_penalty_percent=4.0,
                avg_selected_rank=1.5,
            ),
            fairval.FairMetricTableRow(
                split="plateau-holdout",
                requested_metric="critical_ratio",
                effective_metric="CR",
                cases=4,
                match_percent=25.0,
                unavailable_count=0,
                avg_penalty_percent=5.0,
                max_penalty_percent=10.0,
                avg_selected_rank=2.0,
            ),
        ]

        prose = fairval.render_fair_metric_prose(rows)

        self.assertIn("plateau-focused hold-out", prose)
        self.assertIn(r"\texttt{profile-total-work}", prose)
        self.assertIn("different effective metric", prose)

    def test_fair_profile_validation_renders_image_stress_label(self):
        rows = [
            fairval.FairMetricTableRow(
                split="image-stress",
                requested_metric="profile_total_work",
                effective_metric="PTW",
                cases=4,
                match_percent=50.0,
                unavailable_count=0,
                avg_penalty_percent=1.0,
                max_penalty_percent=4.0,
                avg_selected_rank=1.5,
            ),
            fairval.FairMetricTableRow(
                split="image-stress",
                requested_metric="critical_ratio",
                effective_metric="CR",
                cases=4,
                match_percent=25.0,
                unavailable_count=0,
                avg_penalty_percent=2.0,
                max_penalty_percent=10.0,
                avg_selected_rank=2.0,
            ),
        ]

        prose = fairval.render_fair_metric_prose(rows)

        self.assertIn("image stress split", prose)
        self.assertNotIn("image-grid hold-out", prose)

    def test_fair_profile_validation_renders_experiment_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--roadmap-cache",
                str(Path(tmpdir) / "roadmap-data"),
                "--summaries-only",
            ])
            commands = fairval.build_all_commands(args)
            full_args = fairval.parse_args([
                "--python-executable",
                "python-test",
                "--output-dir",
                tmpdir,
                "--roadmap-cache",
                str(Path(tmpdir) / "roadmap-data"),
            ])
            manifest = fairval.render_experiment_manifest(
                args=args,
                splits=fairval.default_validation_splits(),
                commands=commands,
                full_validation_commands=fairval.build_all_commands(full_args),
            )

        self.assertIn("# Fair Profile Validation Manifest", manifest)
        self.assertIn("profile_vs_measured_base_validation_multi_metric.csv", manifest)
        self.assertIn("profile_metric_fair_validation_comparison.csv", manifest)
        self.assertIn("profile_metric_fair_validation_table.tex", manifest)
        self.assertIn("Run from the repository workspace root", manifest)
        self.assertIn("benchmark_persistence.py", manifest)
        self.assertIn("compare_profile_gates.py", manifest)
        self.assertIn("Commands Executed by This Invocation", manifest)

    def test_selector_decision_summary_keeps_current_tied_best(self):
        csv_text = (
            "label,family,cases,match_percent,unavailable_count,"
            "avg_penalty_percent,max_penalty_percent,avg_selected_rank\n"
            "base:adaptive,all,10,70,0,1.0,5.0,1.2\n"
            "base:critical_ratio,all,10,70,0,1.0,5.0,1.2\n"
            "base:profile_total_work,all,10,40,0,4.0,12.0,2.1\n"
            "roadmap:adaptive,all,3,100,0,0.0,0.0,1.0\n"
            "roadmap:critical_ratio,all,3,0,0,7.0,11.0,2.0\n"
            "roadmap:profile_total_work,all,3,100,0,0.0,0.0,1.0\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "comparison.csv"
            path.write_text(csv_text)

            rows = seld.summarize_decisions([("report", path)], current_metric="adaptive")
            latex = seld.render_latex_table(rows)
            prose = seld.render_prose(rows)

        self.assertEqual(len(rows), 2)
        self.assertEqual({row.decision for row in rows}, {"keep-current"})
        self.assertEqual(rows[0].source, "report")
        self.assertEqual(rows[0].split, "base")
        self.assertIn("adaptive", rows[0].best_metrics)
        self.assertIn(r"\label{tab:selector-decision-summary}", latex)
        self.assertIn("keeps the current", prose)

    def test_selector_decision_summary_marks_switch_candidate(self):
        csv_text = (
            "label,family,cases,match_percent,unavailable_count,"
            "avg_penalty_percent,max_penalty_percent,avg_selected_rank\n"
            "adaptive,all,8,50,0,2.0,9.0,1.8\n"
            "profile_total_work,all,8,70,0,0.8,3.0,1.2\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "plateau.csv"
            path.write_text(csv_text)

            rows = seld.summarize_decisions(
                [("plateau-holdout", path)],
                current_metric="adaptive",
                promotion_margin=0.5,
            )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].split, "plateau-holdout")
        self.assertEqual(rows[0].best_metrics, "profile_total_work")
        self.assertEqual(rows[0].decision, "switch-candidate")
        self.assertAlmostEqual(rows[0].gap_percent_points, 1.2)

    def test_selector_feature_diagnostic_compares_critical_with_time(self):
        csv_text = (
            "family,name,seed,size,num_simplices,profile_selection_metric,"
            "profile_selected_algorithm,measured_best_algorithm,profile_penalty_percent,"
            "shape_simplex_vertex_ratio,shape_largest_level_size,"
            "shape_largest_level_ratio,shape_level_concentration,"
            "shape_level_entropy_ratio\n"
            "image-grid,a,0,6,100,critical_ratio,f-min,f-max,4.0,"
            "4.0,70,0.70,0.60,0.80\n"
            "image-grid,a,0,6,100,profile_total_work,f-max,f-max,1.0,"
            "4.0,70,0.70,0.60,0.80\n"
            "image-grid,a,0,6,100,profile_seconds,f-min,f-max,2.0,"
            "4.0,70,0.70,0.60,0.80\n"
            "image-grid,b,1,8,160,critical_ratio,f-max,f-max,0.5,"
            "5.0,90,0.56,0.45,0.90\n"
            "image-grid,b,1,8,160,profile_total_work,f-min,f-max,2.0,"
            "5.0,90,0.56,0.45,0.90\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "profiles.csv"
            path.write_text(csv_text)

            cases = selfeat.build_selector_cases([("toy", path)])
            summaries = selfeat.summarize_selector_cases(cases)
            latex = selfeat.render_latex_table(summaries)
            prose = selfeat.render_prose(summaries)

        self.assertEqual(len(cases), 2)
        self.assertEqual([case.winner for case in cases], ["time", "critical"])
        self.assertEqual(cases[0].time_metric, "profile_total_work")
        toy_summary = next(
            row for row in summaries
            if row.split == "toy" and row.family == "image-grid"
        )
        self.assertEqual(toy_summary.time_wins, 1)
        self.assertEqual(toy_summary.critical_wins, 1)
        self.assertAlmostEqual(toy_summary.avg_margin_percent_points, 0.75)
        self.assertAlmostEqual(toy_summary.descriptor_coverage_percent, 100.0)
        self.assertAlmostEqual(toy_summary.avg_simplex_vertex_ratio, 4.5)
        self.assertIn(r"\label{tab:selector-feature-diagnostic}", latex)
        self.assertIn("selector-feature diagnostic", prose)

    def test_synthetic_scale_renderer_includes_core_strategies(self):
        rows = [
            {
                "family": "lower-star",
                "sequence_algorithm": "f-max",
                "num_simplices": "100",
                "num_critical_simplices": "20",
                "morse_seconds": "0.000010",
                "speedup": "2.0",
            },
            {
                "family": "lower-star",
                "sequence_algorithm": "f-min",
                "num_simplices": "100",
                "num_critical_simplices": "30",
                "morse_seconds": "0.000020",
                "speedup": "1.0",
            },
            {
                "family": "lower-star",
                "sequence_algorithm": "saturated",
                "num_simplices": "100",
                "num_critical_simplices": "40",
                "morse_seconds": "0.000030",
                "speedup": "0.8",
            },
            {
                "family": "lower-star",
                "sequence_algorithm": "same-level-reduction",
                "num_simplices": "100",
                "num_critical_simplices": "50",
                "morse_seconds": "0.000040",
                "speedup": "0.7",
            },
        ]

        summaries = synthscale.summarize_rows(
            rows,
            families=["lower-star"],
            strategies=["f-max", "f-min", "saturated", "same-level-reduction"],
        )
        latex = synthscale.render_latex_table(summaries)
        prose = synthscale.render_prose(summaries)

        self.assertEqual([row.strategy for row in summaries], [
            "f-max",
            "f-min",
            "saturated",
            "same-level-reduction",
        ])
        self.assertAlmostEqual(summaries[0].critical_percent, 20.0)
        self.assertAlmostEqual(summaries[0].morse_us, 10.0)
        self.assertIn(r"\texttt{f-max}", latex)
        self.assertIn(r"\texttt{f-min}", latex)
        self.assertIn(r"\slred", latex)
        self.assertIn("fastest core strategy", prose)

    def test_perseus_timing_invokes_executable_and_validates_barcode(self):
        complex_ = bench.mp.FilteredComplex.from_simplices([([0], 0.0)])

        with tempfile.TemporaryDirectory() as tmpdir:
            executable = Path(tmpdir) / "fake-perseus"
            executable.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "import sys\n"
                "prefix = Path(sys.argv[3])\n"
                "Path(f'{prefix}_0.txt').write_text('1 -1\\n')\n"
            )
            executable.chmod(0o755)

            timing = bench.time_perseus_persistence(
                complex_,
                repeats=2,
                executable=executable,
            )

        self.assertIsNotNone(timing)
        self.assertTrue(timing.matches_standard)
        self.assertGreaterEqual(timing.seconds, 0.0)
        self.assertEqual(timing.finite_barcode, ())
        self.assertEqual(timing.essential_barcode, ((0, 0.0),))

    def test_summarize_rows_groups_by_family_and_size(self):
        rows = bench.run_benchmarks(
            families=["plateau"],
            sizes=[5, 6],
            seeds=[0, 1],
            repeats=1,
            plateau_levels=3,
        )

        summaries = bench.summarize_rows(rows)

        self.assertEqual(len(summaries), 2)
        self.assertEqual([summary.family for summary in summaries], ["plateau", "plateau"])
        self.assertEqual([summary.size for summary in summaries], [5, 6])
        for summary in summaries:
            self.assertEqual(summary.cases, 2)
            self.assertGreater(summary.avg_num_simplices, 0.0)
            self.assertGreater(summary.avg_speedup, 0.0)
            self.assertGreaterEqual(summary.min_speedup, 0.0)
            self.assertGreaterEqual(summary.max_speedup, summary.min_speedup)
            measured_share = (
                summary.avg_sequence_share
                + summary.avg_reference_share
                + summary.avg_reduction_share
            )
            self.assertGreater(measured_share, 0.0)
            self.assertLessEqual(measured_share, 1.0)

    def test_summarize_rows_keeps_frame_modes_separate(self):
        fused = bench.run_benchmarks(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            plateau_levels=3,
            frame_mode="fused",
        )
        separate = bench.run_benchmarks(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            plateau_levels=3,
            frame_mode="separate",
        )

        summaries = bench.summarize_rows(fused + separate)

        self.assertEqual(len(summaries), 2)
        self.assertEqual({summary.frame_mode for summary in summaries}, {"fused", "separate"})

    def test_summary_writer_includes_aggregate_columns(self):
        rows = bench.run_benchmarks(
            families=["plateau"],
            sizes=[5],
            seeds=[0],
            repeats=1,
            plateau_levels=3,
        )
        output = StringIO()

        bench.write_summary(rows, output)

        self.assertIn("crit%", output.getvalue())
        self.assertIn("frame", output.getvalue())
        self.assertIn("pers_us", output.getvalue())
        self.assertIn("plateau", output.getvalue())

    def test_named_presets_are_available(self):
        self.assertIn("smoke", bench.PRESETS)
        self.assertIn("grid", bench.PRESETS)
        self.assertIn("scale", bench.PRESETS)
        self.assertIn("terrain", bench.PRESETS)
        self.assertIn("image-grid", bench.PRESETS)
        self.assertIn("cam", bench.PRESETS)
        self.assertIn("roadmap", bench.PRESETS)
        self.assertEqual(bench.PRESETS["smoke"].sizes, (5,))
        self.assertEqual(bench.PRESETS["terrain"].families, ("terrain",))
        self.assertEqual(bench.PRESETS["image-grid"].families, ("image-grid",))


if __name__ == "__main__":
    unittest.main()
