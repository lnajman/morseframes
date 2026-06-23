import random
import unittest
from itertools import combinations

import morseframes as mp


def edge_complex(edge_filtration=1.0):
    return mp.FilteredComplex.from_simplices(
        [
            ([0], 0.0),
            ([1], 0.0),
            ([0, 1], edge_filtration),
        ]
    )


def filled_triangle_complex():
    return mp.FilteredComplex.from_simplices(
        [
            ([0], 0.0),
            ([1], 0.0),
            ([2], 0.0),
            ([0, 1], 1.0),
            ([0, 2], 1.0),
            ([1, 2], 1.0),
            ([0, 1, 2], 2.0),
        ]
    )


def plateau_complex():
    return mp.FilteredComplex.from_simplices(
        [
            ([0], 0.0),
            ([1], 0.0),
            ([2], 0.0),
            ([3], 0.0),
            ([0, 1], 0.0),
            ([0, 2], 1.0),
            ([1, 2], 1.0),
            ([1, 3], 2.0),
            ([2, 3], 1.0),
            ([0, 1, 2], 1.0),
            ([1, 2, 3], 2.0),
        ]
    )


def priority_plateau_complex():
    return mp.FilteredComplex.from_simplices(
        [
            ([0], 1.0),
            ([1], 1.0),
            ([2], 0.0),
            ([3], 0.0),
            ([0, 1], 1.0),
            ([0, 2], 1.0),
            ([0, 3], 1.0),
            ([1, 2], 1.0),
            ([1, 3], 1.0),
            ([2, 3], 1.0),
            ([0, 1, 2], 1.0),
            ([0, 1, 3], 1.0),
            ([0, 2, 3], 1.0),
            ([1, 2, 3], 1.0),
            ([0, 1, 2, 3], 1.0),
        ]
    )


def _gudhi_available():
    try:
        import gudhi  # noqa: F401
    except ImportError:
        return False
    return True


class PythonApiTest(unittest.TestCase):
    def test_cpp_backend_is_stateful_when_available(self):
        if not mp.cpp_backend_available():
            self.skipTest("C++ backend is not built")

        complex_ = edge_complex()
        self.assertTrue(complex_.cpp_backend_active())

        sequence = mp.compute_morse_sequence(complex_)
        self.assertTrue(sequence.cpp_backend_active())
        self.assertIsNotNone(sequence._cpp_sequence)
        self.assertIsNone(sequence._cpp_payload)
        self.assertEqual(sequence.algorithm, mp.SATURATED_SEQUENCE)
        self.assertEqual(mp.compute_reference_map(complex_, sequence), ((0,), (2,), (1,)))
        self.assertEqual(mp.compute_morse_persistence(complex_, sequence).finite_barcode(), ((0, 0.0, 1.0),))

    def test_low_level_cpp_api_when_available(self):
        if not mp.cpp_backend_available():
            self.skipTest("C++ backend is not built")

        cpp_complex = mp.cpp_filtered_complex_from_simplices(
            [
                ([0], 0.0),
                ([1], 0.0),
                ([0, 1], 1.0),
            ]
        )
        self.assertEqual(cpp_complex.size, 3)
        self.assertEqual(tuple(cpp_complex.boundary(1)), (2, 0))

        cpp_sequence = mp.cpp_compute_morse_sequence(cpp_complex, algorithm="saturated")
        cpp_frame = mp.cpp_compute_morse_sequence_and_reference_map_object(
            cpp_complex,
            algorithm="saturated",
        )
        references = mp.cpp_compute_reference_map(cpp_complex, cpp_sequence)
        fused_references = mp.cpp_reference_map_to_tuple(cpp_frame.references)
        reference_object = mp.cpp_compute_reference_map_object(cpp_complex, cpp_sequence)
        field_references = mp.cpp_compute_reference_map_modp(cpp_complex, cpp_sequence, 3)
        field_coreferences = mp.cpp_compute_coreference_map_modp(cpp_complex, cpp_sequence, 3)
        self.assertEqual(
            tuple(simplex for simplex in cpp_frame.sequence.critical_simplices),
            tuple(simplex for simplex in cpp_sequence.critical_simplices),
        )
        self.assertEqual(references, ((0,), (2,), (1,)))
        self.assertEqual(field_references, (((0, 1),), ((2, 1),), ((1, 1),)))
        self.assertEqual(field_coreferences, (((0, 1),), ((2, 1),), ((1, 1),)))
        self.assertEqual(fused_references, references)
        self.assertEqual(reference_object.size, 3)
        self.assertEqual(tuple(reference_object.annotation(1)), (2,))
        self.assertEqual(mp.cpp_reference_map_to_tuple(reference_object), references)
        self.assertEqual(
            mp.cpp_compute_morse_persistence(cpp_complex, cpp_sequence).finite_barcode(),
            ((0, 0.0, 1.0),),
        )
        self.assertEqual(
            mp.cpp_compute_morse_persistence_modp(
                cpp_complex,
                cpp_sequence,
                modulus=3,
            ).finite_barcode(),
            ((0, 0.0, 1.0),),
        )
        with self.assertRaises(ValueError):
            mp.cpp_compute_morse_persistence_modp(cpp_complex, cpp_sequence, modulus=6)
        self.assertEqual(
            mp.cpp_compute_morse_coreference_persistence_modp(
                cpp_complex,
                cpp_sequence,
                modulus=3,
            ).finite_barcode(),
            ((0, 0.0, 1.0),),
        )
        with self.assertRaises(ValueError):
            mp.cpp_compute_morse_coreference_persistence_modp(
                cpp_complex,
                cpp_sequence,
                modulus=6,
            )
        self.assertEqual(
            mp.cpp_reduce_morse_persistence(cpp_complex, cpp_sequence, references).finite_barcode(),
            ((0, 0.0, 1.0),),
        )
        self.assertEqual(
            mp.cpp_reduce_morse_persistence(cpp_complex, cpp_sequence, reference_object).finite_barcode(),
            ((0, 0.0, 1.0),),
        )
        diagram, metrics = mp.cpp_reduce_morse_persistence_with_metrics(
            cpp_complex, cpp_sequence, reference_object
        )
        self.assertEqual(diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(metrics["working_set_size"], 3)
        self.assertGreaterEqual(metrics["boundary_annotation_xors"], 0)
        self.assertIn("xor_candidate_scans", metrics["inverse_store"])
        self.assertEqual(
            mp.cpp_compute_standard_persistence(cpp_complex).finite_barcode(),
            ((0, 0.0, 1.0),),
        )
        self.assertEqual(
            mp.cpp_compute_standard_persistence_modp(cpp_complex, 3).finite_barcode(),
            ((0, 0.0, 1.0),),
        )
        with self.assertRaises(ValueError):
            mp.cpp_compute_standard_persistence_modp(cpp_complex, 6)
        profile = mp.cpp_profile_morse_reference_frame(cpp_complex, algorithm="saturated")
        self.assertEqual(profile.sequence_algorithm, mp.SATURATED_SEQUENCE)
        self.assertEqual(profile.num_simplices, 3)
        self.assertGreater(profile.num_critical_simplices, 0)
        self.assertGreater(profile.reducer_working_set_size, 0)
        self.assertGreaterEqual(profile.estimated_reducer_work, 0)

    def test_cubical_grid_python_api_when_available(self):
        if not mp.cpp_backend_available() or mp.CppCubicalGrid2DComplex is None:
            self.skipTest("C++ cubical backend is not built")

        grid = mp.CubicalGrid2DComplex.from_vertex_values(
            3,
            3,
            [
                0.0,
                1.0,
                0.0,
                1.0,
                2.0,
                1.0,
                0.0,
                1.0,
                0.0,
            ],
        )
        self.assertTrue(grid.cpp_backend_active())
        self.assertEqual(grid.size, 25)
        self.assertEqual(grid.vertex_width, 3)
        self.assertEqual(grid.vertex_height, 3)
        self.assertEqual(grid.square_width, 2)
        self.assertEqual(grid.square_height, 2)
        self.assertEqual(grid.vertex(2, 1), 5)
        self.assertEqual(grid.horizontal_edge(0, 1), 11)
        self.assertEqual(grid.vertical_edge(0, 0), 15)
        self.assertEqual(grid.square(0, 0), 21)

        square = grid.square(0, 0)
        self.assertEqual(grid.cell_type(square), "square")
        self.assertEqual(grid.boundary(square), (16, 15, 11, 9))
        self.assertEqual(
            [grid.boundary_coefficient(square, index, 3) for index in range(4)],
            [1, 2, 2, 1],
        )
        self.assertEqual(grid.simplex(square).cell_type, "square")
        self.assertTrue(grid.contains(grid.vertices(square)))
        self.assertEqual(grid.filtration(square), 2.0)

        standard = mp.compute_standard_persistence(grid, modulus=3)
        algorithms = (
            "saturated",
            "same-level-reduction",
            "f-max",
            "f-min",
            "flooding-max",
            "flooding-min",
            "flooding-minmax",
            "flooding-maxmin",
        )
        for algorithm in algorithms:
            with self.subTest(algorithm=algorithm):
                sequence = mp.compute_morse_sequence(grid, algorithm=algorithm)
                self.assertTrue(sequence.cpp_backend_active())
                morse = mp.compute_morse_persistence(
                    grid,
                    sequence,
                    algorithm=algorithm,
                    modulus=3,
                )
                coreference = mp.compute_morse_coreference_persistence(
                    grid,
                    sequence,
                    algorithm=algorithm,
                    modulus=3,
                )
                self.assertEqual(morse.finite_barcode(), standard.finite_barcode())
                self.assertEqual(morse.essential_barcode(), standard.essential_barcode())
                self.assertEqual(coreference.finite_barcode(), standard.finite_barcode())
                self.assertEqual(coreference.essential_barcode(), standard.essential_barcode())

    def test_low_level_cpp_cubical_grid_when_available(self):
        if not mp.cpp_backend_available() or mp.CppCubicalGrid2DComplex is None:
            self.skipTest("C++ cubical backend is not built")

        cpp_grid = mp.cpp_cubical_grid_2d_from_vertex_values(
            2,
            2,
            [0.0, 1.0, 2.0, 3.0],
        )
        self.assertEqual(cpp_grid.size, 9)
        self.assertEqual(cpp_grid.square(0, 0), 8)
        self.assertEqual(tuple(cpp_grid.boundary(8)), (7, 6, 5, 4))
        self.assertEqual(
            [cpp_grid.boundary_coefficient(8, index, 5) for index in range(4)],
            [1, 4, 4, 1],
        )

        sequence = mp.cpp_compute_morse_sequence(cpp_grid, algorithm="saturated")
        standard = mp.cpp_compute_standard_persistence_modp(cpp_grid, 3)
        morse = mp.cpp_compute_morse_persistence_modp(cpp_grid, sequence, modulus=3)
        self.assertEqual(morse.finite_barcode(), standard.finite_barcode())
        self.assertEqual(morse.essential_barcode(), standard.essential_barcode())

    def test_morse_sequence_algorithm_option_is_explicit(self):
        complex_ = edge_complex()

        sequence = mp.compute_morse_sequence(complex_, algorithm="saturated")
        alias_sequence = mp.compute_morse_sequence(complex_, algorithm="f-sequence")
        greedy_sequence = mp.compute_morse_sequence(complex_, algorithm="plateau")
        coreduction_sequence = mp.compute_morse_sequence(complex_, algorithm="same-level-reduction")
        legacy_coreduction_sequence = mp.compute_morse_sequence(complex_, algorithm="coreduction")
        flooding_sequence = mp.compute_morse_sequence(complex_, algorithm="flooding")
        flooding_minmax_sequence = mp.compute_morse_sequence(complex_, algorithm="minmax")
        flooding_maxmin_sequence = mp.compute_morse_sequence(complex_, algorithm="flooding-maxmin")
        flooding_max_sequence = mp.compute_morse_sequence(complex_, algorithm="flooding-maximal")
        flooding_min_sequence = mp.compute_morse_sequence(complex_, algorithm="minimal-flooding")
        f_max_sequence = mp.compute_morse_sequence(complex_, algorithm="paper-max")
        f_min_sequence = mp.compute_morse_sequence(complex_, algorithm="min-s-f")
        diagram = mp.compute_morse_persistence(complex_, algorithm="saturated")
        greedy_diagram = mp.compute_morse_persistence(complex_, algorithm="plateau-greedy")
        coreduction_diagram = mp.compute_morse_persistence(complex_, algorithm="same-level-reduction")
        flooding_diagram = mp.compute_morse_persistence(complex_, algorithm="flooding-minmax")
        f_max_diagram = mp.compute_morse_persistence(complex_, algorithm="f-max")
        f_min_diagram = mp.compute_morse_persistence(complex_, algorithm="f-min")
        result = mp.compute_persistence_adaptive(complex_, sequence_algorithm="saturated")

        self.assertIn(mp.F_MAX_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.F_MIN_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.PLATEAU_GREEDY_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.SAME_LEVEL_REDUCTION_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.COREDUCTION_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.FLOODING_MAX_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.FLOODING_MIN_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.FLOODING_MINMAX_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.FLOODING_MAXMIN_SEQUENCE, mp.MORSE_SEQUENCE_ALGORITHMS)
        self.assertIn(mp.FLOODING_MAX_SEQUENCE, mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO)
        self.assertIn(mp.FLOODING_MIN_SEQUENCE, mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO)
        self.assertIn(mp.FLOODING_MINMAX_SEQUENCE, mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO)
        self.assertIn(mp.FLOODING_MAXMIN_SEQUENCE, mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO)
        self.assertIn(mp.F_MAX_SEQUENCE, mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO)
        self.assertIn(mp.F_MIN_SEQUENCE, mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO)
        self.assertEqual(sequence.algorithm, mp.SATURATED_SEQUENCE)
        self.assertEqual(alias_sequence.algorithm, mp.SATURATED_SEQUENCE)
        self.assertEqual(greedy_sequence.algorithm, mp.PLATEAU_GREEDY_SEQUENCE)
        self.assertEqual(coreduction_sequence.algorithm, mp.COREDUCTION_SEQUENCE)
        self.assertEqual(legacy_coreduction_sequence.algorithm, mp.SAME_LEVEL_REDUCTION_SEQUENCE)
        self.assertEqual(flooding_sequence.algorithm, mp.FLOODING_MINMAX_SEQUENCE)
        self.assertEqual(flooding_minmax_sequence.algorithm, mp.FLOODING_MINMAX_SEQUENCE)
        self.assertEqual(flooding_maxmin_sequence.algorithm, mp.FLOODING_MAXMIN_SEQUENCE)
        self.assertEqual(flooding_max_sequence.algorithm, mp.FLOODING_MAX_SEQUENCE)
        self.assertEqual(flooding_min_sequence.algorithm, mp.FLOODING_MIN_SEQUENCE)
        self.assertEqual(f_max_sequence.algorithm, mp.F_MAX_SEQUENCE)
        self.assertEqual(f_min_sequence.algorithm, mp.F_MIN_SEQUENCE)
        self.assertEqual(sequence.steps, alias_sequence.steps)
        self.assertEqual(diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(greedy_diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(coreduction_diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(flooding_diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(f_max_diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(f_min_diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(result.sequence.algorithm, mp.SATURATED_SEQUENCE)

        plateau = priority_plateau_complex()
        minmax_plateau = mp.compute_morse_sequence(plateau, algorithm="flooding-minmax")
        maxmin_plateau = mp.compute_morse_sequence(plateau, algorithm="flooding-maxmin")
        self.assertNotEqual(minmax_plateau.steps, maxmin_plateau.steps)
        self.assertEqual(
            mp.compute_morse_persistence(plateau, algorithm="flooding-minmax").finite_barcode(),
            mp.compute_morse_persistence(plateau, algorithm="flooding-maxmin").finite_barcode(),
        )

        if mp.cpp_backend_available() and complex_.cpp_backend_active():
            benchmark = mp.benchmark_persistence(
                complex_,
                repeats=1,
                sequence_algorithm="plateau-greedy",
                materialize_barcodes=False,
            )
            self.assertEqual(benchmark.sequence_algorithm, mp.PLATEAU_GREEDY_SEQUENCE)
            self.assertEqual(benchmark.validation_mode, "core")
            coreduction_benchmark = mp.benchmark_persistence(
                complex_,
                repeats=1,
                sequence_algorithm="same-level-reduction",
                materialize_barcodes=False,
            )
            self.assertEqual(coreduction_benchmark.sequence_algorithm, mp.COREDUCTION_SEQUENCE)
            self.assertEqual(coreduction_benchmark.validation_mode, "core")
            flooding_profile = mp.profile_morse_reference_frame(
                complex_,
                algorithm="flooding-minmax",
            )
            self.assertEqual(flooding_profile.sequence_algorithm, mp.FLOODING_MINMAX_SEQUENCE)
            self.assertGreaterEqual(flooding_profile.estimated_reducer_work, 0)

    def test_morse_sequence_algorithm_rejects_unknown_or_reserved_names(self):
        complex_ = edge_complex()

        with self.assertRaises(ValueError):
            mp.compute_morse_sequence(complex_, algorithm="not-an-algorithm")
        with self.assertRaises(NotImplementedError):
            mp.compute_morse_sequence(complex_, algorithm="stack-flooding")

    def test_fused_morse_reference_frame_matches_separate_construction(self):
        complex_ = plateau_complex()
        sequence = mp.compute_morse_sequence(complex_)
        references = mp.compute_reference_map(complex_, sequence)

        frame = mp.compute_morse_sequence_and_reference_map(complex_)

        self.assertEqual(frame.sequence.steps, sequence.steps)
        self.assertEqual(frame.references, references)
        self.assertEqual(
            mp.compute_morse_persistence(
                complex_,
                frame.sequence,
                frame.references,
            ).finite_barcode(),
            mp.compute_morse_persistence(complex_, sequence, references).finite_barcode(),
        )
        if mp.cpp_backend_available():
            self.assertTrue(frame.cpp_reference_map_active())

    def test_coreduction_coreference_frame_matches_separate_construction(self):
        complex_ = plateau_complex()
        frame = mp.compute_morse_sequence_and_coreference_map(
            complex_,
            algorithm="same-level-reduction",
        )
        separate = mp.compute_coreference_map(complex_, frame.sequence)

        self.assertEqual(frame.sequence.algorithm, mp.COREDUCTION_SEQUENCE)
        self.assertEqual(frame.coreferences, separate)
        self.assertEqual(
            mp.compute_morse_coreference_persistence(
                complex_,
                frame.sequence,
                frame.coreferences,
            ).finite_barcode(),
            mp.compute_standard_persistence(complex_).finite_barcode(),
        )
        self.assertEqual(
            mp.compute_morse_coreference_persistence(
                complex_,
                algorithm="same-level-reduction",
            ).essential_barcode(),
            mp.compute_standard_persistence(complex_).essential_barcode(),
        )
        if mp.cpp_backend_available():
            self.assertTrue(frame.cpp_coreference_map_active())

    def test_coreduction_direction_benchmark_compares_reference_and_coreference(self):
        rows = mp.benchmark_same_level_reduction_directions(plateau_complex(), repeats=1)

        self.assertEqual({row.direction for row in rows}, {"reference", "coreference"})
        for row in rows:
            self.assertGreater(row.num_simplices, 0)
            self.assertGreater(row.num_critical_simplices, 0)
            self.assertGreaterEqual(row.critical_ratio, 0.0)
            self.assertGreaterEqual(row.frame_seconds, 0.0)
            self.assertGreaterEqual(row.morse_reduction_seconds, 0.0)
            self.assertGreaterEqual(row.morse_seconds, 0.0)
            self.assertGreaterEqual(row.standard_seconds, 0.0)
            self.assertGreater(row.speedup, 0.0)
            self.assertGreaterEqual(row.reducer_working_set_size, 0)
            self.assertGreaterEqual(row.reducer_initial_total_annotation_size, 0)
            self.assertGreaterEqual(row.reducer_boundary_annotation_xors, 0)

    def test_inspection_views_materialize_simplices(self):
        complex_ = edge_complex()

        self.assertEqual(complex_.simplex_list(), ((0,), (0, 1), (1,)))
        self.assertEqual(complex_.as_simplices(), complex_.simplex_list())
        self.assertEqual(complex_.simplex_records(), tuple(complex_.simplices()))
        self.assertEqual(complex_.filtration_list(), (((0,), 0.0), ((1,), 0.0), ((0, 1), 1.0)))
        self.assertTrue(complex_.contains([1, 0]))
        self.assertIn((0, 1), complex_)
        self.assertEqual(complex_.simplex_id([1, 0]), complex_.find_simplex([0, 1]))
        self.assertEqual(complex_.filtration_of([1, 0]), 1.0)
        self.assertEqual(set(complex_.boundary_simplices([0, 1])), {(0,), (1,)})
        self.assertEqual(complex_.coboundary_simplices([0]), ((0, 1),))

        frame = mp.compute_morse_sequence_and_reference_map(complex_, algorithm="f-max")
        sequence = frame.sequence
        step_simplices = sequence.steps_as_simplices(complex_)
        self.assertEqual(step_simplices, sequence.as_simplices(complex_))
        self.assertEqual(step_simplices, mp.morse_sequence_as_simplices(complex_, sequence))
        self.assertEqual(len(step_simplices), len(sequence.steps))
        self.assertEqual(step_simplices[0].sigma, complex_.vertices(sequence.steps[0].sigma))
        self.assertEqual(
            sequence.critical_simplices_as_simplices(complex_),
            tuple(complex_.vertices(simplex) for simplex in sequence.critical_simplices),
        )
        self.assertEqual(len(sequence.paired_with_as_simplices(complex_)), complex_.size)
        for step in step_simplices:
            self.assertIn(step.sigma, complex_)
            if step.tau is not None:
                self.assertIn(step.tau, complex_)

        reference_complex = frame.reference_complex(complex_)
        self.assertEqual(reference_complex, frame.references_as_simplices(complex_))
        self.assertEqual(reference_complex, mp.reference_map_as_simplices(complex_, sequence, frame.references))
        self.assertEqual(reference_complex, mp.compute_reference_complex(complex_, sequence, frame.references))
        self.assertEqual(reference_complex[0], mp.annotation_as_simplices(complex_, sequence, frame.references[0]))
        for annotation in reference_complex:
            for simplex in annotation:
                self.assertIn(simplex, complex_)

        core_frame = mp.compute_morse_sequence_and_coreference_map(
            complex_,
            algorithm="same-level-reduction",
        )
        coreference_complex = core_frame.coreference_complex(complex_)
        self.assertEqual(coreference_complex, core_frame.coreferences_as_simplices(complex_))
        self.assertEqual(
            coreference_complex,
            mp.coreference_map_as_simplices(
                complex_,
                core_frame.sequence,
                core_frame.coreferences,
            ),
        )
        self.assertEqual(
            coreference_complex,
            mp.compute_coreference_complex(
                complex_,
                core_frame.sequence,
                core_frame.coreferences,
            ),
        )

        morse_complex = mp.compute_morse_complex(complex_, sequence, frame.references)
        self.assertEqual(morse_complex.as_simplices(complex_), sequence.critical_simplices_as_simplices(complex_))
        self.assertEqual(
            morse_complex.boundaries_as_simplices(complex_),
            {
                critical_simplex: morse_complex.boundary_as_simplices(complex_, critical_id)
                for critical_id, critical_simplex in enumerate(morse_complex.as_simplices(complex_))
            },
        )

        diagram = mp.compute_morse_persistence(complex_, sequence, frame.references)
        intervals = diagram.intervals_as_simplices(complex_)
        self.assertEqual(intervals["finite_pairs"], diagram.finite_pairs_as_simplices(complex_))
        self.assertEqual(intervals["essential"], diagram.essential_as_simplices(complex_))
        for pair in diagram.finite_pairs_as_simplices(complex_):
            self.assertIn(pair.birth, complex_)
            self.assertIn(pair.death, complex_)
        for interval in diagram.essential_as_simplices(complex_):
            self.assertIn(interval.birth, complex_)

    def test_later_edge_barcode(self):
        complex_ = edge_complex()
        sequence = mp.compute_morse_sequence(complex_)
        references = mp.compute_reference_map(complex_, sequence)
        morse_complex = mp.compute_morse_complex(complex_, sequence, references)
        diagram = mp.compute_morse_persistence(complex_, sequence, references)

        self.assertEqual(diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(diagram.essential_barcode(), ((0, 0.0),))
        self.assertGreaterEqual(len(sequence.steps), 1)
        self.assertEqual(morse_complex.boundary(0), ())
        mp.assert_matches_standard(complex_)

    def test_same_level_edge_has_no_off_diagonal_pair(self):
        complex_ = edge_complex(edge_filtration=0.0)
        diagram = mp.compute_morse_persistence(complex_)

        self.assertEqual(diagram.finite_barcode(), ())
        self.assertEqual(diagram.essential_barcode(), ((0, 0.0),))
        mp.assert_matches_standard(complex_)

    def test_filled_triangle_matches_standard(self):
        complex_ = filled_triangle_complex()
        morse = mp.compute_morse_persistence(complex_)
        standard = mp.compute_standard_persistence(complex_)

        self.assertEqual(morse.finite_barcode(), standard.finite_barcode())
        self.assertEqual(morse.essential_barcode(), standard.essential_barcode())
        self.assertIn((1, 1.0, 2.0), morse.finite_barcode())

    def test_standard_persistence_supports_prime_field_moduli(self):
        complex_ = filled_triangle_complex()
        z2 = mp.compute_standard_persistence(complex_)

        for modulus in (3, 5, 7):
            with self.subTest(modulus=modulus):
                diagram = mp.compute_standard_persistence(complex_, modulus=modulus)
                self.assertEqual(diagram.finite_barcode(), z2.finite_barcode())
                self.assertEqual(diagram.essential_barcode(), z2.essential_barcode())
                if _gudhi_available():
                    self.assertEqual(
                        (
                            diagram.finite_barcode(include_zero=True),
                            diagram.essential_barcode(),
                        ),
                        mp.gudhi_barcode(complex_, include_zero=True, modulus=modulus),
                    )

    def test_standard_persistence_modp_wrapper(self):
        complex_ = edge_complex()

        diagram = mp.compute_standard_persistence_modp(complex_, 3)

        self.assertEqual(diagram.finite_barcode(), ((0, 0.0, 1.0),))
        self.assertEqual(diagram.essential_barcode(), ((0, 0.0),))

    def test_standard_persistence_rejects_non_field_moduli(self):
        complex_ = edge_complex()

        for modulus in (0, 1, 4, 6, 9):
            with self.subTest(modulus=modulus):
                with self.assertRaises(ValueError):
                    mp.compute_standard_persistence(complex_, modulus=modulus)

    def test_morse_persistence_supports_prime_field_moduli(self):
        complex_ = plateau_complex()
        algorithms = (
            mp.SATURATED_SEQUENCE,
            mp.PLATEAU_GREEDY_SEQUENCE,
            mp.SAME_LEVEL_REDUCTION_SEQUENCE,
            mp.F_MAX_SEQUENCE,
            mp.F_MIN_SEQUENCE,
            mp.FLOODING_MAX_SEQUENCE,
            mp.FLOODING_MIN_SEQUENCE,
            mp.FLOODING_MINMAX_SEQUENCE,
            mp.FLOODING_MAXMIN_SEQUENCE,
        )

        for algorithm in algorithms:
            for modulus in (3, 5):
                with self.subTest(algorithm=algorithm, modulus=modulus):
                    morse = mp.compute_morse_persistence(
                        complex_,
                        algorithm=algorithm,
                        modulus=modulus,
                    )
                    standard = mp.compute_standard_persistence(complex_, modulus=modulus)
                    self.assertEqual(morse.finite_barcode(), standard.finite_barcode())
                    self.assertEqual(morse.essential_barcode(), standard.essential_barcode())

    def test_morse_reference_map_modp_keeps_signed_coefficients(self):
        complex_ = mp.FilteredComplex.from_lower_star(
            [(0, 1, 2), (0, 2, 3)],
            {0: 1.0, 1: 0.0, 2: 1.0, 3: 0.0},
        )

        references = mp.compute_reference_map_modp(
            complex_,
            algorithm=mp.SAME_LEVEL_REDUCTION_SEQUENCE,
            modulus=3,
        )
        morse = mp.compute_morse_persistence_modp(
            complex_,
            references=references,
            algorithm=mp.SAME_LEVEL_REDUCTION_SEQUENCE,
            modulus=3,
        )
        standard = mp.compute_standard_persistence(complex_, modulus=3)

        self.assertTrue(any(coefficient == 2 for annotation in references for _, coefficient in annotation))
        self.assertEqual(morse.finite_barcode(), standard.finite_barcode())
        self.assertEqual(morse.essential_barcode(), standard.essential_barcode())

    def test_morse_persistence_rejects_non_field_moduli(self):
        complex_ = edge_complex()

        for modulus in (0, 1, 4, 6, 9):
            with self.subTest(modulus=modulus):
                with self.assertRaises(ValueError):
                    mp.compute_morse_persistence(complex_, modulus=modulus)

    def test_morse_coreference_persistence_supports_prime_field_moduli(self):
        complex_ = plateau_complex()
        algorithms = (
            mp.SATURATED_SEQUENCE,
            mp.PLATEAU_GREEDY_SEQUENCE,
            mp.SAME_LEVEL_REDUCTION_SEQUENCE,
            mp.F_MAX_SEQUENCE,
            mp.F_MIN_SEQUENCE,
            mp.FLOODING_MAX_SEQUENCE,
            mp.FLOODING_MIN_SEQUENCE,
            mp.FLOODING_MINMAX_SEQUENCE,
            mp.FLOODING_MAXMIN_SEQUENCE,
        )

        for algorithm in algorithms:
            for modulus in (3, 5):
                with self.subTest(algorithm=algorithm, modulus=modulus):
                    morse = mp.compute_morse_coreference_persistence(
                        complex_,
                        algorithm=algorithm,
                        modulus=modulus,
                    )
                    standard = mp.compute_standard_persistence(complex_, modulus=modulus)
                    self.assertEqual(morse.finite_barcode(), standard.finite_barcode())
                    self.assertEqual(morse.essential_barcode(), standard.essential_barcode())

    def test_morse_coreference_map_modp_keeps_signed_coefficients(self):
        complex_ = mp.FilteredComplex.from_lower_star(
            [(0, 1, 2), (0, 2, 3)],
            {0: 1.0, 1: 0.0, 2: 1.0, 3: 0.0},
        )

        coreferences = mp.compute_coreference_map_modp(
            complex_,
            algorithm=mp.SATURATED_SEQUENCE,
            modulus=3,
        )
        morse = mp.compute_morse_coreference_persistence_modp(
            complex_,
            coreferences=coreferences,
            algorithm=mp.SATURATED_SEQUENCE,
            modulus=3,
        )
        standard = mp.compute_standard_persistence(complex_, modulus=3)

        self.assertTrue(any(coefficient == 2 for annotation in coreferences for _, coefficient in annotation))
        self.assertEqual(morse.finite_barcode(), standard.finite_barcode())
        self.assertEqual(morse.essential_barcode(), standard.essential_barcode())

    def test_morse_coreference_persistence_rejects_non_field_moduli(self):
        complex_ = edge_complex()

        for modulus in (0, 1, 4, 6, 9):
            with self.subTest(modulus=modulus):
                with self.assertRaises(ValueError):
                    mp.compute_morse_coreference_persistence(complex_, modulus=modulus)

    def test_plateau_complex_matches_oracles(self):
        complex_ = plateau_complex()
        sequence = mp.compute_morse_sequence(complex_)

        self.assertEqual(complex_.num_levels, 3)
        self.assertTrue(
            any(
                step.type == mp.REGULAR_PAIR
                and step.tau is not None
                and complex_.filtration(step.sigma) == complex_.filtration(step.tau)
                for step in sequence.steps
            )
        )
        mp.assert_matches_standard(complex_)
        if _gudhi_available():
            mp.assert_matches_gudhi(complex_)

    def test_adaptive_persistence_keeps_compressed_complex_on_morse_path(self):
        complex_ = mp.FilteredComplex.from_lower_star(
            [[0, 1, 2]],
            {0: 0.0, 1: 2.0, 2: 1.0},
            dimension_offset=0.0,
        )

        result = mp.compute_persistence_adaptive(complex_)

        self.assertEqual(result.method, "morse")
        self.assertIsNone(result.fallback_reason)
        self.assertLessEqual(result.critical_ratio, result.max_critical_ratio)
        self.assertEqual(result.num_simplices, complex_.size)
        self.assertEqual(result.num_critical_simplices, len(result.sequence.critical_simplices))
        self.assertGreaterEqual(result.sequence_seconds, 0.0)
        self.assertGreaterEqual(result.persistence_seconds, 0.0)
        self.assertGreaterEqual(result.total_seconds, result.sequence_seconds)
        self.assertEqual(result.finite_barcode(), mp.compute_standard_persistence(complex_).finite_barcode())
        self.assertEqual(result.essential_barcode(), mp.compute_standard_persistence(complex_).essential_barcode())

    def test_adaptive_persistence_falls_back_when_critical_ratio_is_high(self):
        complex_ = mp.FilteredComplex.from_rips_distance_matrix(
            [
                [0.0, 1.0, 2.0, 4.0],
                [1.0, 0.0, 1.0, 3.0],
                [2.0, 1.0, 0.0, 1.0],
                [4.0, 3.0, 1.0, 0.0],
            ],
            max_edge_length=2.0,
            max_dimension=2,
        )

        result = mp.compute_persistence_adaptive(complex_, max_critical_ratio=0.4)

        self.assertEqual(result.method, "standard")
        self.assertIsNotNone(result.fallback_reason)
        self.assertGreater(result.critical_ratio, result.max_critical_ratio)
        self.assertEqual(result.finite_barcode(), mp.compute_standard_persistence(complex_).finite_barcode())
        self.assertEqual(result.essential_barcode(), mp.compute_standard_persistence(complex_).essential_barcode())

    def test_adaptive_persistence_rejects_invalid_threshold(self):
        with self.assertRaises(ValueError):
            mp.compute_persistence_adaptive(edge_complex(), max_critical_ratio=-0.01)
        with self.assertRaises(ValueError):
            mp.compute_persistence_adaptive(edge_complex(), max_critical_ratio=1.01)

    def test_sequence_algorithm_portfolio_helpers(self):
        complex_ = edge_complex()

        profiles = mp.profile_morse_sequence_algorithms(
            complex_,
            algorithms=("saturated", "same-level-reduction"),
        )
        selected_profile = mp.select_morse_sequence_profile(
            complex_,
            algorithms=("saturated", "same-level-reduction"),
            selection_metric="estimated_reducer_work",
        )
        selected_total_profile = mp.select_morse_sequence_profile(
            complex_,
            algorithms=("saturated", "same-level-reduction"),
            selection_metric="profile_total_work",
        )
        candidates = mp.benchmark_sequence_algorithms(
            complex_,
            algorithms=("saturated", "same-level-reduction"),
            repeats=1,
            materialize_barcodes=False,
        )
        selected = mp.select_morse_sequence_algorithm(
            complex_,
            algorithms=("saturated", "same-level-reduction"),
            repeats=1,
            selection_metric="reducer_work",
            materialize_barcodes=False,
        )
        auto = mp.compute_persistence_adaptive(
            complex_,
            sequence_algorithm=mp.AUTO_MORSE_SEQUENCE_ALGORITHM,
            candidate_algorithms=("saturated", "same-level-reduction"),
            selection_repeats=1,
        )
        benchmark_auto = mp.compute_persistence_adaptive(
            complex_,
            sequence_algorithm=mp.AUTO_MORSE_SEQUENCE_ALGORITHM,
            candidate_algorithms=("saturated", "same-level-reduction"),
            selection_repeats=1,
            selection_mode=mp.BENCHMARK_SELECTION_MODE,
            selection_metric="morse_seconds",
        )

        self.assertEqual(len(profiles), 2)
        self.assertEqual(
            {profile.sequence_algorithm for profile in profiles},
            {mp.SATURATED_SEQUENCE, mp.COREDUCTION_SEQUENCE},
        )
        self.assertIn(selected_profile.sequence_algorithm, {profile.sequence_algorithm for profile in profiles})
        self.assertIn(
            selected_total_profile.sequence_algorithm,
            {profile.sequence_algorithm for profile in profiles},
        )
        self.assertGreaterEqual(selected_profile.estimated_reducer_work, 0)
        self.assertEqual(len(candidates), 2)
        self.assertEqual(
            {candidate.sequence_algorithm for candidate in candidates},
            {mp.SATURATED_SEQUENCE, mp.COREDUCTION_SEQUENCE},
        )
        self.assertIn(selected.sequence_algorithm, {candidate.sequence_algorithm for candidate in candidates})
        self.assertEqual(auto.selection_metric, "estimated_reducer_work")
        self.assertEqual(auto.selection_mode, mp.PROFILE_SELECTION_MODE)
        self.assertEqual(len(auto.candidate_profiles), 2)
        self.assertEqual(auto.candidate_benchmarks, ())
        self.assertIn(auto.method, {"morse", "standard"})
        self.assertIn(auto.sequence.algorithm, {mp.SATURATED_SEQUENCE, mp.COREDUCTION_SEQUENCE})
        self.assertEqual(auto.finite_barcode(), mp.compute_standard_persistence(complex_).finite_barcode())
        self.assertEqual(benchmark_auto.selection_metric, "morse_seconds")
        self.assertEqual(benchmark_auto.selection_mode, mp.BENCHMARK_SELECTION_MODE)
        self.assertEqual(len(benchmark_auto.candidate_benchmarks), 2)

    def test_invalid_nonclosed_complex_is_rejected(self):
        complex_ = mp.FilteredComplex()
        complex_.insert([0, 1], 1.0)

        with self.assertRaises(ValueError):
            complex_.finalize()

    @unittest.skipUnless(_gudhi_available(), "GUDHI is not importable")
    def test_matches_gudhi_oracle(self):
        mp.assert_matches_gudhi(edge_complex())
        mp.assert_matches_gudhi(filled_triangle_complex())

    @unittest.skipUnless(_gudhi_available(), "GUDHI is not importable")
    def test_gudhi_simplex_tree_import_export(self):
        import gudhi

        simplex_tree = gudhi.SimplexTree()
        simplex_tree.insert([0, 1], filtration=1.0)
        simplex_tree.insert([1, 2], filtration=2.0)

        complex_ = mp.FilteredComplex.from_gudhi_simplex_tree(simplex_tree)

        self.assertEqual(complex_.size, simplex_tree.num_simplices())
        self.assertIsNotNone(complex_.find_simplex([0, 1]))
        mp.assert_matches_gudhi(complex_)

        exported = (
            mp.SimplexTreeBuilder()
            .insert([0, 1, 2], 2.0)
            .to_gudhi_simplex_tree()
        )

        self.assertEqual(exported.num_simplices(), 7)
        exported_simplices = {tuple(simplex) for simplex, _ in exported.get_filtration()}
        self.assertIn((0, 1, 2), exported_simplices)

    def test_facets_constructor_closes_faces(self):
        complex_ = mp.FilteredComplex.from_facets(
            [[0, 1, 2]],
            simplex_filtration=lambda simplex: float(len(simplex) - 1),
        )

        self.assertEqual(complex_.size, 7)
        self.assertIsNotNone(complex_.find_simplex([0, 1]))
        self.assertEqual(complex_.filtration(complex_.find_simplex([0, 1, 2])), 2.0)
        mp.assert_matches_standard(complex_)

    def test_lower_star_constructor(self):
        complex_ = mp.FilteredComplex.from_lower_star(
            [[0, 1, 2]],
            {0: 0.0, 1: 2.0, 2: 1.0},
            dimension_offset=0.25,
        )

        self.assertEqual(complex_.filtration(complex_.find_simplex([0, 2])), 1.25)
        self.assertEqual(complex_.filtration(complex_.find_simplex([0, 1, 2])), 2.5)
        mp.assert_matches_standard(complex_)

    def test_simplex_tree_builder_finalize_matches_facets(self):
        builder = mp.SimplexTreeBuilder()
        self.assertIs(
            builder.insert_facet(
                [0, 1, 2],
                simplex_filtration=lambda simplex: float(len(simplex) - 1),
            ),
            builder,
        )

        self.assertEqual(builder.size, 7)
        self.assertEqual(builder.max_dimension, 2)
        self.assertTrue(builder.contains([0, 1]))
        self.assertIn([0, 1, 2], builder)
        self.assertEqual(builder.filtration([0, 1, 2]), 2.0)

        complex_ = builder.finalize()

        self.assertEqual(len(builder), 0)
        self.assertEqual(complex_.size, 7)
        self.assertIsNotNone(complex_.find_simplex([0, 1]))
        self.assertEqual(complex_.filtration(complex_.find_simplex([0, 1, 2])), 2.0)
        mp.assert_matches_standard(complex_)

    def test_simplex_tree_builder_can_be_kept_after_finalize(self):
        builder = mp.SimplexTreeBuilder()
        builder.insert([0, 1], 1.0, include_faces=True)

        complex_ = builder.finalize(clear=False)

        self.assertEqual(len(builder), 3)
        self.assertTrue(builder.contains([0, 1]))
        self.assertIsNotNone(complex_.find_simplex([0, 1]))
        mp.assert_matches_standard(complex_)

    def test_simplex_tree_builder_insert_closes_faces_by_default(self):
        builder = mp.SimplexTreeBuilder()
        builder.insert([0, 1], 1.0)

        self.assertEqual(builder.num_simplices(), 3)
        self.assertEqual(builder.num_vertices(), 2)
        self.assertTrue(builder.find_simplex([0]))
        self.assertTrue(builder.find_simplex([1]))
        self.assertTrue(builder.find_simplex([0, 1]))
        self.assertEqual(builder.simplex_filtration([0, 1]), 1.0)

        complex_ = builder.finalize()
        self.assertEqual(complex_.size, 3)
        mp.assert_matches_standard(complex_)

    def test_simplex_tree_builder_min_merge_lowers_repeated_faces(self):
        builder = mp.SimplexTreeBuilder()
        builder.insert([0], 5.0)
        builder.insert([0, 1], 1.0, include_faces=True)

        self.assertEqual(builder.filtration([0]), 1.0)
        self.assertEqual(builder.filtration([1]), 1.0)
        self.assertEqual(builder.filtration([0, 1]), 1.0)

        complex_ = builder.finalize()
        mp.assert_matches_standard(complex_)

    def test_simplex_tree_builder_strict_duplicate_rejects_conflict(self):
        builder = mp.FilteredComplexBuilder(merge="strict")
        builder.insert([0], 1.0)

        with self.assertRaises(ValueError):
            builder.insert([0], 2.0)

    def test_simplex_tree_builder_lower_star_facet(self):
        complex_ = (
            mp.SimplexTreeBuilder.from_lower_star(
                [[0, 1, 2]],
                {0: 0.0, 1: 2.0, 2: 1.0},
                dimension_offset=0.25,
            )
            .finalize()
        )

        self.assertEqual(complex_.filtration(complex_.find_simplex([0, 2])), 1.25)
        self.assertEqual(complex_.filtration(complex_.find_simplex([0, 1, 2])), 2.5)
        mp.assert_matches_standard(complex_)

    def test_graph_constructor(self):
        complex_ = mp.FilteredComplex.from_graph(
            [(0, 1), (1, 2)],
            vertices=[3],
            vertex_filtration={0: 0.0, 1: 0.5, 2: 1.0, 3: 3.0},
            edge_filtration={(0, 1): 1.5, (1, 2): 2.0},
        )

        self.assertEqual(complex_.size, 6)
        self.assertIsNotNone(complex_.find_simplex([3]))
        self.assertIsNone(complex_.find_simplex([0, 2]))
        mp.assert_matches_standard(complex_)

    def test_rips_distance_matrix_constructor(self):
        complex_ = mp.FilteredComplex.from_rips_distance_matrix(
            [
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 3.0],
                [1.0, 3.0, 0.0],
            ],
            max_edge_length=1.1,
            max_dimension=2,
        )

        self.assertIsNotNone(complex_.find_simplex([0, 1]))
        self.assertIsNotNone(complex_.find_simplex([0, 2]))
        self.assertIsNone(complex_.find_simplex([1, 2]))
        self.assertIsNone(complex_.find_simplex([0, 1, 2]))
        self.assertEqual(mp.compute_morse_persistence(complex_).finite_barcode(), ((0, 0.0, 1.0), (0, 0.0, 1.0)))

    def test_random_lower_star_complexes_match_oracles(self):
        for seed in range(20):
            with self.subTest(seed=seed):
                rng = random.Random(seed)
                n = rng.randint(3, 7)
                facets = [[vertex] for vertex in range(n)]
                for _ in range(rng.randint(1, 6)):
                    size = rng.randint(2, min(4, n))
                    facets.append(rng.sample(range(n), size))
                vertex_values = {vertex: float(rng.randint(0, 5)) for vertex in range(n)}
                complex_ = mp.FilteredComplex.from_lower_star(
                    facets,
                    vertex_values,
                    dimension_offset=0.125,
                )

                mp.assert_matches_standard(complex_)
                if _gudhi_available():
                    mp.assert_matches_gudhi(complex_)

    def test_random_graphs_match_oracles(self):
        for seed in range(20, 40):
            with self.subTest(seed=seed):
                rng = random.Random(seed)
                n = rng.randint(3, 8)
                vertices = list(range(n))
                edges = []
                edge_values = {}
                vertex_values = {vertex: float(rng.randint(0, 3)) for vertex in vertices}
                for i in range(n):
                    for j in range(i + 1, n):
                        if rng.random() < 0.35:
                            edge = (i, j)
                            edges.append(edge)
                            edge_values[edge] = max(vertex_values[i], vertex_values[j]) + float(rng.randint(0, 4))

                complex_ = mp.FilteredComplex.from_graph(
                    edges,
                    vertices=vertices,
                    vertex_filtration=vertex_values,
                    edge_filtration=edge_values,
                )

                mp.assert_matches_standard(complex_)
                if _gudhi_available():
                    mp.assert_matches_gudhi(complex_)

    def test_random_rips_complexes_match_oracles(self):
        for seed in range(40, 50):
            with self.subTest(seed=seed):
                rng = random.Random(seed)
                n = rng.randint(3, 6)
                distances = [[0.0 for _ in range(n)] for _ in range(n)]
                for i in range(n):
                    for j in range(i + 1, n):
                        value = float(rng.randint(1, 5))
                        distances[i][j] = value
                        distances[j][i] = value

                complex_ = mp.FilteredComplex.from_rips_distance_matrix(
                    distances,
                    max_edge_length=3.0,
                    max_dimension=2,
                )

                mp.assert_matches_standard(complex_)
                if _gudhi_available():
                    mp.assert_matches_gudhi(complex_)

    def test_random_plateau_complexes_match_oracles(self):
        for seed in range(50, 60):
            with self.subTest(seed=seed):
                rng = random.Random(seed)
                n = rng.randint(3, 7)
                levels = 3
                facets = [[vertex] for vertex in range(n)]
                for _ in range(rng.randint(1, 6)):
                    size = rng.randint(2, min(4, n))
                    facets.append(rng.sample(range(n), size))

                simplices = set()
                for facet in facets:
                    vertices = tuple(sorted(facet))
                    for dimension in range(1, len(vertices) + 1):
                        simplices.update(combinations(vertices, dimension))

                filtration = {}
                for simplex in sorted(simplices, key=lambda simplex: (len(simplex), simplex)):
                    if len(simplex) == 1:
                        lower_bound = 0
                    else:
                        lower_bound = max(
                            int(filtration[face])
                            for face in combinations(simplex, len(simplex) - 1)
                        )
                    filtration[simplex] = float(max(lower_bound, rng.randrange(levels)))

                complex_ = mp.FilteredComplex.from_simplices(filtration.items())

                mp.assert_matches_standard(complex_)
                if _gudhi_available():
                    mp.assert_matches_gudhi(complex_)

    def test_benchmark_persistence_lower_star(self):
        complex_ = mp.FilteredComplex.from_lower_star(
            [[0, 1, 2, 3], [1, 3, 4]],
            {0: 0.0, 1: 2.0, 2: 1.0, 3: 3.0, 4: 1.5},
            dimension_offset=0.1,
        )

        result = mp.benchmark_persistence(complex_, repeats=2)
        self.assertEqual(result.num_simplices, complex_.size)
        self.assertEqual(result.num_levels, complex_.num_levels)
        self.assertGreater(result.num_critical_simplices, 0)
        self.assertEqual(result.sequence_algorithm, mp.SATURATED_SEQUENCE)
        self.assertEqual(result.frame_mode, mp.FUSED_FRAME)
        self.assertTrue(result.barcodes_materialized)
        self.assertEqual(result.validation_mode, "materialized")
        self.assertGreaterEqual(result.reducer_working_set_size, 0)
        self.assertEqual(result.reference_final_live_nonempty_annotations, 0)
        self.assertEqual(result.reference_final_live_total_annotation_size, 0)
        self.assertEqual(result.reference_peak_live_nonempty_annotations, 0)
        self.assertEqual(result.reference_peak_live_total_annotation_size, 0)
        self.assertEqual(result.reference_released_annotations, 0)
        self.assertEqual(result.reference_released_total_annotation_size, 0)
        self.assertGreaterEqual(result.reducer_initial_nonempty_annotations, 0)
        self.assertGreaterEqual(result.reducer_initial_total_annotation_size, 0)
        self.assertGreaterEqual(result.reducer_initial_max_annotation_size, 0)
        self.assertGreaterEqual(result.reducer_initial_inverse_list_entries, 0)
        self.assertEqual(
            result.reducer_initial_total_annotation_size,
            result.reducer_initial_inverse_list_entries,
        )
        self.assertGreaterEqual(result.reducer_boundary_plan_face_scans, 0)
        self.assertGreaterEqual(result.reducer_boundary_annotation_candidate_criticals, 0)
        self.assertGreaterEqual(result.reducer_boundary_annotation_zero_skipped_criticals, 0)
        self.assertGreaterEqual(result.reducer_boundary_annotation_zero_skipped_faces, 0)
        self.assertGreaterEqual(result.reducer_boundary_annotation_xors, 0)
        self.assertGreaterEqual(result.reducer_boundary_annotation_total_output_size, 0)
        self.assertGreaterEqual(result.reducer_boundary_annotation_max_output_size, 0)
        self.assertGreaterEqual(result.reducer_pivot_eliminations, 0)
        self.assertGreaterEqual(result.reducer_remove_candidate_scans, 0)
        self.assertGreaterEqual(result.reducer_remove_total_annotation_size, 0)
        self.assertGreaterEqual(result.reducer_xor_candidate_scans, 0)
        self.assertGreaterEqual(result.reducer_xor_applied, 0)
        self.assertGreaterEqual(result.reducer_xor_total_input_size, 0)
        self.assertGreaterEqual(result.reducer_xor_total_output_size, 0)
        self.assertGreaterEqual(result.reducer_xor_inserted_labels, 0)
        self.assertGreaterEqual(result.reducer_xor_removed_labels, 0)
        self.assertGreaterEqual(result.sequence_seconds, 0.0)
        self.assertGreaterEqual(result.reference_seconds, 0.0)
        self.assertGreaterEqual(result.morse_reduction_seconds, 0.0)
        self.assertGreaterEqual(result.reducer_setup_seconds, 0.0)
        self.assertGreaterEqual(result.reducer_compute_seconds, 0.0)
        self.assertGreaterEqual(result.morse_seconds, 0.0)
        self.assertGreaterEqual(result.standard_seconds, 0.0)
        self.assertAlmostEqual(
            result.morse_seconds,
            result.sequence_seconds + result.reference_seconds + result.morse_reduction_seconds,
            places=6,
        )
        self.assertEqual(result.finite_barcode, mp.compute_standard_persistence(complex_).finite_barcode())
        self.assertEqual(result.finite_interval_count, len(result.finite_barcode))
        self.assertEqual(result.essential_interval_count, len(result.essential_barcode))

        if mp.cpp_backend_available() and complex_.cpp_backend_active():
            core_result = mp.benchmark_persistence(
                complex_,
                repeats=1,
                materialize_barcodes=False,
            )
            self.assertFalse(core_result.barcodes_materialized)
            self.assertEqual(core_result.validation_mode, "core")
            self.assertEqual(core_result.finite_barcode, ())
            self.assertEqual(core_result.essential_barcode, ())
            self.assertEqual(core_result.finite_interval_count, len(result.finite_barcode))
            self.assertEqual(core_result.essential_interval_count, len(result.essential_barcode))
            self.assertEqual(core_result.reference_final_live_nonempty_annotations, 0)
            self.assertEqual(core_result.reference_final_live_total_annotation_size, 0)
            self.assertGreaterEqual(core_result.reference_peak_live_nonempty_annotations, 0)
            self.assertGreaterEqual(core_result.reference_peak_live_total_annotation_size, 0)
            self.assertGreaterEqual(core_result.reference_released_annotations, 0)
            self.assertGreaterEqual(core_result.reference_released_total_annotation_size, 0)

    def test_benchmark_persistence_can_force_frame_modes(self):
        complex_ = plateau_complex()

        fused = mp.benchmark_persistence(complex_, repeats=1, frame_mode="fused")
        separate = mp.benchmark_persistence(complex_, repeats=1, frame_mode="separate")

        self.assertEqual(fused.frame_mode, mp.FUSED_FRAME)
        self.assertEqual(separate.frame_mode, mp.SEPARATE_FRAME)
        self.assertEqual(fused.finite_barcode, separate.finite_barcode)
        self.assertEqual(fused.essential_barcode, separate.essential_barcode)
        self.assertEqual(fused.reference_seconds, 0.0)
        self.assertGreaterEqual(separate.reference_seconds, 0.0)

        with self.assertRaises(ValueError):
            mp.benchmark_persistence(complex_, frame_mode="not-a-mode")

    def test_benchmark_persistence_rips(self):
        complex_ = mp.FilteredComplex.from_rips_distance_matrix(
            [
                [0.0, 1.0, 2.0, 4.0],
                [1.0, 0.0, 1.0, 3.0],
                [2.0, 1.0, 0.0, 1.0],
                [4.0, 3.0, 1.0, 0.0],
            ],
            max_edge_length=2.0,
            max_dimension=2,
        )

        result = mp.benchmark_persistence(complex_, repeats=2)
        self.assertEqual(result.num_simplices, complex_.size)
        self.assertGreater(result.speedup, 0.0)
        self.assertEqual(result.essential_barcode, mp.compute_standard_persistence(complex_).essential_barcode())


if __name__ == "__main__":
    unittest.main()
