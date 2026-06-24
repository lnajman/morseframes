"""Python interface for MorseFrames.

This module intentionally mirrors the C++ core in ``include/morseframes``.
When the nanobind extension is available, finalized complexes and Morse
sequences carry stateful C++ objects; otherwise the same API falls back to
the pure-Python implementation.
"""

from __future__ import annotations

import heapq
import os
from collections import deque
from dataclasses import dataclass, field
from itertools import combinations
from math import inf, isclose
from time import perf_counter
from typing import Callable, Iterable, Iterator, Mapping, Sequence

__version__ = "0.1.0a1"

if os.environ.get("MORSEFRAMES_DISABLE_CPP_BACKEND"):
    _morse_core = None
else:
    try:
        from . import _morse_core
    except ImportError:
        _morse_core = None

CppFilteredComplex = getattr(_morse_core, "FilteredComplex", None)
CppCubicalGrid2DComplex = getattr(_morse_core, "CubicalGrid2DComplex", None)
CppMorseSequence = getattr(_morse_core, "MorseSequence", None)
CppReferenceMap = getattr(_morse_core, "ReferenceMap", None)
CppMorseReferenceFrame = getattr(_morse_core, "MorseReferenceFrame", None)
CppMorseCoreferenceFrame = getattr(_morse_core, "MorseCoreferenceFrame", None)
CppSimplexTreeBuilder = getattr(_morse_core, "SimplexTreeBuilder", None)

Simplex = tuple[int, ...]
Annotation = tuple[int, ...]
FieldAnnotation = tuple[tuple[int, int], ...]
FiltrationSpec = float | Mapping[Simplex, float] | Callable[[Simplex], float]
VertexFiltrationSpec = float | Mapping[int, float] | Sequence[float] | Callable[[int], float]
EdgeFiltrationSpec = float | Mapping[tuple[int, int], float] | Callable[[tuple[int, int]], float]

CRITICAL = "critical"
REGULAR_PAIR = "regular_pair"

# Strategy taxonomy: saturated/same-level/plateau-greedy/flooding-* are
# filtration-monotone flooding constructions; f-max and f-min are global
# seed-and-expand F-sequence builders that still keep every regular pair inside
# one filtration level.
SATURATED_SEQUENCE = "saturated"
PLATEAU_GREEDY_SEQUENCE = "plateau-greedy"
SAME_LEVEL_REDUCTION_SEQUENCE = "same-level-reduction"
COREDUCTION_SEQUENCE = SAME_LEVEL_REDUCTION_SEQUENCE
F_MAX_SEQUENCE = "f-max"
F_MIN_SEQUENCE = "f-min"
FLOODING_MAX_SEQUENCE = "flooding-max"
FLOODING_MIN_SEQUENCE = "flooding-min"
FLOODING_MINMAX_SEQUENCE = "flooding-minmax"
FLOODING_MAXMIN_SEQUENCE = "flooding-maxmin"
AUTO_MORSE_SEQUENCE_ALGORITHM = "auto"
DEFAULT_MORSE_SEQUENCE_ALGORITHM = SATURATED_SEQUENCE
MORSE_SEQUENCE_ALGORITHMS = (
    SATURATED_SEQUENCE,
    PLATEAU_GREEDY_SEQUENCE,
    SAME_LEVEL_REDUCTION_SEQUENCE,
    F_MAX_SEQUENCE,
    F_MIN_SEQUENCE,
    FLOODING_MAX_SEQUENCE,
    FLOODING_MIN_SEQUENCE,
    FLOODING_MINMAX_SEQUENCE,
    FLOODING_MAXMIN_SEQUENCE,
)
DEFAULT_MORSE_ALGORITHM_PORTFOLIO = (
    SATURATED_SEQUENCE,
    F_MAX_SEQUENCE,
    F_MIN_SEQUENCE,
    SAME_LEVEL_REDUCTION_SEQUENCE,
    FLOODING_MAX_SEQUENCE,
    FLOODING_MIN_SEQUENCE,
    FLOODING_MINMAX_SEQUENCE,
    FLOODING_MAXMIN_SEQUENCE,
    PLATEAU_GREEDY_SEQUENCE,
)
MORSE_ALGORITHM_SELECTION_METRICS = ("morse_seconds", "critical_ratio", "reducer_work")
MORSE_PROFILE_SELECTION_METRICS = (
    "estimated_reducer_work",
    "profile_total_work",
    "critical_ratio",
    "profile_seconds",
    "working_set_size",
    "initial_annotation_size",
    "boundary_annotation_work",
)
PROFILE_SELECTION_MODE = "profile"
BENCHMARK_SELECTION_MODE = "benchmark"
MORSE_ALGORITHM_SELECTION_MODES = (PROFILE_SELECTION_MODE, BENCHMARK_SELECTION_MODE)
RESERVED_MORSE_SEQUENCE_ALGORITHMS = ("stack-flooding",)
FUSED_FRAME = "fused"
SEPARATE_FRAME = "separate"
DEFAULT_MORSE_FRAME_MODE = FUSED_FRAME
MORSE_FRAME_MODES = (FUSED_FRAME, SEPARATE_FRAME)


@dataclass(frozen=True)
class SimplexRecord:
    id: int
    vertices: Simplex
    dimension: int
    level: int
    filtration: float
    boundary: tuple[int, ...]
    coboundary: tuple[int, ...]


@dataclass(frozen=True)
class CubicalCellRecord:
    id: int
    vertices: tuple[int, ...]
    dimension: int
    level: int
    filtration: float
    boundary: tuple[int, ...]
    coboundary: tuple[int, ...]
    cell_type: str


@dataclass(frozen=True)
class MorseStep:
    type: str
    sigma: int
    tau: int | None
    level: int


@dataclass(frozen=True)
class MorseStepSimplices:
    type: str
    sigma: Simplex
    tau: Simplex | None
    level: int


@dataclass(frozen=True)
class PersistencePair:
    birth: int
    death: int
    dimension: int
    birth_value: float
    death_value: float


@dataclass(frozen=True)
class PersistencePairSimplices:
    birth: Simplex
    death: Simplex
    dimension: int
    birth_value: float
    death_value: float


@dataclass(frozen=True)
class EssentialInterval:
    birth: int
    dimension: int
    birth_value: float


@dataclass(frozen=True)
class EssentialIntervalSimplices:
    birth: Simplex
    dimension: int
    birth_value: float


@dataclass(frozen=True)
class PersistenceDiagram:
    finite_pairs: tuple[PersistencePair, ...]
    essential: tuple[EssentialInterval, ...]

    def finite_barcode(self, *, include_zero: bool = False) -> tuple[tuple[int, float, float], ...]:
        pairs = (
            (pair.dimension, pair.birth_value, pair.death_value)
            for pair in self.finite_pairs
            if include_zero or pair.birth_value < pair.death_value
        )
        return tuple(sorted(pairs))

    def essential_barcode(self) -> tuple[tuple[int, float], ...]:
        return tuple(sorted((bar.dimension, bar.birth_value) for bar in self.essential))

    def finite_pairs_as_simplices(self, complex_: "FilteredComplex") -> tuple[PersistencePairSimplices, ...]:
        return tuple(
            PersistencePairSimplices(
                birth=complex_.vertices(pair.birth),
                death=complex_.vertices(pair.death),
                dimension=pair.dimension,
                birth_value=pair.birth_value,
                death_value=pair.death_value,
            )
            for pair in self.finite_pairs
        )

    def essential_as_simplices(self, complex_: "FilteredComplex") -> tuple[EssentialIntervalSimplices, ...]:
        return tuple(
            EssentialIntervalSimplices(
                birth=complex_.vertices(interval.birth),
                dimension=interval.dimension,
                birth_value=interval.birth_value,
            )
            for interval in self.essential
        )

    def intervals_as_simplices(
        self,
        complex_: "FilteredComplex",
    ) -> dict[str, tuple[PersistencePairSimplices, ...] | tuple[EssentialIntervalSimplices, ...]]:
        return {
            "finite_pairs": self.finite_pairs_as_simplices(complex_),
            "essential": self.essential_as_simplices(complex_),
        }


@dataclass(frozen=True)
class PersistenceBenchmark:
    num_simplices: int
    num_levels: int
    num_critical_simplices: int
    sequence_algorithm: str
    frame_mode: str
    reference_final_live_nonempty_annotations: int
    reference_final_live_total_annotation_size: int
    reference_peak_live_nonempty_annotations: int
    reference_peak_live_total_annotation_size: int
    reference_released_annotations: int
    reference_released_total_annotation_size: int
    reducer_working_set_size: int
    reducer_initial_nonempty_annotations: int
    reducer_initial_total_annotation_size: int
    reducer_initial_max_annotation_size: int
    reducer_initial_inverse_list_entries: int
    reducer_boundary_plan_face_scans: int
    reducer_boundary_annotation_candidate_criticals: int
    reducer_boundary_annotation_zero_skipped_criticals: int
    reducer_boundary_annotation_zero_skipped_faces: int
    reducer_boundary_annotation_xors: int
    reducer_boundary_annotation_total_input_size: int
    reducer_boundary_annotation_total_output_size: int
    reducer_boundary_annotation_max_size: int
    reducer_boundary_annotation_max_output_size: int
    reducer_pivot_eliminations: int
    reducer_remove_candidate_scans: int
    reducer_remove_applied: int
    reducer_remove_total_annotation_size: int
    reducer_remove_max_annotation_size: int
    reducer_xor_candidate_scans: int
    reducer_xor_applied: int
    reducer_xor_changed_labels: int
    reducer_xor_total_input_size: int
    reducer_xor_total_output_size: int
    reducer_xor_max_input_size: int
    reducer_xor_max_output_size: int
    reducer_xor_inserted_labels: int
    reducer_xor_removed_labels: int
    reducer_inverse_list_appends: int
    sequence_seconds: float
    reference_seconds: float
    morse_reduction_seconds: float
    reducer_setup_seconds: float
    reducer_compute_seconds: float
    morse_seconds: float
    standard_seconds: float
    repeats: int
    finite_interval_count: int
    essential_interval_count: int
    barcodes_materialized: bool
    validation_mode: str
    finite_barcode: tuple[tuple[int, float, float], ...]
    essential_barcode: tuple[tuple[int, float], ...]

    @property
    def speedup(self) -> float:
        if self.morse_seconds == 0.0:
            return inf
        return self.standard_seconds / self.morse_seconds


@dataclass(frozen=True)
class MorseReferenceProfile:
    num_simplices: int
    num_levels: int
    num_critical_simplices: int
    num_regular_pairs: int
    sequence_algorithm: str
    frame_mode: str
    reference_final_live_nonempty_annotations: int
    reference_final_live_total_annotation_size: int
    reference_peak_live_nonempty_annotations: int
    reference_peak_live_total_annotation_size: int
    reference_released_annotations: int
    reference_released_total_annotation_size: int
    reducer_working_set_size: int
    reducer_initial_nonempty_annotations: int
    reducer_initial_total_annotation_size: int
    reducer_initial_max_annotation_size: int
    reducer_initial_inverse_list_entries: int
    reducer_boundary_plan_face_scans: int
    reducer_boundary_annotation_candidate_criticals: int
    reducer_boundary_annotation_zero_skipped_criticals: int
    reducer_boundary_annotation_zero_skipped_faces: int
    reducer_boundary_annotation_xors: int
    reducer_boundary_annotation_total_input_size: int
    reducer_boundary_annotation_total_output_size: int
    reducer_boundary_annotation_max_size: int
    reducer_boundary_annotation_max_output_size: int
    estimated_reducer_work: int
    profile_seconds: float
    metrics: dict[str, object] = field(default_factory=dict, repr=False, compare=False)

    @property
    def critical_ratio(self) -> float:
        if self.num_simplices == 0:
            return 0.0
        return float(self.num_critical_simplices) / float(self.num_simplices)

    @property
    def boundary_annotation_work(self) -> int:
        return (
            self.reducer_boundary_annotation_xors
            + self.reducer_boundary_annotation_total_input_size
            + self.reducer_boundary_annotation_total_output_size
        )


@dataclass(frozen=True)
class CoreductionDirectionBenchmark:
    direction: str
    num_simplices: int
    num_levels: int
    num_critical_simplices: int
    frame_seconds: float
    morse_reduction_seconds: float
    reducer_setup_seconds: float
    reducer_compute_seconds: float
    morse_seconds: float
    standard_seconds: float
    repeats: int
    finite_interval_count: int
    essential_interval_count: int
    reducer_working_set_size: int
    reducer_initial_total_annotation_size: int
    reducer_boundary_annotation_xors: int
    reducer_boundary_annotation_total_input_size: int
    reducer_boundary_annotation_total_output_size: int
    reducer_pivot_eliminations: int
    reducer_xor_applied: int
    reducer_xor_total_input_size: int
    reducer_xor_total_output_size: int
    metrics: dict[str, object] = field(default_factory=dict, repr=False, compare=False)

    @property
    def critical_ratio(self) -> float:
        if self.num_simplices == 0:
            return 0.0
        return float(self.num_critical_simplices) / float(self.num_simplices)

    @property
    def speedup(self) -> float:
        if self.morse_seconds == 0.0:
            return inf
        return self.standard_seconds / self.morse_seconds


@dataclass(frozen=True)
class AdaptivePersistenceResult:
    diagram: PersistenceDiagram
    method: str
    sequence: MorseSequence
    num_simplices: int
    num_critical_simplices: int
    critical_ratio: float
    max_critical_ratio: float
    sequence_seconds: float
    persistence_seconds: float
    total_seconds: float
    fallback_reason: str | None = None
    candidate_profiles: tuple[MorseReferenceProfile, ...] = ()
    candidate_benchmarks: tuple[PersistenceBenchmark, ...] = ()
    selection_metric: str | None = None
    selection_mode: str | None = None

    def finite_barcode(self, *, include_zero: bool = False) -> tuple[tuple[int, float, float], ...]:
        return self.diagram.finite_barcode(include_zero=include_zero)

    def essential_barcode(self) -> tuple[tuple[int, float], ...]:
        return self.diagram.essential_barcode()


@dataclass(frozen=True)
class MorseSequence:
    steps: tuple[MorseStep, ...]
    critical_simplices: tuple[int, ...]
    critical_index_of_simplex: tuple[int, ...]
    paired_with: tuple[int | None, ...]
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM
    _cpp_payload: dict[str, object] | None = field(default=None, repr=False, compare=False)
    _cpp_sequence: object | None = field(default=None, repr=False, compare=False)
    _complex_signature: tuple[tuple[Simplex, float], ...] | None = field(
        default=None, repr=False, compare=False
    )

    def critical_index(self, simplex_id: int) -> int:
        return self.critical_index_of_simplex[simplex_id]

    def is_critical(self, simplex_id: int) -> bool:
        return self.critical_index(simplex_id) >= 0

    def steps_as_simplices(self, complex_: "FilteredComplex") -> tuple[MorseStepSimplices, ...]:
        return tuple(
            MorseStepSimplices(
                type=step.type,
                sigma=complex_.vertices(step.sigma),
                tau=None if step.tau is None else complex_.vertices(step.tau),
                level=step.level,
            )
            for step in self.steps
        )

    def critical_simplices_as_simplices(self, complex_: "FilteredComplex") -> tuple[Simplex, ...]:
        return tuple(complex_.vertices(simplex) for simplex in self.critical_simplices)

    def pairs_as_simplices(self, complex_: "FilteredComplex") -> tuple[tuple[Simplex, Simplex], ...]:
        return tuple(
            (complex_.vertices(step.sigma), complex_.vertices(step.tau))
            for step in self.steps
            if step.type == REGULAR_PAIR and step.tau is not None
        )

    def paired_with_as_simplices(self, complex_: "FilteredComplex") -> tuple[Simplex | None, ...]:
        return tuple(
            None if simplex is None else complex_.vertices(simplex)
            for simplex in self.paired_with
        )

    def as_simplices(self, complex_: "FilteredComplex") -> tuple[MorseStepSimplices, ...]:
        return self.steps_as_simplices(complex_)

    def cpp_backend_active(self) -> bool:
        return self._cpp_sequence is not None


@dataclass
class MorseReferenceFrame:
    sequence: MorseSequence
    _references: tuple[Annotation, ...] | None = field(default=None, repr=False)
    _cpp_references: object | None = field(default=None, repr=False)
    _cpp_frame: object | None = field(default=None, repr=False)

    @property
    def references(self) -> tuple[Annotation, ...]:
        if self._references is None:
            cpp_references = self._cpp_references
            if cpp_references is None and self._cpp_frame is not None:
                cpp_references = self._cpp_frame.references
                self._cpp_references = cpp_references
            if cpp_references is None:
                raise RuntimeError("Morse reference frame does not contain references.")
            self._references = cpp_reference_map_to_tuple(cpp_references)
        return self._references

    def cpp_reference_map_active(self) -> bool:
        return self._cpp_references is not None or self._cpp_frame is not None

    def references_as_simplices(self, complex_: "FilteredComplex") -> tuple[tuple[Simplex, ...], ...]:
        return reference_map_as_simplices(complex_, self.sequence, self.references)

    def reference_complex(self, complex_: "FilteredComplex") -> tuple[tuple[Simplex, ...], ...]:
        return self.references_as_simplices(complex_)


@dataclass
class MorseCoreferenceFrame:
    sequence: MorseSequence
    _coreferences: tuple[Annotation, ...] | None = field(default=None, repr=False)
    _cpp_coreferences: object | None = field(default=None, repr=False)
    _cpp_frame: object | None = field(default=None, repr=False)

    @property
    def coreferences(self) -> tuple[Annotation, ...]:
        if self._coreferences is None:
            cpp_coreferences = self._cpp_coreferences
            if cpp_coreferences is None and self._cpp_frame is not None:
                cpp_coreferences = self._cpp_frame.coreferences
                self._cpp_coreferences = cpp_coreferences
            if cpp_coreferences is None:
                raise RuntimeError("Morse coreference frame does not contain coreferences.")
            self._coreferences = cpp_reference_map_to_tuple(cpp_coreferences)
        return self._coreferences

    def cpp_coreference_map_active(self) -> bool:
        return self._cpp_coreferences is not None or self._cpp_frame is not None

    def coreferences_as_simplices(self, complex_: "FilteredComplex") -> tuple[tuple[Simplex, ...], ...]:
        return reference_map_as_simplices(complex_, self.sequence, self.coreferences)

    def coreference_complex(self, complex_: "FilteredComplex") -> tuple[tuple[Simplex, ...], ...]:
        return self.coreferences_as_simplices(complex_)


@dataclass(frozen=True)
class MorseComplex:
    critical_simplices: tuple[int, ...]
    boundaries: dict[int, Annotation]
    simplex_to_critical: dict[int, int]

    def boundary(self, critical_id: int) -> Annotation:
        return self.boundaries[critical_id]

    def as_simplices(self, complex_: "FilteredComplex") -> tuple[Simplex, ...]:
        return tuple(complex_.vertices(simplex) for simplex in self.critical_simplices)

    def boundary_as_simplices(
        self,
        complex_: "FilteredComplex",
        critical_id: int,
    ) -> tuple[Simplex, ...]:
        return tuple(
            complex_.vertices(self.critical_simplices[boundary_critical_id])
            for boundary_critical_id in self.boundary(critical_id)
        )

    def boundaries_as_simplices(self, complex_: "FilteredComplex") -> dict[Simplex, tuple[Simplex, ...]]:
        return {
            complex_.vertices(simplex): self.boundary_as_simplices(complex_, critical_id)
            for critical_id, simplex in enumerate(self.critical_simplices)
        }


class FilteredComplex:
    """Finite filtered abstract simplicial complex.

    Simplices are inserted explicitly. ``finalize`` checks closure under faces,
    validates filtration monotonicity, assigns simplex ids, and builds boundary
    and coboundary incidence.
    """

    def __init__(self) -> None:
        self._pending: dict[Simplex, float] = {}
        self._records: list[SimplexRecord] = []
        self._simplex_to_id: dict[Simplex, int] = {}
        self._level_values: tuple[float, ...] = ()
        self._level_buckets: tuple[tuple[int, ...], ...] = ()
        self._filtration_order: tuple[int, ...] = ()
        self._cpp: object | None = None
        self._finalized = False

    @classmethod
    def from_simplices(
        cls, simplices: Iterable[tuple[Sequence[int], float]], *, finalize: bool = True
    ) -> "FilteredComplex":
        complex_ = cls()
        for vertices, filtration in simplices:
            complex_.insert(vertices, filtration)
        if finalize:
            complex_.finalize()
        return complex_

    @classmethod
    def from_gudhi_simplex_tree(cls, simplex_tree: object, *, finalize: bool = True) -> "FilteredComplex":
        get_filtration = getattr(simplex_tree, "get_filtration", None)
        if get_filtration is None:
            raise TypeError("Expected a GUDHI-like SimplexTree with get_filtration().")
        return cls.from_simplices(get_filtration(), finalize=finalize)

    @classmethod
    def from_facets(
        cls,
        facets: Iterable[Sequence[int]],
        *,
        simplex_filtration: FiltrationSpec = 0.0,
        finalize: bool = True,
    ) -> "FilteredComplex":
        complex_ = cls()
        for facet in facets:
            for simplex in _all_nonempty_faces(_canonical_simplex(facet)):
                complex_.insert(simplex, _simplex_filtration_value(simplex_filtration, simplex))
        if finalize:
            complex_.finalize()
        return complex_

    @classmethod
    def from_lower_star(
        cls,
        facets: Iterable[Sequence[int]],
        vertex_values: VertexFiltrationSpec,
        *,
        dimension_offset: float = 0.0,
        finalize: bool = True,
    ) -> "FilteredComplex":
        def filtration(simplex: Simplex) -> float:
            return max(_vertex_filtration_value(vertex_values, vertex) for vertex in simplex) + (
                float(dimension_offset) * float(len(simplex) - 1)
            )

        return cls.from_facets(facets, simplex_filtration=filtration, finalize=finalize)

    @classmethod
    def from_graph(
        cls,
        edges: Iterable[Sequence[int]],
        *,
        vertices: Iterable[int] | None = None,
        vertex_filtration: VertexFiltrationSpec = 0.0,
        edge_filtration: EdgeFiltrationSpec | None = None,
        finalize: bool = True,
    ) -> "FilteredComplex":
        complex_ = cls()
        vertex_set = {int(vertex) for vertex in vertices or ()}
        canonical_edges: list[tuple[int, int]] = []
        for edge in edges:
            simplex = _canonical_simplex(edge)
            if len(simplex) != 2:
                raise ValueError("Graph edges must contain exactly two distinct vertices.")
            canonical_edges.append((simplex[0], simplex[1]))
            vertex_set.update(simplex)

        for vertex in sorted(vertex_set):
            complex_.insert([vertex], _vertex_filtration_value(vertex_filtration, vertex))

        for edge in canonical_edges:
            edge_value = (
                max(_vertex_filtration_value(vertex_filtration, vertex) for vertex in edge)
                if edge_filtration is None
                else _edge_filtration_value(edge_filtration, edge)
            )
            complex_.insert(edge, edge_value)

        if finalize:
            complex_.finalize()
        return complex_

    @classmethod
    def from_rips_distance_matrix(
        cls,
        distances: Sequence[Sequence[float]],
        *,
        max_edge_length: float = inf,
        max_dimension: int = 1,
        vertex_filtration: VertexFiltrationSpec = 0.0,
        finalize: bool = True,
    ) -> "FilteredComplex":
        if max_dimension < 0:
            raise ValueError("max_dimension must be non-negative.")
        n = len(distances)
        if any(len(row) != n for row in distances):
            raise ValueError("Distance matrix must be square.")

        complex_ = cls()
        for vertex in range(n):
            complex_.insert([vertex], _vertex_filtration_value(vertex_filtration, vertex))

        for size in range(2, max_dimension + 2):
            for simplex in combinations(range(n), size):
                edge_lengths = [
                    float(distances[i][j])
                    for i, j in combinations(simplex, 2)
                ]
                if not edge_lengths or max(edge_lengths) <= max_edge_length:
                    value = max(
                        max(edge_lengths, default=0.0),
                        max(_vertex_filtration_value(vertex_filtration, vertex) for vertex in simplex),
                    )
                    complex_.insert(simplex, value)

        if finalize:
            complex_.finalize()
        return complex_

    def insert(self, vertices: Sequence[int], filtration: float) -> "FilteredComplex":
        simplex = _canonical_simplex(vertices)
        value = float(filtration)
        old_value = self._pending.get(simplex)
        if old_value is not None and not isclose(old_value, value, abs_tol=1e-12):
            raise ValueError("Duplicate simplex inserted with a different filtration value.")
        self._pending[simplex] = value
        self._cpp = None
        self._finalized = False
        return self

    add_simplex = insert

    def finalize(self) -> "FilteredComplex":
        if not self._pending:
            raise ValueError("Cannot finalize an empty complex.")

        simplex_to_id = {simplex: index for index, simplex in enumerate(sorted(self._pending))}
        level_values = tuple(sorted(set(self._pending.values())))
        level_of = {value: index for index, value in enumerate(level_values)}

        partial: list[dict[str, object]] = []
        for simplex, simplex_id in simplex_to_id.items():
            filtration = self._pending[simplex]
            boundary: list[int] = []
            for face in _codimension_one_faces(simplex):
                try:
                    face_id = simplex_to_id[face]
                except KeyError as exc:
                    raise ValueError("Input is not closed under faces.") from exc
                if self._pending[face] > filtration + 1e-12:
                    raise ValueError("Filtration is not monotone on faces.")
                boundary.append(face_id)

            partial.append(
                {
                    "id": simplex_id,
                    "vertices": simplex,
                    "dimension": len(simplex) - 1,
                    "level": level_of[filtration],
                    "filtration": filtration,
                    "boundary": tuple(boundary),
                    "coboundary": [],
                }
            )

        for record in partial:
            simplex_id = int(record["id"])
            for face_id in record["boundary"]:  # type: ignore[union-attr]
                partial[face_id]["coboundary"].append(simplex_id)  # type: ignore[index,union-attr]

        records = [
            SimplexRecord(
                id=int(record["id"]),
                vertices=record["vertices"],  # type: ignore[arg-type]
                dimension=int(record["dimension"]),
                level=int(record["level"]),
                filtration=float(record["filtration"]),
                boundary=record["boundary"],  # type: ignore[arg-type]
                coboundary=tuple(record["coboundary"]),  # type: ignore[arg-type]
            )
            for record in partial
        ]

        order = tuple(
            sorted(
                range(len(records)),
                key=lambda simplex_id: (
                    records[simplex_id].level,
                    records[simplex_id].dimension,
                    records[simplex_id].vertices,
                ),
            )
        )

        buckets: list[list[int]] = [[] for _ in level_values]
        for simplex_id in order:
            buckets[records[simplex_id].level].append(simplex_id)

        self._records = records
        self._simplex_to_id = simplex_to_id
        self._level_values = level_values
        self._level_buckets = tuple(tuple(bucket) for bucket in buckets)
        self._filtration_order = order
        self._finalized = True
        self._refresh_cpp_backend()
        return self

    def __len__(self) -> int:
        return self.size

    @property
    def size(self) -> int:
        self._ensure_finalized()
        return len(self._records)

    @property
    def num_levels(self) -> int:
        self._ensure_finalized()
        return len(self._level_values)

    @property
    def level_values(self) -> tuple[float, ...]:
        self._ensure_finalized()
        return self._level_values

    @property
    def filtration_order(self) -> tuple[int, ...]:
        self._ensure_finalized()
        return self._filtration_order

    def simplices(self) -> Iterator[SimplexRecord]:
        self._ensure_finalized()
        return iter(self._records)

    def simplex_records(self) -> tuple[SimplexRecord, ...]:
        self._ensure_finalized()
        return tuple(self._records)

    def simplex_list(self) -> tuple[Simplex, ...]:
        self._ensure_finalized()
        return tuple(record.vertices for record in self._records)

    def as_simplices(self) -> tuple[Simplex, ...]:
        return self.simplex_list()

    def filtration_list(self) -> tuple[tuple[Simplex, float], ...]:
        self._ensure_finalized()
        return tuple(
            (self._records[simplex_id].vertices, self._records[simplex_id].filtration)
            for simplex_id in self._filtration_order
        )

    def simplex(self, simplex_id: int) -> SimplexRecord:
        self._ensure_finalized()
        return self._records[simplex_id]

    def vertices(self, simplex_id: int) -> Simplex:
        return self.simplex(simplex_id).vertices

    def dimension(self, simplex_id: int) -> int:
        return self.simplex(simplex_id).dimension

    def level(self, simplex_id: int) -> int:
        return self.simplex(simplex_id).level

    def filtration(self, simplex_id: int) -> float:
        return self.simplex(simplex_id).filtration

    def boundary(self, simplex_id: int) -> tuple[int, ...]:
        return self.simplex(simplex_id).boundary

    def coboundary(self, simplex_id: int) -> tuple[int, ...]:
        return self.simplex(simplex_id).coboundary

    def simplices_of_level(self, level: int) -> tuple[int, ...]:
        self._ensure_finalized()
        return self._level_buckets[level]

    def find_simplex(self, vertices: Sequence[int]) -> int | None:
        self._ensure_finalized()
        return self._simplex_to_id.get(_canonical_simplex(vertices))

    def simplex_id(self, vertices: Sequence[int]) -> int | None:
        return self.find_simplex(vertices)

    def contains(self, vertices: Sequence[int]) -> bool:
        return self.find_simplex(vertices) is not None

    def require_simplex_id(self, vertices: Sequence[int]) -> int:
        simplex = _canonical_simplex(vertices)
        simplex_id = self.find_simplex(simplex)
        if simplex_id is None:
            raise KeyError(f"Simplex {simplex} is not present in the complex.")
        return simplex_id

    def filtration_of(self, vertices: Sequence[int]) -> float:
        return self.filtration(self.require_simplex_id(vertices))

    def boundary_simplices(self, simplex_or_vertices: int | Sequence[int]) -> tuple[Simplex, ...]:
        simplex_id = (
            simplex_or_vertices
            if isinstance(simplex_or_vertices, int)
            else self.require_simplex_id(simplex_or_vertices)
        )
        return tuple(self.vertices(face) for face in self.boundary(simplex_id))

    def coboundary_simplices(self, simplex_or_vertices: int | Sequence[int]) -> tuple[Simplex, ...]:
        simplex_id = (
            simplex_or_vertices
            if isinstance(simplex_or_vertices, int)
            else self.require_simplex_id(simplex_or_vertices)
        )
        return tuple(self.vertices(coface) for coface in self.coboundary(simplex_id))

    def __contains__(self, vertices: object) -> bool:
        if isinstance(vertices, (str, bytes)):
            return False
        try:
            return self.contains(vertices)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return False

    def cpp_backend_active(self) -> bool:
        return self._cpp is not None

    def to_cpp(self) -> object:
        self._ensure_finalized()
        if self._cpp is None:
            raise RuntimeError("The C++ backend is not available for this complex.")
        return self._cpp

    def _ensure_finalized(self) -> None:
        if not self._finalized:
            raise RuntimeError("FilteredComplex.finalize() must be called before querying.")

    def _refresh_cpp_backend(self) -> None:
        if _morse_core is None:
            self._cpp = None
            return
        core = _morse_core.FilteredComplex()
        for simplex, filtration in self._pending.items():
            core.add_simplex(list(simplex), filtration)
        core.finalize()
        self._cpp = core


class CubicalGrid2DComplex:
    """Native-backed filtered two-dimensional cubical grid.

    The grid is built from row-major vertex values. Each edge and square receives
    the maximum filtration value of its vertices, so plateaus are represented
    directly without a lower-star subdivision.
    """

    _simplicial_payload_supported = False

    def __init__(
        self,
        vertex_width: int,
        vertex_height: int,
        vertex_values: Sequence[float],
        *,
        _cpp: object | None = None,
    ) -> None:
        if CppCubicalGrid2DComplex is None:
            raise RuntimeError("The C++ backend is required for cubical complexes.")
        self._cpp = (
            _cpp
            if _cpp is not None
            else CppCubicalGrid2DComplex(
                int(vertex_width),
                int(vertex_height),
                [float(value) for value in vertex_values],
            )
        )
        self._records = self._make_records()
        self._cell_key_to_id = {record.vertices: record.id for record in self._records}

    @classmethod
    def from_vertex_values(
        cls,
        vertex_width: int,
        vertex_height: int,
        vertex_values: Sequence[float],
    ) -> "CubicalGrid2DComplex":
        if CppCubicalGrid2DComplex is None:
            raise RuntimeError("The C++ backend is required for cubical complexes.")
        constructor = getattr(CppCubicalGrid2DComplex, "from_vertex_values", None)
        if constructor is None:
            return cls(vertex_width, vertex_height, vertex_values)
        return cls(
            int(vertex_width),
            int(vertex_height),
            vertex_values,
            _cpp=constructor(
                int(vertex_width),
                int(vertex_height),
                [float(value) for value in vertex_values],
            ),
        )

    def _make_records(self) -> tuple[CubicalCellRecord, ...]:
        records: list[CubicalCellRecord] = []
        for cell in range(int(self._cpp.size)):
            records.append(
                CubicalCellRecord(
                    id=cell,
                    vertices=tuple(int(vertex) for vertex in self._cpp.vertices(cell)),
                    dimension=int(self._cpp.dimension(cell)),
                    level=int(self._cpp.level(cell)),
                    filtration=float(self._cpp.filtration(cell)),
                    boundary=tuple(int(face) for face in self._cpp.boundary(cell)),
                    coboundary=tuple(int(coface) for coface in self._cpp.coboundary(cell)),
                    cell_type=str(self._cpp.cell_type(cell)),
                )
            )
        return tuple(records)

    def __len__(self) -> int:
        return self.size

    @property
    def size(self) -> int:
        return int(self._cpp.size)

    @property
    def num_levels(self) -> int:
        return int(self._cpp.num_levels)

    @property
    def level_values(self) -> tuple[float, ...]:
        return tuple(float(value) for value in self._cpp.level_values)

    @property
    def filtration_order(self) -> tuple[int, ...]:
        return tuple(int(cell) for cell in self._cpp.filtration_order)

    @property
    def vertex_width(self) -> int:
        return int(self._cpp.vertex_width)

    @property
    def vertex_height(self) -> int:
        return int(self._cpp.vertex_height)

    @property
    def square_width(self) -> int:
        return int(self._cpp.square_width)

    @property
    def square_height(self) -> int:
        return int(self._cpp.square_height)

    def simplices(self) -> Iterator[CubicalCellRecord]:
        return iter(self._records)

    def simplex_records(self) -> tuple[CubicalCellRecord, ...]:
        return self._records

    def simplex_list(self) -> tuple[tuple[int, ...], ...]:
        return tuple(record.vertices for record in self._records)

    def as_simplices(self) -> tuple[tuple[int, ...], ...]:
        return self.simplex_list()

    def filtration_list(self) -> tuple[tuple[tuple[int, ...], float], ...]:
        return tuple(
            (self._records[cell].vertices, self._records[cell].filtration)
            for cell in self.filtration_order
        )

    def simplex(self, simplex_id: int) -> CubicalCellRecord:
        return self.cell(simplex_id)

    def cell(self, cell_id: int) -> CubicalCellRecord:
        return self._records[int(cell_id)]

    def vertices(self, cell_id: int) -> tuple[int, ...]:
        return self.cell(cell_id).vertices

    def dimension(self, cell_id: int) -> int:
        return self.cell(cell_id).dimension

    def level(self, cell_id: int) -> int:
        return self.cell(cell_id).level

    def filtration(self, cell_id: int) -> float:
        return self.cell(cell_id).filtration

    def boundary(self, cell_id: int) -> tuple[int, ...]:
        return self.cell(cell_id).boundary

    def coboundary(self, cell_id: int) -> tuple[int, ...]:
        return self.cell(cell_id).coboundary

    def simplices_of_level(self, level: int) -> tuple[int, ...]:
        return tuple(int(cell) for cell in self._cpp.simplices_of_level(int(level)))

    def cell_type(self, cell_id: int) -> str:
        return self.cell(cell_id).cell_type

    def boundary_coefficient(self, cell_id: int, boundary_index: int, modulus: int) -> int:
        return int(self._cpp.boundary_coefficient(int(cell_id), int(boundary_index), int(modulus)))

    def vertex(self, x: int, y: int) -> int:
        return int(self._cpp.vertex(int(x), int(y)))

    def horizontal_edge(self, x: int, y: int) -> int:
        return int(self._cpp.horizontal_edge(int(x), int(y)))

    def vertical_edge(self, x: int, y: int) -> int:
        return int(self._cpp.vertical_edge(int(x), int(y)))

    def square(self, x: int, y: int) -> int:
        return int(self._cpp.square(int(x), int(y)))

    def find_simplex(self, vertices: Sequence[int]) -> int | None:
        return self._cell_key_to_id.get(tuple(sorted(int(vertex) for vertex in vertices)))

    def simplex_id(self, vertices: Sequence[int]) -> int | None:
        return self.find_simplex(vertices)

    def contains(self, vertices: Sequence[int]) -> bool:
        return self.find_simplex(vertices) is not None

    def require_simplex_id(self, vertices: Sequence[int]) -> int:
        cell_key = tuple(sorted(int(vertex) for vertex in vertices))
        cell_id = self.find_simplex(cell_key)
        if cell_id is None:
            raise KeyError(f"Cell with vertices {cell_key} is not present in the cubical grid.")
        return cell_id

    def filtration_of(self, vertices: Sequence[int]) -> float:
        return self.filtration(self.require_simplex_id(vertices))

    def boundary_simplices(self, simplex_or_vertices: int | Sequence[int]) -> tuple[tuple[int, ...], ...]:
        cell_id = (
            simplex_or_vertices
            if isinstance(simplex_or_vertices, int)
            else self.require_simplex_id(simplex_or_vertices)
        )
        return tuple(self.vertices(face) for face in self.boundary(int(cell_id)))

    def coboundary_simplices(self, simplex_or_vertices: int | Sequence[int]) -> tuple[tuple[int, ...], ...]:
        cell_id = (
            simplex_or_vertices
            if isinstance(simplex_or_vertices, int)
            else self.require_simplex_id(simplex_or_vertices)
        )
        return tuple(self.vertices(coface) for coface in self.coboundary(int(cell_id)))

    def __contains__(self, vertices: object) -> bool:
        if isinstance(vertices, (str, bytes)):
            return False
        try:
            return self.contains(vertices)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return False

    def cpp_backend_active(self) -> bool:
        return True

    def to_cpp(self) -> object:
        return self._cpp

    def finalize(self) -> "CubicalGrid2DComplex":
        return self

    def _ensure_finalized(self) -> None:
        return None


class SimplexTreeBuilder:
    """Mutable construction helper that finalizes into a ``FilteredComplex``.

    The builder is intentionally a front-end only: it deduplicates incremental
    insertions and offers simplex-tree-style membership queries, then releases
    to the compact filtered complex used by the Morse and persistence kernels.
    """

    _FILTRATION_KEY = object()

    def __init__(self, *, merge: str = "min", use_cpp: bool = True) -> None:
        if merge not in {"min", "strict"}:
            raise ValueError("merge must be either 'min' or 'strict'.")
        self._merge = merge
        self._cpp = CppSimplexTreeBuilder(merge) if use_cpp and CppSimplexTreeBuilder is not None else None
        self._root: dict[object, object] = {}
        self._size = 0
        self._max_dimension = -1

    @classmethod
    def from_facets(
        cls,
        facets: Iterable[Sequence[int]],
        *,
        simplex_filtration: FiltrationSpec = 0.0,
        merge: str = "min",
    ) -> "SimplexTreeBuilder":
        builder = cls(merge=merge)
        for facet in facets:
            builder.insert_facet(facet, simplex_filtration=simplex_filtration)
        return builder

    @classmethod
    def from_lower_star(
        cls,
        facets: Iterable[Sequence[int]],
        vertex_values: VertexFiltrationSpec,
        *,
        dimension_offset: float = 0.0,
        merge: str = "min",
    ) -> "SimplexTreeBuilder":
        builder = cls(merge=merge)
        for facet in facets:
            builder.insert_lower_star_facet(
                facet,
                vertex_values,
                dimension_offset=dimension_offset,
            )
        return builder

    def insert(
        self,
        vertices: Sequence[int],
        filtration: float = 0.0,
        *,
        include_faces: bool = True,
    ) -> "SimplexTreeBuilder":
        simplex = _canonical_simplex(vertices)
        value = float(filtration)
        if self._cpp is not None:
            self._cpp.insert(list(simplex), value, include_faces=include_faces)
            return self
        if include_faces:
            for face in _all_nonempty_faces(simplex):
                self._insert_canonical(face, value)
        else:
            self._insert_canonical(simplex, value)
        return self

    add_simplex = insert
    insert_simplex = insert

    def insert_simplex_only(
        self,
        vertices: Sequence[int],
        filtration: float = 0.0,
    ) -> "SimplexTreeBuilder":
        simplex = _canonical_simplex(vertices)
        value = float(filtration)
        if self._cpp is not None:
            self._cpp.insert_simplex_only(list(simplex), value)
        else:
            self._insert_canonical(simplex, value)
        return self

    def insert_facet(
        self,
        facet: Sequence[int],
        *,
        simplex_filtration: FiltrationSpec = 0.0,
    ) -> "SimplexTreeBuilder":
        for simplex in _all_nonempty_faces(_canonical_simplex(facet)):
            self._insert_canonical(simplex, _simplex_filtration_value(simplex_filtration, simplex))
        return self

    def insert_lower_star_facet(
        self,
        facet: Sequence[int],
        vertex_values: VertexFiltrationSpec,
        *,
        dimension_offset: float = 0.0,
    ) -> "SimplexTreeBuilder":
        def filtration(simplex: Simplex) -> float:
            return max(_vertex_filtration_value(vertex_values, vertex) for vertex in simplex) + (
                float(dimension_offset) * float(len(simplex) - 1)
            )

        return self.insert_facet(facet, simplex_filtration=filtration)

    def contains(self, vertices: Sequence[int]) -> bool:
        if self._cpp is not None:
            return bool(self._cpp.contains(list(_canonical_simplex(vertices))))
        node = self._find_node(_canonical_simplex(vertices))
        return node is not None and self._FILTRATION_KEY in node

    find_simplex = contains

    def __contains__(self, vertices: object) -> bool:
        if isinstance(vertices, (str, bytes)):
            return False
        try:
            return self.contains(vertices)  # type: ignore[arg-type]
        except (TypeError, ValueError):
            return False

    def filtration(self, vertices: Sequence[int]) -> float:
        simplex = _canonical_simplex(vertices)
        if self._cpp is not None:
            return float(self._cpp.simplex_filtration(list(simplex)))
        node = self._find_node(simplex)
        if node is None or self._FILTRATION_KEY not in node:
            raise KeyError(f"Simplex {simplex} is not present in the builder.")
        return float(node[self._FILTRATION_KEY])

    simplex_filtration = filtration

    def simplices(self) -> tuple[tuple[Simplex, float], ...]:
        if self._cpp is not None:
            return tuple((tuple(simplex), float(value)) for simplex, value in self._cpp.simplices())
        return tuple(self._iter_simplices(self._root, ()))

    def get_filtration(self) -> tuple[tuple[Simplex, float], ...]:
        if self._cpp is not None:
            return tuple((tuple(simplex), float(value)) for simplex, value in self._cpp.get_filtration())
        return tuple(
            sorted(
                self.simplices(),
                key=lambda item: (item[1], len(item[0]), item[0]),
            )
        )

    def to_gudhi_simplex_tree(self) -> object:
        try:
            import gudhi
        except ImportError as exc:
            raise ImportError("GUDHI is required to export a SimplexTree.") from exc

        simplex_tree = gudhi.SimplexTree()
        for simplex, filtration in self.get_filtration():
            simplex_tree.insert(list(simplex), filtration=filtration)
        return simplex_tree

    @property
    def size(self) -> int:
        if self._cpp is not None:
            return int(self._cpp.size)
        return self._size

    def num_simplices(self) -> int:
        return self.size

    def num_vertices(self) -> int:
        if self._cpp is not None:
            return int(self._cpp.num_vertices())
        return len({vertex for simplex, _ in self.simplices() for vertex in simplex})

    @property
    def max_dimension(self) -> int:
        if self._cpp is not None:
            return int(self._cpp.max_dimension)
        return self._max_dimension

    def __len__(self) -> int:
        return self.size

    def clear(self) -> None:
        if self._cpp is not None:
            self._cpp.clear()
        self._root.clear()
        self._size = 0
        self._max_dimension = -1

    def to_filtered_complex(self, *, finalize: bool = True, clear: bool = False) -> FilteredComplex:
        complex_ = FilteredComplex.from_simplices(self.simplices(), finalize=finalize)
        if clear:
            self.clear()
        return complex_

    def finalize(self, *, clear: bool = True) -> FilteredComplex:
        return self.to_filtered_complex(finalize=True, clear=clear)

    def _insert_canonical(self, simplex: Simplex, filtration: float) -> None:
        if self._cpp is not None:
            self._cpp.insert_simplex_only(list(simplex), filtration)
            return
        node = self._ensure_node(simplex)
        if self._FILTRATION_KEY in node:
            old_value = float(node[self._FILTRATION_KEY])
            if self._merge == "strict" and not isclose(old_value, filtration, abs_tol=1e-12):
                raise ValueError("Duplicate simplex inserted with a different filtration value.")
            node[self._FILTRATION_KEY] = min(old_value, filtration)
            return

        node[self._FILTRATION_KEY] = filtration
        self._size += 1
        self._max_dimension = max(self._max_dimension, len(simplex) - 1)

    def _ensure_node(self, simplex: Simplex) -> dict[object, object]:
        node = self._root
        for vertex in simplex:
            child = node.get(vertex)
            if child is None:
                child = {}
                node[vertex] = child
            if not isinstance(child, dict):
                raise RuntimeError("Corrupt simplex-tree builder node.")
            node = child
        return node

    def _find_node(self, simplex: Simplex) -> dict[object, object] | None:
        node = self._root
        for vertex in simplex:
            child = node.get(vertex)
            if not isinstance(child, dict):
                return None
            node = child
        return node

    def _iter_simplices(
        self,
        node: dict[object, object],
        prefix: Simplex,
    ) -> Iterator[tuple[Simplex, float]]:
        if self._FILTRATION_KEY in node:
            yield prefix, float(node[self._FILTRATION_KEY])
        for vertex in sorted(key for key in node if isinstance(key, int)):
            child = node[vertex]
            if not isinstance(child, dict):
                raise RuntimeError("Corrupt simplex-tree builder node.")
            yield from self._iter_simplices(child, prefix + (vertex,))


FilteredComplexBuilder = SimplexTreeBuilder


def compute_morse_sequence(
    complex_: FilteredComplex,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> MorseSequence:
    algorithm = _normalize_morse_sequence_algorithm(algorithm)
    cpp_sequence = _compute_cpp_sequence(complex_, algorithm)
    if cpp_sequence is not None:
        return _sequence_from_cpp_sequence(complex_, cpp_sequence, algorithm=algorithm)
    if algorithm == PLATEAU_GREEDY_SEQUENCE:
        return _compute_plateau_greedy_morse_sequence_python(complex_, algorithm=algorithm)
    if algorithm == COREDUCTION_SEQUENCE:
        return _compute_coreduction_morse_sequence_python(complex_, algorithm=algorithm)
    if algorithm in {F_MAX_SEQUENCE, F_MIN_SEQUENCE}:
        return _compute_paper_f_morse_sequence_python(complex_, algorithm=algorithm)
    if algorithm in {
        FLOODING_MAX_SEQUENCE,
        FLOODING_MIN_SEQUENCE,
        FLOODING_MINMAX_SEQUENCE,
        FLOODING_MAXMIN_SEQUENCE,
    }:
        return _compute_flooding_morse_sequence_python(complex_, algorithm=algorithm)

    n = complex_.size
    steps: list[MorseStep] = []
    critical_simplices: list[int] = []
    critical_index_of_simplex = [-1] * n
    paired_with: list[int | None] = [None] * n
    inserted = [False] * n
    missing_count = [len(complex_.boundary(simplex)) for simplex in range(n)]
    unique_missing = [None] * n
    remaining_by_level = [0] * complex_.num_levels

    for simplex in range(n):
        if missing_count[simplex] == 1:
            unique_missing[simplex] = complex_.boundary(simplex)[0]
        remaining_by_level[complex_.level(simplex)] += 1

    def refresh_unique_missing(simplex: int) -> None:
        unique_missing[simplex] = None
        if missing_count[simplex] != 1:
            return
        for face in complex_.boundary(simplex):
            if not inserted[face]:
                unique_missing[simplex] = face
                return
        raise RuntimeError("Missing-face count is inconsistent.")

    def insert_simplex(simplex: int) -> None:
        if inserted[simplex]:
            raise RuntimeError("Tried to insert a simplex twice.")
        inserted[simplex] = True
        remaining_by_level[complex_.level(simplex)] -= 1
        for coface in complex_.coboundary(simplex):
            if inserted[coface]:
                continue
            if missing_count[coface] == 0:
                raise RuntimeError("Missing-face count underflow.")
            missing_count[coface] -= 1
            refresh_unique_missing(coface)

    for level in range(complex_.num_levels):
        bucket = complex_.simplices_of_level(level)
        while remaining_by_level[level] > 0:
            inserted_pair = False
            for tau in bucket:
                if inserted[tau] or complex_.dimension(tau) == 0 or missing_count[tau] != 1:
                    continue
                sigma = unique_missing[tau]
                if sigma is not None and not inserted[sigma] and complex_.level(sigma) == level:
                    steps.append(MorseStep(REGULAR_PAIR, sigma, tau, level))
                    paired_with[sigma] = tau
                    paired_with[tau] = sigma
                    insert_simplex(sigma)
                    insert_simplex(tau)
                    inserted_pair = True
                    break

            if inserted_pair:
                continue

            fillable = None
            for simplex in bucket:
                if not inserted[simplex] and missing_count[simplex] == 0:
                    fillable = simplex
                    break

            if fillable is None:
                raise RuntimeError("No valid F-sequence step found.")

            critical_id = len(critical_simplices)
            steps.append(MorseStep(CRITICAL, fillable, None, level))
            critical_simplices.append(fillable)
            critical_index_of_simplex[fillable] = critical_id
            insert_simplex(fillable)

    return MorseSequence(
        steps=tuple(steps),
        critical_simplices=tuple(critical_simplices),
        critical_index_of_simplex=tuple(critical_index_of_simplex),
        paired_with=tuple(paired_with),
        algorithm=algorithm,
    )


def _compute_plateau_greedy_morse_sequence_python(
    complex_: FilteredComplex,
    *,
    algorithm: str,
) -> MorseSequence:
    n = complex_.size
    steps: list[MorseStep] = []
    critical_simplices: list[int] = []
    critical_index_of_simplex = [-1] * n
    paired_with: list[int | None] = [None] * n
    inserted = [False] * n
    missing_count = [len(complex_.boundary(simplex)) for simplex in range(n)]
    missing_xor = [0] * n
    remaining_by_level = [0] * complex_.num_levels

    for simplex in range(n):
        xor_value = 0
        for face in complex_.boundary(simplex):
            xor_value ^= face
        missing_xor[simplex] = xor_value
        remaining_by_level[complex_.level(simplex)] += 1

    def is_fillable(simplex: int) -> bool:
        return not inserted[simplex] and missing_count[simplex] == 0

    def is_pairable(tau: int, level: int) -> bool:
        if inserted[tau] or complex_.dimension(tau) == 0 or missing_count[tau] != 1:
            return False
        sigma = missing_xor[tau]
        return sigma < n and not inserted[sigma] and complex_.level(sigma) == level

    for level in range(complex_.num_levels):
        bucket = complex_.simplices_of_level(level)
        pair_candidates: list[tuple[int, int]] = []
        fillable_candidates: list[tuple[int, int, int]] = []

        def count_unlocks(simplex: int) -> int:
            unlocks = 0
            for coface in complex_.coboundary(simplex):
                if (
                    inserted[coface]
                    or complex_.level(coface) != level
                    or missing_count[coface] != 2
                ):
                    continue
                other_missing = missing_xor[coface] ^ simplex
                if (
                    other_missing < n
                    and not inserted[other_missing]
                    and complex_.level(other_missing) == level
                ):
                    unlocks += 1
            return unlocks

        def enqueue_pair_candidate(simplex: int) -> None:
            if complex_.level(simplex) == level and is_pairable(simplex, level):
                heapq.heappush(pair_candidates, (-complex_.dimension(simplex), simplex))

        def enqueue_fillable_candidate(simplex: int) -> None:
            if complex_.level(simplex) == level and is_fillable(simplex):
                heapq.heappush(
                    fillable_candidates,
                    (-count_unlocks(simplex), complex_.dimension(simplex), simplex),
                )

        def enqueue_current_level_candidate(simplex: int) -> None:
            enqueue_pair_candidate(simplex)
            enqueue_fillable_candidate(simplex)

        def enqueue_missing_faces_for_unlocks(simplex: int) -> None:
            if (
                complex_.level(simplex) != level
                or inserted[simplex]
                or missing_count[simplex] != 2
            ):
                return
            for face in complex_.boundary(simplex):
                if not inserted[face]:
                    enqueue_fillable_candidate(face)

        for simplex in bucket:
            enqueue_current_level_candidate(simplex)

        def insert_simplex(simplex: int) -> None:
            if inserted[simplex]:
                raise RuntimeError("Tried to insert a simplex twice.")
            inserted[simplex] = True
            remaining_by_level[complex_.level(simplex)] -= 1
            for coface in complex_.coboundary(simplex):
                if inserted[coface]:
                    continue
                if missing_count[coface] == 0:
                    raise RuntimeError("Missing-face count underflow.")
                missing_count[coface] -= 1
                missing_xor[coface] ^= simplex
                enqueue_current_level_candidate(coface)
                enqueue_missing_faces_for_unlocks(coface)

        while remaining_by_level[level] > 0:
            inserted_pair = False
            while pair_candidates:
                _, tau = heapq.heappop(pair_candidates)
                if not is_pairable(tau, level):
                    continue
                sigma = missing_xor[tau]
                steps.append(MorseStep(REGULAR_PAIR, sigma, tau, level))
                paired_with[sigma] = tau
                paired_with[tau] = sigma
                insert_simplex(sigma)
                insert_simplex(tau)
                inserted_pair = True
                break

            if inserted_pair:
                continue

            fillable = None
            while fillable_candidates:
                neg_unlocks, _, simplex = heapq.heappop(fillable_candidates)
                if not is_fillable(simplex):
                    continue
                current_unlocks = count_unlocks(simplex)
                if current_unlocks != -neg_unlocks:
                    heapq.heappush(
                        fillable_candidates,
                        (-current_unlocks, complex_.dimension(simplex), simplex),
                    )
                    continue
                fillable = simplex
                break

            if fillable is None:
                raise RuntimeError("No valid F-sequence step found.")

            critical_id = len(critical_simplices)
            steps.append(MorseStep(CRITICAL, fillable, None, level))
            critical_simplices.append(fillable)
            critical_index_of_simplex[fillable] = critical_id
            insert_simplex(fillable)

    return MorseSequence(
        steps=tuple(steps),
        critical_simplices=tuple(critical_simplices),
        critical_index_of_simplex=tuple(critical_index_of_simplex),
        paired_with=tuple(paired_with),
        algorithm=algorithm,
    )


def _compute_coreduction_morse_sequence_python(
    complex_: FilteredComplex,
    *,
    algorithm: str,
) -> MorseSequence:
    n = complex_.size
    steps: list[MorseStep] = []
    critical_simplices: list[int] = []
    critical_index_of_simplex = [-1] * n
    paired_with: list[int | None] = [None] * n
    inserted = [False] * n
    active = [False] * n
    active_coface_count = [0] * n

    def emit_critical(simplex: int, level: int) -> None:
        if inserted[simplex]:
            raise RuntimeError("Tried to insert a simplex twice.")
        for face in complex_.boundary(simplex):
            if not inserted[face]:
                raise RuntimeError(
                    "Same-level reduction critical has a missing boundary face."
                )
        critical_id = len(critical_simplices)
        steps.append(MorseStep(CRITICAL, simplex, None, level))
        critical_simplices.append(simplex)
        critical_index_of_simplex[simplex] = critical_id
        inserted[simplex] = True

    def emit_pair(sigma: int, tau: int, level: int) -> None:
        if inserted[sigma] or inserted[tau]:
            raise RuntimeError("Tried to insert a regular pair twice.")
        for face in complex_.boundary(tau):
            if face != sigma and not inserted[face]:
                raise RuntimeError("Same-level reduction pair has a missing boundary face.")
        steps.append(MorseStep(REGULAR_PAIR, sigma, tau, level))
        paired_with[sigma] = tau
        paired_with[tau] = sigma
        inserted[sigma] = True
        inserted[tau] = True

    for level in range(complex_.num_levels):
        bucket = complex_.simplices_of_level(level)
        free_faces: list[int] = []
        collapse_pairs: list[tuple[int, int]] = []

        for simplex in bucket:
            active[simplex] = True
            active_coface_count[simplex] = 0

        def count_active_same_level_cofaces(simplex: int) -> int:
            count = 0
            for coface in complex_.coboundary(simplex):
                if complex_.level(coface) == level and active[coface]:
                    count += 1
            return count

        for simplex in bucket:
            active_coface_count[simplex] = count_active_same_level_cofaces(simplex)
            if active_coface_count[simplex] == 1:
                free_faces.append(simplex)

        def unique_active_same_level_coface(simplex: int) -> int | None:
            result = None
            for coface in complex_.coboundary(simplex):
                if complex_.level(coface) != level or not active[coface]:
                    continue
                if result is not None:
                    return None
                result = coface
            return result

        def decrement_active_faces_of(simplex: int) -> None:
            for face in complex_.boundary(simplex):
                if complex_.level(face) != level or not active[face]:
                    continue
                if active_coface_count[face] == 0:
                    raise RuntimeError("Active coface count underflow.")
                active_coface_count[face] -= 1
                if active_coface_count[face] == 1:
                    free_faces.append(face)

        cursor = 0
        while cursor < len(free_faces):
            sigma = free_faces[cursor]
            cursor += 1
            if not active[sigma] or active_coface_count[sigma] != 1:
                continue
            tau = unique_active_same_level_coface(sigma)
            if tau is None:
                continue
            collapse_pairs.append((sigma, tau))
            active[sigma] = False
            active[tau] = False
            decrement_active_faces_of(tau)
            decrement_active_faces_of(sigma)

        for simplex in bucket:
            if not active[simplex]:
                continue
            active[simplex] = False
            emit_critical(simplex, level)

        for sigma, tau in reversed(collapse_pairs):
            emit_pair(sigma, tau, level)

    return MorseSequence(
        steps=tuple(steps),
        critical_simplices=tuple(critical_simplices),
        critical_index_of_simplex=tuple(critical_index_of_simplex),
        paired_with=tuple(paired_with),
        algorithm=algorithm,
    )


def _compute_paper_f_morse_sequence_python(
    complex_: FilteredComplex,
    *,
    algorithm: str,
) -> MorseSequence:
    n = complex_.size
    steps: list[MorseStep] = []
    critical_simplices: list[int] = []
    critical_index_of_simplex = [-1] * n
    paired_with: list[int | None] = [None] * n

    def emit_critical(simplex: int) -> None:
        critical_id = len(critical_simplices)
        steps.append(MorseStep(CRITICAL, simplex, None, complex_.level(simplex)))
        critical_simplices.append(simplex)
        critical_index_of_simplex[simplex] = critical_id

    def emit_pair(sigma: int, tau: int) -> None:
        if complex_.level(sigma) != complex_.level(tau):
            raise RuntimeError("A paper F-sequence pair crosses filtration levels.")
        steps.append(MorseStep(REGULAR_PAIR, sigma, tau, complex_.level(tau)))
        paired_with[sigma] = tau
        paired_with[tau] = sigma

    order = complex_.filtration_order

    if algorithm == F_MAX_SEQUENCE:
        inserted = [False] * n
        remaining_boundary_count = [len(complex_.boundary(simplex)) for simplex in range(n)]
        candidates: deque[int] = deque(
            simplex for simplex in order if remaining_boundary_count[simplex] == 1
        )

        def unique_remaining_face(tau: int) -> int | None:
            result = None
            for face in complex_.boundary(tau):
                if inserted[face]:
                    continue
                if result is not None:
                    return None
                result = face
            return result

        def decrement_boundary_counts(simplex: int) -> None:
            for coface in complex_.coboundary(simplex):
                if inserted[coface]:
                    continue
                if remaining_boundary_count[coface] == 0:
                    raise RuntimeError("F-Max boundary count underflow.")
                remaining_boundary_count[coface] -= 1
                if remaining_boundary_count[coface] == 1:
                    candidates.append(coface)

        def insert_critical(simplex: int) -> None:
            if inserted[simplex]:
                raise RuntimeError("Tried to insert an F-Max critical twice.")
            for face in complex_.boundary(simplex):
                if not inserted[face]:
                    raise RuntimeError("F-Max critical has a missing boundary face.")
            emit_critical(simplex)
            inserted[simplex] = True
            decrement_boundary_counts(simplex)

        def insert_pair(sigma: int, tau: int) -> None:
            if inserted[sigma] or inserted[tau]:
                raise RuntimeError("Tried to insert an F-Max pair twice.")
            for face in complex_.boundary(tau):
                if face != sigma and not inserted[face]:
                    raise RuntimeError("F-Max pair has a missing boundary face.")
            emit_pair(sigma, tau)
            inserted[sigma] = True
            inserted[tau] = True
            decrement_boundary_counts(sigma)
            decrement_boundary_counts(tau)

        order_index = 0
        while order_index < len(order):
            while candidates:
                tau = candidates.popleft()
                if inserted[tau] or remaining_boundary_count[tau] != 1:
                    continue
                sigma = unique_remaining_face(tau)
                if (
                    sigma is None
                    or inserted[sigma]
                    or complex_.level(sigma) != complex_.level(tau)
                ):
                    continue
                insert_pair(sigma, tau)

            while order_index < len(order) and inserted[order[order_index]]:
                order_index += 1
            if order_index < len(order):
                insert_critical(order[order_index])

    elif algorithm == F_MIN_SEQUENCE:
        removed = [False] * n
        remaining_coboundary_count = [
            len(complex_.coboundary(simplex)) for simplex in range(n)
        ]
        candidates: deque[int] = deque(
            simplex
            for simplex in reversed(order)
            if remaining_coboundary_count[simplex] == 1
        )
        decreasing_events: list[tuple[str, int, int | None]] = []

        def unique_remaining_coface(sigma: int) -> int | None:
            result = None
            for coface in complex_.coboundary(sigma):
                if removed[coface]:
                    continue
                if result is not None:
                    return None
                result = coface
            return result

        def decrement_coboundary_counts(simplex: int) -> None:
            for face in complex_.boundary(simplex):
                if removed[face]:
                    continue
                if remaining_coboundary_count[face] == 0:
                    raise RuntimeError("F-Min coboundary count underflow.")
                remaining_coboundary_count[face] -= 1
                if remaining_coboundary_count[face] == 1:
                    candidates.append(face)

        def remove_critical(simplex: int) -> None:
            if removed[simplex]:
                raise RuntimeError("Tried to remove an F-Min critical twice.")
            decreasing_events.append((CRITICAL, simplex, None))
            removed[simplex] = True
            decrement_coboundary_counts(simplex)

        def remove_pair(sigma: int, tau: int) -> None:
            if removed[sigma] or removed[tau]:
                raise RuntimeError("Tried to remove an F-Min pair twice.")
            if complex_.level(sigma) != complex_.level(tau):
                raise RuntimeError("F-Min pair crosses filtration levels.")
            decreasing_events.append((REGULAR_PAIR, sigma, tau))
            removed[sigma] = True
            removed[tau] = True
            decrement_coboundary_counts(sigma)
            decrement_coboundary_counts(tau)

        reverse_index = len(order)
        while reverse_index > 0:
            while candidates:
                sigma = candidates.popleft()
                if removed[sigma] or remaining_coboundary_count[sigma] != 1:
                    continue
                tau = unique_remaining_coface(sigma)
                if (
                    tau is None
                    or removed[tau]
                    or complex_.level(sigma) != complex_.level(tau)
                ):
                    continue
                remove_pair(sigma, tau)

            while reverse_index > 0 and removed[order[reverse_index - 1]]:
                reverse_index -= 1
            if reverse_index > 0:
                remove_critical(order[reverse_index - 1])

        for kind, sigma, tau in reversed(decreasing_events):
            if kind == CRITICAL:
                emit_critical(sigma)
            else:
                if tau is None:
                    raise RuntimeError("F-Min pair is missing its upper simplex.")
                emit_pair(sigma, tau)
    else:
        raise ValueError(f"Unsupported paper F-sequence algorithm: {algorithm!r}")

    return MorseSequence(
        steps=tuple(steps),
        critical_simplices=tuple(critical_simplices),
        critical_index_of_simplex=tuple(critical_index_of_simplex),
        paired_with=tuple(paired_with),
        algorithm=algorithm,
    )


def _compute_flooding_morse_sequence_python(
    complex_: FilteredComplex,
    *,
    algorithm: str,
) -> MorseSequence:
    n = complex_.size
    steps: list[MorseStep] = []
    critical_simplices: list[int] = []
    critical_index_of_simplex = [-1] * n
    paired_with: list[int | None] = [None] * n
    inserted = [False] * n
    active = [False] * n
    active_boundary_count = [0] * n
    active_coboundary_count = [0] * n
    active_boundary_xor = [0] * n
    active_coboundary_xor = [0] * n

    def emit_critical(simplex: int, level: int) -> None:
        if inserted[simplex]:
            raise RuntimeError("Tried to insert a flooding critical twice.")
        for face in complex_.boundary(simplex):
            if not inserted[face]:
                raise RuntimeError("Flooding critical has a missing boundary face.")
        critical_id = len(critical_simplices)
        steps.append(MorseStep(CRITICAL, simplex, None, level))
        critical_simplices.append(simplex)
        critical_index_of_simplex[simplex] = critical_id
        inserted[simplex] = True

    def emit_pair(sigma: int, tau: int, level: int) -> None:
        if inserted[sigma] or inserted[tau]:
            raise RuntimeError("Tried to insert a flooding regular pair twice.")
        for face in complex_.boundary(tau):
            if face != sigma and not inserted[face]:
                raise RuntimeError("Flooding regular pair has a missing boundary face.")
        steps.append(MorseStep(REGULAR_PAIR, sigma, tau, level))
        paired_with[sigma] = tau
        paired_with[tau] = sigma
        inserted[sigma] = True
        inserted[tau] = True

    for level in range(complex_.num_levels):
        bucket = complex_.simplices_of_level(level)
        level_order = {simplex: index for index, simplex in enumerate(bucket)}
        coreduction_candidates: list[tuple[int, int, int]] = []
        reduction_candidates: list[tuple[int, int, int]] = []
        coperforation_candidates: list[tuple[int, int, int]] = []
        perforation_candidates: list[tuple[int, int, int]] = []
        left_events: list[tuple[str, int, int | None]] = []
        right_events: list[tuple[str, int, int | None]] = []
        remaining = len(bucket)

        for simplex in bucket:
            active[simplex] = True
            active_boundary_count[simplex] = 0
            active_coboundary_count[simplex] = 0
            active_boundary_xor[simplex] = 0
            active_coboundary_xor[simplex] = 0

        for simplex in bucket:
            for face in complex_.boundary(simplex):
                if complex_.level(face) == level and active[face]:
                    active_boundary_count[simplex] += 1
                    active_boundary_xor[simplex] ^= face
            for coface in complex_.coboundary(simplex):
                if complex_.level(coface) == level and active[coface]:
                    active_coboundary_count[simplex] += 1
                    active_coboundary_xor[simplex] ^= coface

        def is_coreduction(tau: int) -> bool:
            sigma = active_boundary_xor[tau]
            return active[tau] and active_boundary_count[tau] == 1 and sigma < n and active[sigma]

        def is_reduction(sigma: int) -> bool:
            tau = active_coboundary_xor[sigma]
            return active[sigma] and active_coboundary_count[sigma] == 1 and tau < n and active[tau]

        def is_coperforation(simplex: int) -> bool:
            return active[simplex] and active_boundary_count[simplex] == 0

        def is_perforation(simplex: int) -> bool:
            return active[simplex] and active_coboundary_count[simplex] == 0

        def maximal_priority(simplex: int) -> tuple[int, int, int]:
            return (-complex_.dimension(simplex), -level_order[simplex], simplex)

        def minimal_priority(simplex: int) -> tuple[int, int, int]:
            return (complex_.dimension(simplex), level_order[simplex], simplex)

        def enqueue_if_candidate(simplex: int) -> None:
            if not active[simplex]:
                return
            if is_coreduction(simplex):
                heapq.heappush(coreduction_candidates, maximal_priority(simplex))
            if is_reduction(simplex):
                heapq.heappush(reduction_candidates, minimal_priority(simplex))
            if is_coperforation(simplex):
                heapq.heappush(coperforation_candidates, maximal_priority(simplex))
            if is_perforation(simplex):
                heapq.heappush(perforation_candidates, minimal_priority(simplex))

        for simplex in bucket:
            enqueue_if_candidate(simplex)

        def remove_simplex(simplex: int) -> None:
            nonlocal remaining
            if not active[simplex]:
                raise RuntimeError("Tried to remove an inactive flooding simplex.")
            active[simplex] = False
            remaining -= 1

            for coface in complex_.coboundary(simplex):
                if complex_.level(coface) != level or not active[coface]:
                    continue
                if active_boundary_count[coface] == 0:
                    raise RuntimeError("Active boundary count underflow.")
                active_boundary_count[coface] -= 1
                active_boundary_xor[coface] ^= simplex
                enqueue_if_candidate(coface)

            for face in complex_.boundary(simplex):
                if complex_.level(face) != level or not active[face]:
                    continue
                if active_coboundary_count[face] == 0:
                    raise RuntimeError("Active coboundary count underflow.")
                active_coboundary_count[face] -= 1
                active_coboundary_xor[face] ^= simplex
                enqueue_if_candidate(face)

        def take_coreduction() -> bool:
            while coreduction_candidates:
                tau = heapq.heappop(coreduction_candidates)[2]
                if not is_coreduction(tau):
                    continue
                sigma = active_boundary_xor[tau]
                left_events.append((REGULAR_PAIR, sigma, tau))
                remove_simplex(sigma)
                remove_simplex(tau)
                return True
            return False

        def take_reduction() -> bool:
            while reduction_candidates:
                sigma = heapq.heappop(reduction_candidates)[2]
                if not is_reduction(sigma):
                    continue
                tau = active_coboundary_xor[sigma]
                right_events.append((REGULAR_PAIR, sigma, tau))
                remove_simplex(sigma)
                remove_simplex(tau)
                return True
            return False

        def take_coperforation() -> bool:
            while coperforation_candidates:
                simplex = heapq.heappop(coperforation_candidates)[2]
                if not is_coperforation(simplex):
                    continue
                left_events.append((CRITICAL, simplex, None))
                remove_simplex(simplex)
                return True
            return False

        def take_perforation() -> bool:
            while perforation_candidates:
                simplex = heapq.heappop(perforation_candidates)[2]
                if not is_perforation(simplex):
                    continue
                right_events.append((CRITICAL, simplex, None))
                remove_simplex(simplex)
                return True
            return False

        while remaining > 0:
            if algorithm == FLOODING_MAX_SEQUENCE:
                removed = take_coreduction() or take_coperforation()
            elif algorithm == FLOODING_MIN_SEQUENCE:
                removed = take_reduction() or take_perforation()
            elif algorithm == FLOODING_MINMAX_SEQUENCE:
                removed = (
                    take_reduction()
                    or take_coreduction()
                    or take_perforation()
                    or take_coperforation()
                )
            else:
                removed = (
                    take_coreduction()
                    or take_reduction()
                    or take_coperforation()
                    or take_perforation()
                )
            if not removed:
                raise RuntimeError("No valid flooding operation found.")

        for kind, sigma, tau in left_events:
            if kind == CRITICAL:
                emit_critical(sigma, level)
            else:
                if tau is None:
                    raise RuntimeError("Regular pair is missing its upper simplex.")
                emit_pair(sigma, tau, level)
        for kind, sigma, tau in reversed(right_events):
            if kind == CRITICAL:
                emit_critical(sigma, level)
            else:
                if tau is None:
                    raise RuntimeError("Regular pair is missing its upper simplex.")
                emit_pair(sigma, tau, level)

    return MorseSequence(
        steps=tuple(steps),
        critical_simplices=tuple(critical_simplices),
        critical_index_of_simplex=tuple(critical_index_of_simplex),
        paired_with=tuple(paired_with),
        algorithm=algorithm,
    )


def compute_morse_sequence_and_reference_map(
    complex_: FilteredComplex,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> MorseReferenceFrame:
    algorithm = _normalize_morse_sequence_algorithm(algorithm)
    cpp_frame = _compute_cpp_reference_frame(complex_, algorithm)
    if cpp_frame is not None:
        cpp_sequence = cpp_frame.sequence
        return MorseReferenceFrame(
            sequence=_sequence_from_cpp_sequence(complex_, cpp_sequence, algorithm=algorithm),
            _cpp_frame=cpp_frame,
        )
    if algorithm == PLATEAU_GREEDY_SEQUENCE:
        sequence = _compute_plateau_greedy_morse_sequence_python(complex_, algorithm=algorithm)
        return MorseReferenceFrame(
            sequence=sequence,
            _references=compute_reference_map(complex_, sequence, algorithm=algorithm),
        )
    if algorithm == COREDUCTION_SEQUENCE:
        sequence = _compute_coreduction_morse_sequence_python(complex_, algorithm=algorithm)
        return MorseReferenceFrame(
            sequence=sequence,
            _references=compute_reference_map(complex_, sequence, algorithm=algorithm),
        )
    if algorithm in {F_MAX_SEQUENCE, F_MIN_SEQUENCE}:
        sequence = _compute_paper_f_morse_sequence_python(complex_, algorithm=algorithm)
        return MorseReferenceFrame(
            sequence=sequence,
            _references=compute_reference_map(complex_, sequence, algorithm=algorithm),
        )
    if algorithm in {
        FLOODING_MAX_SEQUENCE,
        FLOODING_MIN_SEQUENCE,
        FLOODING_MINMAX_SEQUENCE,
        FLOODING_MAXMIN_SEQUENCE,
    }:
        sequence = _compute_flooding_morse_sequence_python(complex_, algorithm=algorithm)
        return MorseReferenceFrame(
            sequence=sequence,
            _references=compute_reference_map(complex_, sequence, algorithm=algorithm),
        )
    return _compute_saturated_morse_reference_frame_python(complex_, algorithm=algorithm)


def compute_morse_sequence_and_coreference_map(
    complex_: FilteredComplex,
    *,
    algorithm: str = COREDUCTION_SEQUENCE,
) -> MorseCoreferenceFrame:
    algorithm = _normalize_morse_sequence_algorithm(algorithm)
    cpp_frame = _compute_cpp_coreference_frame(complex_, algorithm)
    if cpp_frame is not None:
        cpp_sequence = cpp_frame.sequence
        return MorseCoreferenceFrame(
            sequence=_sequence_from_cpp_sequence(complex_, cpp_sequence, algorithm=algorithm),
            _cpp_frame=cpp_frame,
        )

    sequence = compute_morse_sequence(complex_, algorithm=algorithm)
    return MorseCoreferenceFrame(
        sequence=sequence,
        _coreferences=compute_coreference_map(complex_, sequence, algorithm=algorithm),
    )


def compute_reference_map(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> tuple[Annotation, ...]:
    if sequence is None:
        return compute_morse_sequence_and_reference_map(complex_, algorithm=algorithm).references

    cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
    if cpp_sequence is not None and complex_._cpp is not None:
        return _references_from_cpp_annotations(
            _morse_core.compute_reference_map(complex_._cpp, cpp_sequence)
        )

    references: list[Annotation] = [() for _ in range(complex_.size)]

    for step in sequence.steps:
        if step.type == CRITICAL:
            critical_id = sequence.critical_index(step.sigma)
            if critical_id < 0:
                raise RuntimeError("Expected a critical simplex.")
            references[step.sigma] = (critical_id,)
            continue

        if step.tau is None:
            raise RuntimeError("Regular pair is missing its upper simplex.")
        references[step.tau] = ()
        lower_reference: Annotation = ()
        for face in complex_.boundary(step.tau):
            if face != step.sigma:
                lower_reference = _xor_sorted(lower_reference, references[face])
        references[step.sigma] = lower_reference

    return tuple(references)


def compute_reference_map_modp(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    *,
    modulus: int,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> tuple[FieldAnnotation, ...]:
    """Compute the Morse reference map over the prime field ``F_modulus``."""

    modulus = _validate_prime_modulus(modulus)
    if modulus == 2:
        z2_references = compute_reference_map(complex_, sequence, algorithm=algorithm)
        return tuple(_field_annotation_from_labels(annotation) for annotation in z2_references)

    cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
    runner = (
        getattr(_morse_core, "compute_reference_map_modp", None)
        if cpp_sequence is not None and complex_._cpp is not None
        else None
    )
    if runner is not None:
        return _field_references_from_cpp_annotations(
            runner(complex_._cpp, cpp_sequence, modulus)
        )

    sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
    references: list[dict[int, int]] = [{} for _ in range(complex_.size)]

    for step in sequence.steps:
        if step.type == CRITICAL:
            critical_id = sequence.critical_index(step.sigma)
            if critical_id < 0:
                raise RuntimeError("Expected a critical simplex.")
            references[step.sigma] = {critical_id: 1}
            continue

        if step.tau is None:
            raise RuntimeError("Regular pair is missing its upper simplex.")
        references[step.tau] = {}

        paired_face_coefficient: int | None = None
        lower_reference: dict[int, int] = {}
        for removed_index, face in enumerate(complex_.boundary(step.tau)):
            coefficient = _boundary_coefficient_for(
                complex_,
                step.tau,
                removed_index,
                modulus,
            )
            if face == step.sigma:
                paired_face_coefficient = coefficient
                continue
            _modp_add_scaled(lower_reference, references[face], coefficient, modulus)

        if paired_face_coefficient is None:
            raise RuntimeError("Regular pair is not a face/coface pair.")

        scale = (-_modp_inverse(paired_face_coefficient, modulus)) % modulus
        _modp_scale_in_place(lower_reference, scale, modulus)
        references[step.sigma] = lower_reference

    return tuple(_field_annotation_from_dict(annotation) for annotation in references)


def compute_coreference_map(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> tuple[Annotation, ...]:
    if sequence is None:
        return compute_morse_sequence_and_coreference_map(
            complex_,
            algorithm=algorithm,
        ).coreferences

    cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
    if cpp_sequence is not None and complex_._cpp is not None:
        core_runner = getattr(_morse_core, "compute_coreference_map", None)
        if core_runner is not None:
            return _references_from_cpp_annotations(core_runner(complex_._cpp, cpp_sequence))

    coreferences: list[Annotation] = [() for _ in range(complex_.size)]

    for step in reversed(sequence.steps):
        if step.type == CRITICAL:
            critical_id = sequence.critical_index(step.sigma)
            if critical_id < 0:
                raise RuntimeError("Expected a critical simplex.")
            coreferences[step.sigma] = (critical_id,)
            continue

        if step.tau is None:
            raise RuntimeError("Regular pair is missing its upper simplex.")
        coreferences[step.sigma] = ()
        upper_coreference: Annotation = ()
        for coface in complex_.coboundary(step.sigma):
            if coface != step.tau:
                upper_coreference = _xor_sorted(upper_coreference, coreferences[coface])
        coreferences[step.tau] = upper_coreference

    return tuple(coreferences)


def compute_coreference_map_modp(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    *,
    modulus: int,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> tuple[FieldAnnotation, ...]:
    """Compute the Morse coreference map over the prime field ``F_modulus``."""

    modulus = _validate_prime_modulus(modulus)
    if modulus == 2:
        z2_coreferences = compute_coreference_map(complex_, sequence, algorithm=algorithm)
        return tuple(_field_annotation_from_labels(annotation) for annotation in z2_coreferences)

    cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
    runner = (
        getattr(_morse_core, "compute_coreference_map_modp", None)
        if cpp_sequence is not None and complex_._cpp is not None
        else None
    )
    if runner is not None:
        return _field_references_from_cpp_annotations(
            runner(complex_._cpp, cpp_sequence, modulus)
        )

    sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
    coreferences: list[dict[int, int]] = [{} for _ in range(complex_.size)]

    for step in reversed(sequence.steps):
        if step.type == CRITICAL:
            critical_id = sequence.critical_index(step.sigma)
            if critical_id < 0:
                raise RuntimeError("Expected a critical simplex.")
            coreferences[step.sigma] = {critical_id: 1}
            continue

        if step.tau is None:
            raise RuntimeError("Regular pair is missing its upper simplex.")
        coreferences[step.sigma] = {}

        paired_coface_coefficient: int | None = None
        upper_coreference: dict[int, int] = {}
        for coface in complex_.coboundary(step.sigma):
            coefficient = _boundary_coefficient_of_face(
                complex_,
                coface=coface,
                face=step.sigma,
                modulus=modulus,
            )
            if coface == step.tau:
                paired_coface_coefficient = coefficient
                continue
            _modp_add_scaled(upper_coreference, coreferences[coface], coefficient, modulus)

        if paired_coface_coefficient is None:
            raise RuntimeError("Regular pair is not a face/coface pair.")

        scale = (-_modp_inverse(paired_coface_coefficient, modulus)) % modulus
        _modp_scale_in_place(upper_coreference, scale, modulus)
        coreferences[step.tau] = upper_coreference

    return tuple(_field_annotation_from_dict(annotation) for annotation in coreferences)


def morse_sequence_as_simplices(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> tuple[MorseStepSimplices, ...]:
    sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
    return sequence.steps_as_simplices(complex_)


def annotation_as_simplices(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    annotation: Sequence[int],
) -> tuple[Simplex, ...]:
    return tuple(
        complex_.vertices(sequence.critical_simplices[critical_id])
        for critical_id in annotation
    )


def reference_map_as_simplices(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    references: Sequence[Annotation] | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> tuple[tuple[Simplex, ...], ...]:
    if sequence is None or references is None:
        frame = compute_morse_sequence_and_reference_map(complex_, algorithm=algorithm)
        if sequence is None:
            sequence = frame.sequence
        if references is None:
            references = frame.references
    return tuple(annotation_as_simplices(complex_, sequence, annotation) for annotation in references)


def coreference_map_as_simplices(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    coreferences: Sequence[Annotation] | None = None,
    *,
    algorithm: str = COREDUCTION_SEQUENCE,
) -> tuple[tuple[Simplex, ...], ...]:
    if sequence is None or coreferences is None:
        frame = compute_morse_sequence_and_coreference_map(complex_, algorithm=algorithm)
        if sequence is None:
            sequence = frame.sequence
        if coreferences is None:
            coreferences = frame.coreferences
    return tuple(
        annotation_as_simplices(complex_, sequence, annotation)
        for annotation in coreferences
    )


def compute_reference_complex(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    references: Sequence[Annotation] | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> tuple[tuple[Simplex, ...], ...]:
    return reference_map_as_simplices(
        complex_,
        sequence,
        references,
        algorithm=algorithm,
    )


def compute_coreference_complex(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    coreferences: Sequence[Annotation] | None = None,
    *,
    algorithm: str = COREDUCTION_SEQUENCE,
) -> tuple[tuple[Simplex, ...], ...]:
    return coreference_map_as_simplices(
        complex_,
        sequence,
        coreferences,
        algorithm=algorithm,
    )


def compute_morse_complex(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    references: Sequence[Annotation] | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> MorseComplex:
    if sequence is None and references is None:
        frame = compute_morse_sequence_and_reference_map(complex_, algorithm=algorithm)
        sequence = frame.sequence
        references = frame.references
    else:
        sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
        references = references or compute_reference_map(complex_, sequence, algorithm=algorithm)
    boundaries: dict[int, Annotation] = {}
    for critical_id, sigma in enumerate(sequence.critical_simplices):
        boundary_annotation: Annotation = ()
        for face in complex_.boundary(sigma):
            boundary_annotation = _xor_sorted(boundary_annotation, references[face])
        boundaries[critical_id] = boundary_annotation
    return MorseComplex(
        critical_simplices=sequence.critical_simplices,
        boundaries=boundaries,
        simplex_to_critical={
            simplex: critical_id for critical_id, simplex in enumerate(sequence.critical_simplices)
        },
    )


def compute_morse_persistence(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    references: Sequence[Annotation] | Sequence[FieldAnnotation] | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
    modulus: int = 2,
) -> PersistenceDiagram:
    modulus = _validate_prime_modulus(modulus)
    if modulus != 2:
        cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
        runner = (
            getattr(_morse_core, "compute_morse_persistence_modp", None)
            if references is None and cpp_sequence is not None and complex_._cpp is not None
            else None
        )
        if runner is not None:
            return _diagram_from_cpp_dict(runner(complex_._cpp, cpp_sequence, modulus))
        return _compute_morse_persistence_modp(
            complex_,
            sequence=sequence,
            references=references,  # type: ignore[arg-type]
            algorithm=algorithm,
            modulus=modulus,
        )

    if sequence is None and references is None:
        frame = compute_morse_sequence_and_reference_map(complex_, algorithm=algorithm)
        if (
            _morse_core is not None
            and complex_._cpp is not None
            and frame.sequence._cpp_sequence is not None
            and (frame._cpp_frame is not None or frame._cpp_references is not None)
        ):
            if frame._cpp_frame is not None:
                return _diagram_from_cpp_dict(
                    _morse_core.reduce_morse_reference_frame_object(
                        complex_._cpp,
                        frame._cpp_frame,
                    )
                )
            return _diagram_from_cpp_dict(
                _morse_core.reduce_morse_persistence_object(
                    complex_._cpp,
                    frame.sequence._cpp_sequence,
                    frame._cpp_references,
                )
            )
        sequence = frame.sequence
        references = frame.references

    cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
    if cpp_sequence is not None and complex_._cpp is not None:
        if references is None:
            return _diagram_from_cpp_dict(
                _morse_core.compute_morse_persistence(complex_._cpp, cpp_sequence)
            )
        return _diagram_from_cpp_dict(
            _morse_core.reduce_morse_persistence(
                complex_._cpp,
                cpp_sequence,
                [list(annotation) for annotation in references],
            )
        )

    sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
    references = list(references or compute_reference_map(complex_, sequence, algorithm=algorithm))
    critical_simplices = sequence.critical_simplices
    active = [False] * len(critical_simplices)
    killed = [False] * len(critical_simplices)
    finite_pairs: list[PersistencePair] = []

    for sigma in critical_simplices:
        sigma_critical_id = sequence.critical_index(sigma)
        if sigma_critical_id < 0:
            raise RuntimeError("Expected a critical simplex.")

        boundary_annotation: Annotation = ()
        for face in complex_.boundary(sigma):
            boundary_annotation = _xor_sorted(boundary_annotation, references[face])

        if not boundary_annotation:
            active[sigma_critical_id] = True
            continue

        pivot = boundary_annotation[-1]
        birth = critical_simplices[pivot]
        finite_pairs.append(
            PersistencePair(
                birth=birth,
                death=sigma,
                dimension=complex_.dimension(birth),
                birth_value=complex_.filtration(birth),
                death_value=complex_.filtration(sigma),
            )
        )
        killed[pivot] = True
        _remove_label_from_all(references, sigma_critical_id)
        _xor_into_all_containing(references, pivot, boundary_annotation)

    essential = [
        EssentialInterval(
            birth=birth,
            dimension=complex_.dimension(birth),
            birth_value=complex_.filtration(birth),
        )
        for critical_id, birth in enumerate(critical_simplices)
        if active[critical_id] and not killed[critical_id]
    ]
    return PersistenceDiagram(tuple(finite_pairs), tuple(essential))


def compute_morse_persistence_modp(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    references: Sequence[FieldAnnotation] | None = None,
    *,
    modulus: int,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> PersistenceDiagram:
    """Compute Morse persistence over the prime field ``F_modulus``."""

    return compute_morse_persistence(
        complex_,
        sequence=sequence,
        references=references,
        algorithm=algorithm,
        modulus=modulus,
    )


def _compute_morse_persistence_modp(
    complex_: FilteredComplex,
    *,
    sequence: MorseSequence | None,
    references: Sequence[FieldAnnotation] | Sequence[Annotation] | None,
    algorithm: str,
    modulus: int,
) -> PersistenceDiagram:
    sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
    if references is None:
        reference_columns = [
            _field_annotation_to_dict(annotation, modulus)
            for annotation in compute_reference_map_modp(
                complex_,
                sequence,
                algorithm=algorithm,
                modulus=modulus,
            )
        ]
    else:
        reference_columns = [
            _field_annotation_to_dict(annotation, modulus)
            for annotation in references
        ]

    critical_simplices = sequence.critical_simplices
    active = [False] * len(critical_simplices)
    killed = [False] * len(critical_simplices)
    finite_pairs: list[PersistencePair] = []

    for sigma in critical_simplices:
        sigma_critical_id = sequence.critical_index(sigma)
        if sigma_critical_id < 0:
            raise RuntimeError("Expected a critical simplex.")

        boundary_annotation: dict[int, int] = {}
        _add_oriented_boundary_references(
            boundary_annotation,
            complex_,
            sigma,
            reference_columns,
            modulus,
        )

        pivot = _modp_low(boundary_annotation)
        if pivot is None:
            active[sigma_critical_id] = True
            continue

        birth = critical_simplices[pivot]
        finite_pairs.append(
            PersistencePair(
                birth=birth,
                death=sigma,
                dimension=complex_.dimension(birth),
                birth_value=complex_.filtration(birth),
                death_value=complex_.filtration(sigma),
            )
        )
        killed[pivot] = True

        for annotation in reference_columns:
            annotation.pop(sigma_critical_id, None)

        pivot_coefficient = boundary_annotation[pivot]
        inverse_pivot = _modp_inverse(pivot_coefficient, modulus)
        for annotation in reference_columns:
            coefficient = annotation.get(pivot)
            if coefficient is None:
                continue
            scale = (-coefficient * inverse_pivot) % modulus
            _modp_add_scaled(annotation, boundary_annotation, scale, modulus)

    essential = [
        EssentialInterval(
            birth=birth,
            dimension=complex_.dimension(birth),
            birth_value=complex_.filtration(birth),
        )
        for critical_id, birth in enumerate(critical_simplices)
        if active[critical_id] and not killed[critical_id]
    ]
    return PersistenceDiagram(tuple(finite_pairs), tuple(essential))


def compute_morse_coreference_persistence(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    coreferences: Sequence[Annotation] | Sequence[FieldAnnotation] | None = None,
    *,
    algorithm: str = COREDUCTION_SEQUENCE,
    modulus: int = 2,
) -> PersistenceDiagram:
    modulus = _validate_prime_modulus(modulus)
    if modulus != 2:
        cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
        runner = (
            getattr(_morse_core, "compute_morse_coreference_persistence_modp", None)
            if coreferences is None and cpp_sequence is not None and complex_._cpp is not None
            else None
        )
        if runner is not None:
            return _diagram_from_cpp_dict(runner(complex_._cpp, cpp_sequence, modulus))
        return _compute_morse_coreference_persistence_modp(
            complex_,
            sequence=sequence,
            coreferences=coreferences,  # type: ignore[arg-type]
            algorithm=algorithm,
            modulus=modulus,
        )

    if sequence is None and coreferences is None:
        frame = compute_morse_sequence_and_coreference_map(complex_, algorithm=algorithm)
        if (
            _morse_core is not None
            and complex_._cpp is not None
            and frame.sequence._cpp_sequence is not None
            and frame._cpp_frame is not None
        ):
            reducer = getattr(_morse_core, "reduce_morse_coreference_frame_object", None)
            if reducer is not None:
                return _diagram_from_cpp_dict(reducer(complex_._cpp, frame._cpp_frame))
        sequence = frame.sequence
        coreferences = frame.coreferences

    cpp_sequence = _cpp_sequence_for(complex_, sequence, algorithm=algorithm)
    if cpp_sequence is not None and complex_._cpp is not None and coreferences is None:
        runner = getattr(_morse_core, "compute_morse_coreference_persistence", None)
        if runner is not None:
            return _diagram_from_cpp_dict(runner(complex_._cpp, cpp_sequence))

    sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
    coreferences = list(
        coreferences or compute_coreference_map(complex_, sequence, algorithm=algorithm)
    )
    critical_simplices = sequence.critical_simplices
    active_dual = [False] * len(critical_simplices)
    killed_dual = [False] * len(critical_simplices)
    finite_pairs: list[PersistencePair] = []

    for sigma in reversed(critical_simplices):
        sigma_critical_id = sequence.critical_index(sigma)
        if sigma_critical_id < 0:
            raise RuntimeError("Expected a critical simplex.")

        coboundary_annotation: Annotation = ()
        for coface in complex_.coboundary(sigma):
            coboundary_annotation = _xor_sorted(
                coboundary_annotation,
                coreferences[coface],
            )

        if not coboundary_annotation:
            active_dual[sigma_critical_id] = True
            continue

        pivot = coboundary_annotation[0]
        death = critical_simplices[pivot]
        finite_pairs.append(
            PersistencePair(
                birth=sigma,
                death=death,
                dimension=complex_.dimension(sigma),
                birth_value=complex_.filtration(sigma),
                death_value=complex_.filtration(death),
            )
        )
        killed_dual[pivot] = True
        _remove_label_from_all(coreferences, sigma_critical_id)
        _xor_into_all_containing(coreferences, pivot, coboundary_annotation)

    essential = [
        EssentialInterval(
            birth=birth,
            dimension=complex_.dimension(birth),
            birth_value=complex_.filtration(birth),
        )
        for critical_id, birth in enumerate(critical_simplices)
        if active_dual[critical_id] and not killed_dual[critical_id]
    ]
    return PersistenceDiagram(tuple(finite_pairs), tuple(essential))


def compute_morse_coreference_persistence_modp(
    complex_: FilteredComplex,
    sequence: MorseSequence | None = None,
    coreferences: Sequence[FieldAnnotation] | None = None,
    *,
    modulus: int,
    algorithm: str = COREDUCTION_SEQUENCE,
) -> PersistenceDiagram:
    """Compute Morse coreference persistence over the prime field ``F_modulus``."""

    return compute_morse_coreference_persistence(
        complex_,
        sequence=sequence,
        coreferences=coreferences,
        algorithm=algorithm,
        modulus=modulus,
    )


def _compute_morse_coreference_persistence_modp(
    complex_: FilteredComplex,
    *,
    sequence: MorseSequence | None,
    coreferences: Sequence[FieldAnnotation] | Sequence[Annotation] | None,
    algorithm: str,
    modulus: int,
) -> PersistenceDiagram:
    sequence = sequence or compute_morse_sequence(complex_, algorithm=algorithm)
    if coreferences is None:
        coreference_columns = [
            _field_annotation_to_dict(annotation, modulus)
            for annotation in compute_coreference_map_modp(
                complex_,
                sequence,
                algorithm=algorithm,
                modulus=modulus,
            )
        ]
    else:
        coreference_columns = [
            _field_annotation_to_dict(annotation, modulus)
            for annotation in coreferences
        ]

    critical_simplices = sequence.critical_simplices
    active_dual = [False] * len(critical_simplices)
    killed_dual = [False] * len(critical_simplices)
    finite_pairs: list[PersistencePair] = []

    for sigma in reversed(critical_simplices):
        sigma_critical_id = sequence.critical_index(sigma)
        if sigma_critical_id < 0:
            raise RuntimeError("Expected a critical simplex.")

        coboundary_annotation: dict[int, int] = {}
        _add_oriented_coboundary_coreferences(
            coboundary_annotation,
            complex_,
            sigma,
            coreference_columns,
            modulus,
        )

        pivot = _modp_min(coboundary_annotation)
        if pivot is None:
            active_dual[sigma_critical_id] = True
            continue

        death = critical_simplices[pivot]
        finite_pairs.append(
            PersistencePair(
                birth=sigma,
                death=death,
                dimension=complex_.dimension(sigma),
                birth_value=complex_.filtration(sigma),
                death_value=complex_.filtration(death),
            )
        )
        killed_dual[pivot] = True

        for annotation in coreference_columns:
            annotation.pop(sigma_critical_id, None)

        pivot_coefficient = coboundary_annotation[pivot]
        inverse_pivot = _modp_inverse(pivot_coefficient, modulus)
        for annotation in coreference_columns:
            coefficient = annotation.get(pivot)
            if coefficient is None:
                continue
            scale = (-coefficient * inverse_pivot) % modulus
            _modp_add_scaled(annotation, coboundary_annotation, scale, modulus)

    essential = [
        EssentialInterval(
            birth=birth,
            dimension=complex_.dimension(birth),
            birth_value=complex_.filtration(birth),
        )
        for critical_id, birth in enumerate(critical_simplices)
        if active_dual[critical_id] and not killed_dual[critical_id]
    ]
    return PersistenceDiagram(tuple(finite_pairs), tuple(essential))


def compute_standard_persistence(
    complex_: FilteredComplex,
    *,
    modulus: int = 2,
) -> PersistenceDiagram:
    """Compute ordinary persistence over the prime field ``F_modulus``.

    The existing C++ path is used for the default ``Z2`` case.  Odd prime
    fields currently use the pure-Python oriented-boundary reducer.
    Composite moduli are rejected because the barcode API assumes a field.
    """

    modulus = _validate_prime_modulus(modulus)
    if modulus != 2:
        runner = (
            getattr(_morse_core, "compute_standard_persistence_modp", None)
            if _morse_core is not None and complex_._cpp is not None
            else None
        )
        if runner is not None:
            return _diagram_from_cpp_dict(runner(complex_._cpp, modulus))
        return _compute_standard_persistence_modp(complex_, modulus)

    if _morse_core is not None and complex_._cpp is not None:
        return _diagram_from_cpp_dict(_morse_core.compute_standard_persistence(complex_._cpp))

    order = complex_.filtration_order
    order_index = {simplex: index for index, simplex in enumerate(order)}
    reduced_columns: list[list[int]] = [[] for _ in order]
    low_to_column: dict[int, int] = {}
    is_birth = [False] * len(order)
    is_killed = [False] * len(order)
    finite_pairs: list[PersistencePair] = []

    for column_index, sigma in enumerate(order):
        column = sorted(order_index[face] for face in complex_.boundary(sigma))
        if column and column[-1] >= column_index:
            raise RuntimeError("Filtration order does not place faces before cofaces.")

        while column:
            low = column[-1]
            reducer = low_to_column.get(low)
            if reducer is None:
                break
            column = list(_xor_sorted(tuple(column), tuple(reduced_columns[reducer])))

        if not column:
            is_birth[column_index] = True
            continue

        low = column[-1]
        low_to_column[low] = column_index
        reduced_columns[column_index] = column
        is_killed[low] = True

        birth = order[low]
        finite_pairs.append(
            PersistencePair(
                birth=birth,
                death=sigma,
                dimension=complex_.dimension(birth),
                birth_value=complex_.filtration(birth),
                death_value=complex_.filtration(sigma),
            )
        )

    essential = [
        EssentialInterval(
            birth=birth,
            dimension=complex_.dimension(birth),
            birth_value=complex_.filtration(birth),
        )
        for column_index, birth in enumerate(order)
        if is_birth[column_index] and not is_killed[column_index]
    ]
    return PersistenceDiagram(tuple(finite_pairs), tuple(essential))


def compute_standard_persistence_modp(
    complex_: FilteredComplex,
    modulus: int,
) -> PersistenceDiagram:
    """Compute ordinary persistence over the prime field ``F_modulus``."""

    return compute_standard_persistence(complex_, modulus=modulus)


def _compute_standard_persistence_modp(
    complex_: FilteredComplex,
    modulus: int,
) -> PersistenceDiagram:
    order = complex_.filtration_order
    order_index = {simplex: index for index, simplex in enumerate(order)}
    reduced_columns: list[dict[int, int]] = [{} for _ in order]
    low_to_column: dict[int, int] = {}
    is_birth = [False] * len(order)
    is_killed = [False] * len(order)
    finite_pairs: list[PersistencePair] = []

    for column_index, sigma in enumerate(order):
        column = _oriented_boundary_column(complex_, sigma, order_index, modulus)
        low = _modp_low(column)
        if low is not None and low >= column_index:
            raise RuntimeError("Filtration order does not place faces before cofaces.")

        while column:
            low = _modp_low(column)
            if low is None:
                break
            reducer = low_to_column.get(low)
            if reducer is None:
                break
            low_coeff = column[low]
            reducer_coeff = reduced_columns[reducer][low]
            scale = (-low_coeff * _modp_inverse(reducer_coeff, modulus)) % modulus
            _modp_add_scaled(column, reduced_columns[reducer], scale, modulus)

        low = _modp_low(column)
        if low is None:
            is_birth[column_index] = True
            continue

        low_to_column[low] = column_index
        reduced_columns[column_index] = dict(column)
        is_killed[low] = True

        birth = order[low]
        finite_pairs.append(
            PersistencePair(
                birth=birth,
                death=sigma,
                dimension=complex_.dimension(birth),
                birth_value=complex_.filtration(birth),
                death_value=complex_.filtration(sigma),
            )
        )

    essential = [
        EssentialInterval(
            birth=birth,
            dimension=complex_.dimension(birth),
            birth_value=complex_.filtration(birth),
        )
        for column_index, birth in enumerate(order)
        if is_birth[column_index] and not is_killed[column_index]
    ]
    return PersistenceDiagram(tuple(finite_pairs), tuple(essential))


def profile_morse_reference_frame(
    complex_: FilteredComplex,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> MorseReferenceProfile:
    """Profile a Morse-reference frame without running the pivot reduction."""

    algorithm = _normalize_morse_sequence_algorithm(algorithm)
    core_profiler = (
        getattr(_morse_core, "profile_morse_reference_frame_core", None)
        if _morse_core is not None and complex_._cpp is not None
        else None
    )
    if core_profiler is not None:
        return _make_morse_reference_profile(dict(core_profiler(complex_._cpp, algorithm)))

    started = perf_counter()
    frame = compute_morse_sequence_and_reference_map(complex_, algorithm=algorithm)
    payload = _profile_morse_reference_frame_python(
        complex_,
        frame.sequence,
        frame.references,
        algorithm=algorithm,
        profile_seconds=perf_counter() - started,
    )
    return _make_morse_reference_profile(payload)


def profile_morse_sequence_algorithms(
    complex_: FilteredComplex,
    *,
    algorithms: Iterable[str] | str | None = DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
    repeats: int = 1,
) -> tuple[MorseReferenceProfile, ...]:
    """Profile candidate Morse sequence constructors without full reductions."""

    if repeats < 1:
        raise ValueError("repeats must be positive.")

    profiles: list[MorseReferenceProfile] = []
    for algorithm in _normalize_morse_algorithm_portfolio(algorithms):
        best: MorseReferenceProfile | None = None
        for _ in range(repeats):
            candidate = profile_morse_reference_frame(complex_, algorithm=algorithm)
            if best is None or candidate.profile_seconds < best.profile_seconds:
                best = candidate
        if best is None:
            raise RuntimeError("Profile did not run.")
        profiles.append(best)
    return tuple(profiles)


def select_morse_sequence_profile(
    complex_: FilteredComplex,
    *,
    algorithms: Iterable[str] | str | None = DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
    repeats: int = 1,
    selection_metric: str = "estimated_reducer_work",
) -> MorseReferenceProfile:
    """Return the best cheap profile among candidate Morse strategies."""

    candidates = profile_morse_sequence_algorithms(
        complex_,
        algorithms=algorithms,
        repeats=repeats,
    )
    return min(candidates, key=lambda profile: _profile_selection_key(profile, selection_metric))


def compute_persistence_profile_portfolio(
    complex_: FilteredComplex,
    *,
    algorithms: Iterable[str] | str | None = DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
    repeats: int = 1,
    selection_metric: str = "estimated_reducer_work",
    max_critical_ratio: float = 1.0,
) -> AdaptivePersistenceResult:
    if repeats < 1:
        raise ValueError("repeats must be positive.")
    if max_critical_ratio < 0.0 or max_critical_ratio > 1.0:
        raise ValueError("max_critical_ratio must be between 0 and 1.")

    started = perf_counter()
    profiles = profile_morse_sequence_algorithms(
        complex_,
        algorithms=algorithms,
        repeats=repeats,
    )
    selected = min(profiles, key=lambda profile: _profile_selection_key(profile, selection_metric))
    after_selection = perf_counter()

    sequence = compute_morse_sequence(complex_, algorithm=selected.sequence_algorithm)
    after_sequence = perf_counter()

    critical_ratio = selected.critical_ratio
    if critical_ratio > max_critical_ratio:
        diagram = compute_standard_persistence(complex_)
        finished = perf_counter()
        return AdaptivePersistenceResult(
            diagram=diagram,
            method="standard",
            sequence=sequence,
            num_simplices=complex_.size,
            num_critical_simplices=selected.num_critical_simplices,
            critical_ratio=critical_ratio,
            max_critical_ratio=max_critical_ratio,
            sequence_seconds=after_sequence - after_selection,
            persistence_seconds=finished - after_sequence,
            total_seconds=finished - started,
            fallback_reason=(
                f"critical ratio {critical_ratio:.3g} exceeds threshold "
                f"{max_critical_ratio:.3g}"
            ),
            candidate_profiles=profiles,
            selection_metric=selection_metric,
            selection_mode=PROFILE_SELECTION_MODE,
        )

    diagram = compute_morse_persistence(
        complex_,
        sequence,
        algorithm=selected.sequence_algorithm,
    )
    finished = perf_counter()
    return AdaptivePersistenceResult(
        diagram=diagram,
        method="morse",
        sequence=sequence,
        num_simplices=complex_.size,
        num_critical_simplices=selected.num_critical_simplices,
        critical_ratio=critical_ratio,
        max_critical_ratio=max_critical_ratio,
        sequence_seconds=after_sequence - after_selection,
        persistence_seconds=finished - after_sequence,
        total_seconds=finished - started,
        candidate_profiles=profiles,
        selection_metric=selection_metric,
        selection_mode=PROFILE_SELECTION_MODE,
    )


def compute_persistence_portfolio(
    complex_: FilteredComplex,
    *,
    algorithms: Iterable[str] | str | None = DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
    repeats: int = 1,
    frame_mode: str = DEFAULT_MORSE_FRAME_MODE,
    selection_metric: str = "morse_seconds",
    include_standard: bool = True,
    max_critical_ratio: float = 1.0,
) -> AdaptivePersistenceResult:
    if repeats < 1:
        raise ValueError("repeats must be positive.")
    if max_critical_ratio < 0.0 or max_critical_ratio > 1.0:
        raise ValueError("max_critical_ratio must be between 0 and 1.")

    started = perf_counter()
    candidates = benchmark_sequence_algorithms(
        complex_,
        algorithms=algorithms,
        repeats=repeats,
        frame_mode=frame_mode,
        materialize_barcodes=False,
    )
    selected = min(
        candidates,
        key=lambda benchmark: _benchmark_selection_key(benchmark, selection_metric),
    )
    after_selection = perf_counter()

    sequence = compute_morse_sequence(complex_, algorithm=selected.sequence_algorithm)
    after_sequence = perf_counter()

    use_standard = include_standard and selected.standard_seconds <= selected.morse_seconds
    if use_standard:
        diagram = compute_standard_persistence(complex_)
        method = "standard"
        fallback_reason = (
            f"measured standard path {selected.standard_seconds:.3g}s is no slower than "
            f"best Morse path {selected.morse_seconds:.3g}s"
        )
    else:
        diagram = compute_morse_persistence(
            complex_,
            sequence,
            algorithm=selected.sequence_algorithm,
        )
        method = "morse"
        fallback_reason = None

    finished = perf_counter()
    num_simplices = complex_.size
    num_critical = len(sequence.critical_simplices)
    critical_ratio = (float(num_critical) / float(num_simplices)) if num_simplices else 0.0
    return AdaptivePersistenceResult(
        diagram=diagram,
        method=method,
        sequence=sequence,
        num_simplices=num_simplices,
        num_critical_simplices=num_critical,
        critical_ratio=critical_ratio,
        max_critical_ratio=max_critical_ratio,
        sequence_seconds=after_sequence - after_selection,
        persistence_seconds=finished - after_sequence,
        total_seconds=finished - started,
        fallback_reason=fallback_reason,
        candidate_benchmarks=candidates,
        selection_metric=selection_metric,
        selection_mode=BENCHMARK_SELECTION_MODE,
    )


def compute_persistence_adaptive(
    complex_: FilteredComplex,
    *,
    max_critical_ratio: float = 0.4,
    sequence_algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
    candidate_algorithms: Iterable[str] | str | None = None,
    selection_repeats: int = 1,
    selection_metric: str = "estimated_reducer_work",
    selection_mode: str = PROFILE_SELECTION_MODE,
    frame_mode: str = DEFAULT_MORSE_FRAME_MODE,
) -> AdaptivePersistenceResult:
    if max_critical_ratio < 0.0 or max_critical_ratio > 1.0:
        raise ValueError("max_critical_ratio must be between 0 and 1.")

    normalized_request = str(sequence_algorithm).lower().replace("_", "-")
    if candidate_algorithms is not None or normalized_request in {
        AUTO_MORSE_SEQUENCE_ALGORITHM,
        "portfolio",
        "all",
    }:
        normalized_selection_mode = _normalize_morse_selection_mode(selection_mode)
        algorithms = (
            candidate_algorithms
            if candidate_algorithms is not None
            else DEFAULT_MORSE_ALGORITHM_PORTFOLIO
        )
        if normalized_selection_mode == BENCHMARK_SELECTION_MODE:
            benchmark_selection_metric = (
                "morse_seconds"
                if selection_metric == "estimated_reducer_work"
                else selection_metric
            )
            return compute_persistence_portfolio(
                complex_,
                algorithms=algorithms,
                repeats=selection_repeats,
                frame_mode=frame_mode,
                selection_metric=benchmark_selection_metric,
                include_standard=True,
                max_critical_ratio=max_critical_ratio,
            )
        return compute_persistence_profile_portfolio(
            complex_,
            algorithms=algorithms,
            repeats=selection_repeats,
            selection_metric=selection_metric,
            max_critical_ratio=max_critical_ratio,
        )

    sequence_algorithm = _normalize_morse_sequence_algorithm(sequence_algorithm)
    started = perf_counter()
    sequence = compute_morse_sequence(complex_, algorithm=sequence_algorithm)
    after_sequence = perf_counter()

    num_simplices = complex_.size
    num_critical = len(sequence.critical_simplices)
    critical_ratio = (float(num_critical) / float(num_simplices)) if num_simplices else 0.0

    if critical_ratio > max_critical_ratio:
        diagram = compute_standard_persistence(complex_)
        finished = perf_counter()
        return AdaptivePersistenceResult(
            diagram=diagram,
            method="standard",
            sequence=sequence,
            num_simplices=num_simplices,
            num_critical_simplices=num_critical,
            critical_ratio=critical_ratio,
            max_critical_ratio=max_critical_ratio,
            sequence_seconds=after_sequence - started,
            persistence_seconds=finished - after_sequence,
            total_seconds=finished - started,
            fallback_reason=(
                f"critical ratio {critical_ratio:.3g} exceeds threshold "
                f"{max_critical_ratio:.3g}"
            ),
        )

    diagram = compute_morse_persistence(complex_, sequence, algorithm=sequence_algorithm)
    finished = perf_counter()
    return AdaptivePersistenceResult(
        diagram=diagram,
        method="morse",
        sequence=sequence,
        num_simplices=num_simplices,
        num_critical_simplices=num_critical,
        critical_ratio=critical_ratio,
        max_critical_ratio=max_critical_ratio,
        sequence_seconds=after_sequence - started,
        persistence_seconds=finished - after_sequence,
        total_seconds=finished - started,
    )


def gudhi_barcode(
    complex_: FilteredComplex,
    *,
    include_zero: bool = False,
    modulus: int = 2,
) -> tuple[tuple[tuple[int, float, float], ...], tuple[tuple[int, float], ...]]:
    """Compute a comparable barcode with an installed/importable GUDHI."""

    modulus = _validate_prime_modulus(modulus)

    try:
        import gudhi  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError("GUDHI is not importable in this Python environment.") from exc

    simplex_tree = gudhi.SimplexTree()
    for simplex_id in complex_.filtration_order:
        simplex_tree.insert(list(complex_.vertices(simplex_id)), filtration=complex_.filtration(simplex_id))

    persistence = simplex_tree.persistence(
        homology_coeff_field=modulus,
        min_persistence=0,
        persistence_dim_max=True,
    )
    finite: list[tuple[int, float, float]] = []
    essential: list[tuple[int, float]] = []
    for dimension, (birth, death) in persistence:
        if death == inf:
            essential.append((int(dimension), float(birth)))
        elif include_zero or birth < death:
            finite.append((int(dimension), float(birth), float(death)))
    return tuple(sorted(finite)), tuple(sorted(essential))


def gudhi_cubical_barcode(
    grid: CubicalGrid2DComplex,
    *,
    include_zero: bool = False,
    modulus: int = 2,
) -> tuple[tuple[tuple[int, float, float], ...], tuple[tuple[int, float], ...]]:
    """Compute a comparable barcode with GUDHI's 2D cubical complex."""

    modulus = _validate_prime_modulus(modulus)

    try:
        import gudhi  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError("GUDHI is not importable in this Python environment.") from exc

    cubical_complex = getattr(gudhi, "CubicalComplex", None)
    if cubical_complex is None:
        raise RuntimeError("This GUDHI installation does not expose CubicalComplex.")

    vertices = [
        [grid.filtration(grid.vertex(x, y)) for x in range(grid.vertex_width)]
        for y in range(grid.vertex_height)
    ]
    cubical = cubical_complex(vertices=vertices)
    persistence = cubical.persistence(homology_coeff_field=modulus, min_persistence=0)

    finite: list[tuple[int, float, float]] = []
    essential: list[tuple[int, float]] = []
    for dimension, (birth, death) in persistence:
        if death == inf:
            essential.append((int(dimension), float(birth)))
        elif include_zero or birth < death:
            finite.append((int(dimension), float(birth), float(death)))
    return tuple(sorted(finite)), tuple(sorted(essential))


def assert_matches_standard(complex_: FilteredComplex, *, include_zero: bool = False) -> None:
    morse = compute_morse_persistence(complex_).finite_barcode(include_zero=include_zero)
    standard = compute_standard_persistence(complex_).finite_barcode(include_zero=include_zero)
    if morse != standard:
        raise AssertionError(f"Morse finite barcode {morse!r} != standard {standard!r}")

    morse_essential = compute_morse_persistence(complex_).essential_barcode()
    standard_essential = compute_standard_persistence(complex_).essential_barcode()
    if morse_essential != standard_essential:
        raise AssertionError(
            f"Morse essential barcode {morse_essential!r} != standard {standard_essential!r}"
        )


def assert_matches_gudhi(complex_: FilteredComplex, *, include_zero: bool = False) -> None:
    morse = compute_morse_persistence(complex_)
    gudhi_finite, gudhi_essential = gudhi_barcode(complex_, include_zero=include_zero)
    if morse.finite_barcode(include_zero=include_zero) != gudhi_finite:
        raise AssertionError(
            f"Morse finite barcode {morse.finite_barcode(include_zero=include_zero)!r} "
            f"!= GUDHI {gudhi_finite!r}"
        )
    if morse.essential_barcode() != gudhi_essential:
        raise AssertionError(
            f"Morse essential barcode {morse.essential_barcode()!r} != GUDHI {gudhi_essential!r}"
        )


def assert_matches_gudhi_cubical(
    grid: CubicalGrid2DComplex,
    *,
    include_zero: bool = False,
    modulus: int = 2,
) -> None:
    morse = compute_morse_persistence(grid, modulus=modulus)
    gudhi_finite, gudhi_essential = gudhi_cubical_barcode(
        grid,
        include_zero=include_zero,
        modulus=modulus,
    )
    if morse.finite_barcode(include_zero=include_zero) != gudhi_finite:
        raise AssertionError(
            f"Morse finite barcode {morse.finite_barcode(include_zero=include_zero)!r} "
            f"!= GUDHI cubical {gudhi_finite!r}"
        )
    if morse.essential_barcode() != gudhi_essential:
        raise AssertionError(
            f"Morse essential barcode {morse.essential_barcode()!r} "
            f"!= GUDHI cubical {gudhi_essential!r}"
        )


def cpp_backend_available() -> bool:
    return _morse_core is not None


def cpp_backend_active(complex_: FilteredComplex) -> bool:
    return complex_.cpp_backend_active()


def cpp_filtered_complex_from_simplices(
    simplices: Iterable[tuple[Sequence[int], float]], *, finalize: bool = True
) -> object:
    _require_cpp_backend()
    complex_ = _morse_core.FilteredComplex()
    for vertices, filtration in simplices:
        complex_.add_simplex(list(_canonical_simplex(vertices)), float(filtration))
    if finalize:
        complex_.finalize()
    return complex_


def cpp_cubical_grid_2d_from_vertex_values(
    vertex_width: int,
    vertex_height: int,
    vertex_values: Sequence[float],
) -> object:
    _require_cpp_backend()
    if CppCubicalGrid2DComplex is None:
        raise RuntimeError("The C++ backend does not expose 2D cubical grids.")
    constructor = getattr(CppCubicalGrid2DComplex, "from_vertex_values", None)
    values = [float(value) for value in vertex_values]
    if constructor is not None:
        return constructor(int(vertex_width), int(vertex_height), values)
    return CppCubicalGrid2DComplex(int(vertex_width), int(vertex_height), values)


def cpp_compute_morse_sequence(
    cpp_complex: object,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> object:
    _require_cpp_backend()
    return _core_compute_morse_sequence(
        cpp_complex,
        _normalize_morse_sequence_algorithm(algorithm),
    )


def cpp_compute_reference_map(cpp_complex: object, cpp_sequence: object) -> tuple[Annotation, ...]:
    _require_cpp_backend()
    return _references_from_cpp_annotations(
        _morse_core.compute_reference_map(cpp_complex, cpp_sequence)
    )


def cpp_compute_reference_map_modp(
    cpp_complex: object,
    cpp_sequence: object,
    modulus: int,
) -> tuple[FieldAnnotation, ...]:
    _require_cpp_backend()
    runner = getattr(_morse_core, "compute_reference_map_modp", None)
    if runner is None:
        raise RuntimeError("The C++ backend does not expose prime-field reference maps.")
    return _field_references_from_cpp_annotations(
        runner(cpp_complex, cpp_sequence, _validate_prime_modulus(modulus))
    )


def cpp_compute_reference_map_object(cpp_complex: object, cpp_sequence: object) -> object:
    _require_cpp_backend()
    return _morse_core.compute_reference_map_object(cpp_complex, cpp_sequence)


def cpp_compute_coreference_map(cpp_complex: object, cpp_sequence: object) -> tuple[Annotation, ...]:
    _require_cpp_backend()
    return _references_from_cpp_annotations(
        _morse_core.compute_coreference_map(cpp_complex, cpp_sequence)
    )


def cpp_compute_coreference_map_modp(
    cpp_complex: object,
    cpp_sequence: object,
    modulus: int,
) -> tuple[FieldAnnotation, ...]:
    _require_cpp_backend()
    runner = getattr(_morse_core, "compute_coreference_map_modp", None)
    if runner is None:
        raise RuntimeError("The C++ backend does not expose prime-field coreference maps.")
    return _field_references_from_cpp_annotations(
        runner(cpp_complex, cpp_sequence, _validate_prime_modulus(modulus))
    )


def cpp_compute_coreference_map_object(cpp_complex: object, cpp_sequence: object) -> object:
    _require_cpp_backend()
    return _morse_core.compute_coreference_map_object(cpp_complex, cpp_sequence)


def cpp_compute_morse_sequence_and_reference_map_object(
    cpp_complex: object,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> object:
    _require_cpp_backend()
    builder = getattr(_morse_core, "compute_morse_sequence_and_reference_map_object", None)
    if builder is None:
        raise RuntimeError("The C++ backend does not expose fused Morse reference frames.")
    return builder(cpp_complex, _normalize_morse_sequence_algorithm(algorithm))


def cpp_compute_morse_sequence_and_coreference_map_object(
    cpp_complex: object,
    *,
    algorithm: str = COREDUCTION_SEQUENCE,
) -> object:
    _require_cpp_backend()
    builder = getattr(_morse_core, "compute_morse_sequence_and_coreference_map_object", None)
    if builder is None:
        raise RuntimeError("The C++ backend does not expose fused Morse coreference frames.")
    return builder(cpp_complex, _normalize_morse_sequence_algorithm(algorithm))


def cpp_profile_morse_reference_frame(
    cpp_complex: object,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> MorseReferenceProfile:
    _require_cpp_backend()
    profiler = getattr(_morse_core, "profile_morse_reference_frame_core", None)
    if profiler is None:
        raise RuntimeError("The C++ backend does not expose Morse reference profiling.")
    return _make_morse_reference_profile(
        dict(profiler(cpp_complex, _normalize_morse_sequence_algorithm(algorithm)))
    )


def cpp_reference_map_to_tuple(cpp_references: object) -> tuple[Annotation, ...]:
    _require_cpp_backend()
    return _references_from_cpp_annotations(cpp_references.annotations())


def cpp_compute_morse_persistence(
    cpp_complex: object,
    cpp_sequence: object | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> PersistenceDiagram:
    _require_cpp_backend()
    if cpp_sequence is None:
        cpp_sequence = _core_compute_morse_sequence(
            cpp_complex,
            _normalize_morse_sequence_algorithm(algorithm),
        )
    return _diagram_from_cpp_dict(_morse_core.compute_morse_persistence(cpp_complex, cpp_sequence))


def cpp_compute_morse_persistence_modp(
    cpp_complex: object,
    cpp_sequence: object | None = None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
    modulus: int,
) -> PersistenceDiagram:
    _require_cpp_backend()
    runner = getattr(_morse_core, "compute_morse_persistence_modp", None)
    if runner is None:
        raise RuntimeError("The C++ backend does not expose prime-field Morse persistence.")
    if cpp_sequence is None:
        cpp_sequence = _core_compute_morse_sequence(
            cpp_complex,
            _normalize_morse_sequence_algorithm(algorithm),
        )
    return _diagram_from_cpp_dict(runner(cpp_complex, cpp_sequence, _validate_prime_modulus(modulus)))


def cpp_compute_morse_coreference_persistence(
    cpp_complex: object,
    cpp_sequence: object | None = None,
    *,
    algorithm: str = COREDUCTION_SEQUENCE,
) -> PersistenceDiagram:
    _require_cpp_backend()
    if cpp_sequence is None:
        cpp_sequence = _core_compute_morse_sequence(
            cpp_complex,
            _normalize_morse_sequence_algorithm(algorithm),
        )
    return _diagram_from_cpp_dict(
        _morse_core.compute_morse_coreference_persistence(cpp_complex, cpp_sequence)
    )


def cpp_compute_morse_coreference_persistence_modp(
    cpp_complex: object,
    cpp_sequence: object | None = None,
    *,
    algorithm: str = COREDUCTION_SEQUENCE,
    modulus: int,
) -> PersistenceDiagram:
    _require_cpp_backend()
    runner = getattr(_morse_core, "compute_morse_coreference_persistence_modp", None)
    if runner is None:
        raise RuntimeError("The C++ backend does not expose prime-field Morse coreference persistence.")
    if cpp_sequence is None:
        cpp_sequence = _core_compute_morse_sequence(
            cpp_complex,
            _normalize_morse_sequence_algorithm(algorithm),
        )
    return _diagram_from_cpp_dict(runner(cpp_complex, cpp_sequence, _validate_prime_modulus(modulus)))


def cpp_reduce_morse_persistence(
    cpp_complex: object,
    cpp_sequence: object,
    references: Sequence[Annotation] | object,
) -> PersistenceDiagram:
    _require_cpp_backend()
    if CppReferenceMap is not None and isinstance(references, CppReferenceMap):
        return _diagram_from_cpp_dict(
            _morse_core.reduce_morse_persistence_object(cpp_complex, cpp_sequence, references)
        )
    return _diagram_from_cpp_dict(
        _morse_core.reduce_morse_persistence(
            cpp_complex,
            cpp_sequence,
            [list(annotation) for annotation in references],
        )
    )


def cpp_reduce_morse_coreference_persistence(
    cpp_complex: object,
    cpp_sequence: object,
    coreferences: object,
) -> PersistenceDiagram:
    _require_cpp_backend()
    if CppReferenceMap is not None and isinstance(coreferences, CppReferenceMap):
        return _diagram_from_cpp_dict(
            _morse_core.reduce_morse_coreference_persistence_object(
                cpp_complex,
                cpp_sequence,
                coreferences,
            )
        )
    raise TypeError("C++ coreference reduction expects a C++ coreference map object.")


def cpp_reduce_morse_persistence_with_metrics(
    cpp_complex: object,
    cpp_sequence: object,
    references: Sequence[Annotation] | object,
) -> tuple[PersistenceDiagram, dict[str, object]]:
    _require_cpp_backend()
    if CppReferenceMap is not None and isinstance(references, CppReferenceMap):
        result = _morse_core.reduce_morse_persistence_object_with_metrics(
            cpp_complex, cpp_sequence, references
        )
    else:
        result = _morse_core.reduce_morse_persistence_with_metrics(
            cpp_complex,
            cpp_sequence,
            [list(annotation) for annotation in references],
        )
    return _diagram_from_cpp_dict(result["diagram"]), result["metrics"]


def cpp_compute_standard_persistence(cpp_complex: object) -> PersistenceDiagram:
    _require_cpp_backend()
    return _diagram_from_cpp_dict(_morse_core.compute_standard_persistence(cpp_complex))


def cpp_compute_standard_persistence_modp(cpp_complex: object, modulus: int) -> PersistenceDiagram:
    _require_cpp_backend()
    runner = getattr(_morse_core, "compute_standard_persistence_modp", None)
    if runner is None:
        raise RuntimeError("The C++ backend does not expose prime-field standard persistence.")
    return _diagram_from_cpp_dict(runner(cpp_complex, _validate_prime_modulus(modulus)))


def _coreduction_direction_benchmark_from_payload(
    payload: dict[str, object],
    *,
    repeats: int,
    standard_seconds: float | None = None,
) -> CoreductionDirectionBenchmark:
    metrics = payload["metrics"]  # type: ignore[assignment]
    inverse_metrics = metrics["inverse_store"]  # type: ignore[index]
    measured_standard_seconds = (
        1.0e-9 * int(payload["standard_nanoseconds"])
        if standard_seconds is None
        else standard_seconds
    )
    return CoreductionDirectionBenchmark(
        direction=str(payload["direction"]),
        num_simplices=int(payload["num_simplices"]),
        num_levels=int(payload["num_levels"]),
        num_critical_simplices=int(payload["num_critical_simplices"]),
        frame_seconds=1.0e-9 * int(payload["frame_nanoseconds"]),
        morse_reduction_seconds=1.0e-9 * int(payload["morse_reduction_nanoseconds"]),
        reducer_setup_seconds=1.0e-9 * int(metrics["reducer_setup_nanoseconds"]),  # type: ignore[index]
        reducer_compute_seconds=1.0e-9 * int(metrics["reducer_compute_nanoseconds"]),  # type: ignore[index]
        morse_seconds=1.0e-9 * int(payload["morse_nanoseconds"]),
        standard_seconds=measured_standard_seconds,
        repeats=repeats,
        finite_interval_count=int(payload["morse_finite_interval_count"]),
        essential_interval_count=int(payload["morse_essential_interval_count"]),
        reducer_working_set_size=int(metrics["working_set_size"]),  # type: ignore[index]
        reducer_initial_total_annotation_size=int(
            inverse_metrics["initial_total_annotation_size"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_xors=int(metrics["boundary_annotation_xors"]),  # type: ignore[index]
        reducer_boundary_annotation_total_input_size=int(
            metrics["boundary_annotation_total_input_size"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_total_output_size=int(
            metrics["boundary_annotation_total_output_size"]  # type: ignore[index]
        ),
        reducer_pivot_eliminations=int(metrics["pivot_eliminations"]),  # type: ignore[index]
        reducer_xor_applied=int(inverse_metrics["xor_applied"]),  # type: ignore[index]
        reducer_xor_total_input_size=int(
            inverse_metrics["xor_total_input_size"]  # type: ignore[index]
        ),
        reducer_xor_total_output_size=int(
            inverse_metrics["xor_total_output_size"]  # type: ignore[index]
        ),
        metrics=metrics,  # type: ignore[arg-type]
    )


def benchmark_coreduction_directions(
    complex_: FilteredComplex,
    *,
    repeats: int = 3,
) -> tuple[CoreductionDirectionBenchmark, ...]:
    if repeats < 1:
        raise ValueError("repeats must be positive.")

    core_runner = (
        getattr(_morse_core, "benchmark_coreduction_directions_core", None)
        if _morse_core is not None and complex_._cpp is not None
        else None
    )
    if core_runner is not None:
        best_payloads: dict[str, dict[str, object]] = {}
        best_morse_seconds: dict[str, float] = {}
        best_standard_seconds = inf
        for _ in range(repeats):
            for payload in core_runner(complex_._cpp):
                payload = dict(payload)
                if (
                    not payload["finite_signature_matches"]
                    or not payload["essential_signature_matches"]
                ):
                    raise AssertionError(
                        "Same-level reduction direction barcode signature != standard signature"
                    )
                direction = str(payload["direction"])
                morse_seconds = 1.0e-9 * int(payload["morse_nanoseconds"])
                standard_seconds = 1.0e-9 * int(payload["standard_nanoseconds"])
                best_standard_seconds = min(best_standard_seconds, standard_seconds)
                if direction not in best_morse_seconds or morse_seconds < best_morse_seconds[direction]:
                    best_morse_seconds[direction] = morse_seconds
                    best_payloads[direction] = payload

        return tuple(
            _coreduction_direction_benchmark_from_payload(
                best_payloads[direction],
                repeats=repeats,
                standard_seconds=best_standard_seconds,
            )
            for direction in sorted(best_payloads)
        )

    rows: list[CoreductionDirectionBenchmark] = []
    standard_best = inf
    standard_diagram: PersistenceDiagram | None = None
    for _ in range(repeats):
        started = perf_counter()
        standard_diagram = compute_standard_persistence(complex_)
        standard_best = min(standard_best, perf_counter() - started)

    if standard_diagram is None:
        raise RuntimeError("Benchmark did not run.")

    for direction in ("reference", "coreference"):
        best_total = inf
        best_payload: dict[str, object] | None = None
        for _ in range(repeats):
            started = perf_counter()
            if direction == "reference":
                frame = compute_morse_sequence_and_reference_map(
                    complex_,
                    algorithm=COREDUCTION_SEQUENCE,
                )
                after_frame = perf_counter()
                diagram = compute_morse_persistence(
                    complex_,
                    frame.sequence,
                    frame.references,
                    algorithm=COREDUCTION_SEQUENCE,
                )
            else:
                frame = compute_morse_sequence_and_coreference_map(
                    complex_,
                    algorithm=COREDUCTION_SEQUENCE,
                )
                after_frame = perf_counter()
                diagram = compute_morse_coreference_persistence(
                    complex_,
                    frame.sequence,
                    frame.coreferences,
                    algorithm=COREDUCTION_SEQUENCE,
                )
            finished = perf_counter()
            if diagram.finite_barcode() != standard_diagram.finite_barcode():
                raise AssertionError(
                    "Same-level reduction direction finite barcode != standard barcode"
                )
            if diagram.essential_barcode() != standard_diagram.essential_barcode():
                raise AssertionError(
                    "Same-level reduction direction essential barcode != standard barcode"
                )
            total = finished - started
            if total < best_total:
                metrics = _empty_reducer_metrics()
                best_total = total
                best_payload = {
                    "direction": direction,
                    "num_simplices": complex_.size,
                    "num_levels": complex_.num_levels,
                    "num_critical_simplices": len(frame.sequence.critical_simplices),
                    "frame_nanoseconds": int(1.0e9 * (after_frame - started)),
                    "morse_reduction_nanoseconds": int(1.0e9 * (finished - after_frame)),
                    "morse_nanoseconds": int(1.0e9 * total),
                    "standard_nanoseconds": int(1.0e9 * standard_best),
                    "morse_finite_interval_count": len(diagram.finite_barcode()),
                    "morse_essential_interval_count": len(diagram.essential_barcode()),
                    "metrics": metrics,
                }
        if best_payload is None:
            raise RuntimeError("Benchmark did not run.")
        rows.append(
            _coreduction_direction_benchmark_from_payload(
                best_payload,
                repeats=repeats,
                standard_seconds=standard_best,
            )
        )
    return tuple(rows)


SameLevelReductionDirectionBenchmark = CoreductionDirectionBenchmark


def benchmark_same_level_reduction_directions(
    complex_: FilteredComplex,
    *,
    repeats: int = 3,
) -> tuple[CoreductionDirectionBenchmark, ...]:
    """Benchmark reference/coreference directions for same-level reduction."""

    return benchmark_coreduction_directions(complex_, repeats=repeats)


def benchmark_persistence(
    complex_: FilteredComplex,
    *,
    repeats: int = 3,
    sequence_algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
    frame_mode: str = DEFAULT_MORSE_FRAME_MODE,
    materialize_barcodes: bool = True,
) -> PersistenceBenchmark:
    if repeats < 1:
        raise ValueError("repeats must be positive.")
    sequence_algorithm = _normalize_morse_sequence_algorithm(sequence_algorithm)
    frame_mode = _normalize_morse_frame_mode(frame_mode)

    morse_total_best = inf
    sequence_best = inf
    reference_best = inf
    morse_reduction_best = inf
    standard_best = inf
    morse_diagram: PersistenceDiagram | None = None
    standard_diagram: PersistenceDiagram | None = None
    sequence: MorseSequence | None = None
    reducer_metrics: dict[str, object] = _empty_reducer_metrics()
    measured_frame_mode = frame_mode

    if not materialize_barcodes and _morse_core is not None and complex_._cpp is not None:
        core_runner = getattr(_morse_core, "benchmark_morse_reference_core", None)
        if core_runner is not None:
            best_core_result: dict[str, object] | None = None
            best_morse_seconds = inf
            best_standard_seconds = inf
            for _ in range(repeats):
                result = core_runner(complex_._cpp, sequence_algorithm, frame_mode)
                if not result["finite_signature_matches"] or not result["essential_signature_matches"]:
                    raise AssertionError("C++ core Morse barcode signature != standard signature")
                morse_seconds = 1.0e-9 * int(result["morse_nanoseconds"])
                standard_seconds = 1.0e-9 * int(result["standard_nanoseconds"])
                if morse_seconds < best_morse_seconds:
                    best_morse_seconds = morse_seconds
                    best_core_result = result
                best_standard_seconds = min(best_standard_seconds, standard_seconds)

            if best_core_result is None:
                raise RuntimeError("Benchmark did not run.")

            reducer_metrics = best_core_result["metrics"]  # type: ignore[assignment]
            return _make_persistence_benchmark(
                complex_,
                num_critical_simplices=int(best_core_result["num_critical_simplices"]),
                sequence_algorithm=sequence_algorithm,
                frame_mode=str(best_core_result["frame_mode"]),
                frame_metrics=best_core_result["frame_metrics"],  # type: ignore[arg-type]
                reducer_metrics=reducer_metrics,
                sequence_seconds=1.0e-9 * int(best_core_result["sequence_nanoseconds"]),
                reference_seconds=1.0e-9 * int(best_core_result["reference_nanoseconds"]),
                morse_reduction_seconds=1.0e-9
                * int(best_core_result["morse_reduction_nanoseconds"]),
                morse_seconds=best_morse_seconds,
                standard_seconds=best_standard_seconds,
                repeats=repeats,
                finite_interval_count=int(best_core_result["morse_finite_interval_count"]),
                essential_interval_count=int(best_core_result["morse_essential_interval_count"]),
                barcodes_materialized=False,
                validation_mode="core",
                finite_barcode=(),
                essential_barcode=(),
            )

    if _morse_core is not None and complex_._cpp is not None:
        best_cpp_sequence = None
        best_cpp_frame = None
        for _ in range(repeats):
            started = perf_counter()
            cpp_frame = (
                _compute_cpp_reference_frame(complex_, sequence_algorithm)
                if frame_mode == FUSED_FRAME
                else None
            )
            if cpp_frame is None:
                cpp_sequence = _core_compute_morse_sequence(complex_._cpp, sequence_algorithm)
                after_sequence = perf_counter()
                references = _morse_core.compute_reference_map_object(complex_._cpp, cpp_sequence)
                after_reference = perf_counter()
                sequence_duration = after_sequence - started
                reference_duration = after_reference - after_sequence
                current_frame_mode = SEPARATE_FRAME
            else:
                after_reference = perf_counter()
                sequence_duration = after_reference - started
                reference_duration = 0.0
                current_frame_mode = FUSED_FRAME
                reduction_result = _morse_core.reduce_morse_reference_frame_object_with_metrics(
                    complex_._cpp,
                    cpp_frame,
                )
                cpp_sequence = None
            if cpp_frame is None:
                reduction_result = _morse_core.reduce_morse_persistence_object_with_metrics(
                    complex_._cpp, cpp_sequence, references
                )
            morse_diagram = _diagram_from_cpp_dict(reduction_result["diagram"])
            finished = perf_counter()

            total = finished - started
            if total < morse_total_best:
                morse_total_best = total
                sequence_best = sequence_duration
                reference_best = reference_duration
                morse_reduction_best = finished - after_reference
                reducer_metrics = reduction_result["metrics"]
                best_cpp_sequence = cpp_sequence
                best_cpp_frame = cpp_frame
                measured_frame_mode = current_frame_mode

        if best_cpp_frame is not None:
            sequence = _sequence_from_cpp_sequence(
                complex_,
                best_cpp_frame.sequence,
                algorithm=sequence_algorithm,
            )
        elif best_cpp_sequence is not None:
            sequence = _sequence_from_cpp_sequence(
                complex_,
                best_cpp_sequence,
                algorithm=sequence_algorithm,
            )
    else:
        for _ in range(repeats):
            started = perf_counter()
            if frame_mode == FUSED_FRAME:
                frame = compute_morse_sequence_and_reference_map(
                    complex_,
                    algorithm=sequence_algorithm,
                )
                after_reference = perf_counter()
                sequence = frame.sequence
                references = frame.references
                sequence_duration = after_reference - started
                reference_duration = 0.0
                current_frame_mode = FUSED_FRAME
            else:
                sequence = compute_morse_sequence(complex_, algorithm=sequence_algorithm)
                after_sequence = perf_counter()
                references = compute_reference_map(complex_, sequence, algorithm=sequence_algorithm)
                after_reference = perf_counter()
                sequence_duration = after_sequence - started
                reference_duration = after_reference - after_sequence
                current_frame_mode = SEPARATE_FRAME
            morse_diagram = compute_morse_persistence(complex_, sequence, references)
            finished = perf_counter()

            total = finished - started
            if total < morse_total_best:
                morse_total_best = total
                sequence_best = sequence_duration
                reference_best = reference_duration
                morse_reduction_best = finished - after_reference
                measured_frame_mode = current_frame_mode

    for _ in range(repeats):
        started = perf_counter()
        standard_diagram = compute_standard_persistence(complex_)
        standard_best = min(standard_best, perf_counter() - started)

    if morse_diagram is None or standard_diagram is None or sequence is None:
        raise RuntimeError("Benchmark did not run.")

    if morse_diagram.finite_barcode() != standard_diagram.finite_barcode():
        raise AssertionError(
            f"Morse finite barcode {morse_diagram.finite_barcode()!r} "
            f"!= standard {standard_diagram.finite_barcode()!r}"
        )
    if morse_diagram.essential_barcode() != standard_diagram.essential_barcode():
        raise AssertionError(
            f"Morse essential barcode {morse_diagram.essential_barcode()!r} "
            f"!= standard {standard_diagram.essential_barcode()!r}"
        )

    finite_barcode = morse_diagram.finite_barcode()
    essential_barcode = morse_diagram.essential_barcode()
    return _make_persistence_benchmark(
        complex_,
        num_critical_simplices=len(sequence.critical_simplices),
        sequence_algorithm=sequence_algorithm,
        frame_mode=measured_frame_mode,
        frame_metrics=_empty_frame_metrics(),
        reducer_metrics=reducer_metrics,
        sequence_seconds=sequence_best,
        reference_seconds=reference_best,
        morse_reduction_seconds=morse_reduction_best,
        morse_seconds=morse_total_best,
        standard_seconds=standard_best,
        repeats=repeats,
        finite_interval_count=len(finite_barcode),
        essential_interval_count=len(essential_barcode),
        barcodes_materialized=True,
        validation_mode="materialized",
        finite_barcode=finite_barcode,
        essential_barcode=essential_barcode,
    )


def _benchmark_selection_key(
    benchmark: PersistenceBenchmark,
    selection_metric: str,
) -> tuple[float, float, float]:
    metric = selection_metric.lower().replace("_", "-")
    critical_ratio = (
        float(benchmark.num_critical_simplices) / float(benchmark.num_simplices)
        if benchmark.num_simplices
        else 0.0
    )
    reducer_work = float(
        benchmark.reducer_initial_total_annotation_size
        + benchmark.reducer_boundary_annotation_total_input_size
        + benchmark.reducer_boundary_annotation_total_output_size
        + benchmark.reducer_xor_total_input_size
        + benchmark.reducer_xor_total_output_size
        + benchmark.reducer_remove_total_annotation_size
    )
    if metric in {"morse-seconds", "morse-seconds-best", "time"}:
        return (benchmark.morse_seconds, reducer_work, critical_ratio)
    if metric in {"critical-ratio", "critical-count", "critical"}:
        return (critical_ratio, benchmark.morse_seconds, reducer_work)
    if metric in {"reducer-work", "annotation-work", "work"}:
        return (reducer_work, benchmark.morse_seconds, critical_ratio)
    supported = ", ".join(MORSE_ALGORITHM_SELECTION_METRICS)
    raise ValueError(f"Unknown algorithm selection metric {selection_metric!r}. Supported: {supported}.")


def benchmark_sequence_algorithms(
    complex_: FilteredComplex,
    *,
    algorithms: Iterable[str] | str | None = DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
    repeats: int = 1,
    frame_mode: str = DEFAULT_MORSE_FRAME_MODE,
    materialize_barcodes: bool = False,
) -> tuple[PersistenceBenchmark, ...]:
    """Benchmark a portfolio of Morse sequence constructors on one complex."""

    return tuple(
        benchmark_persistence(
            complex_,
            repeats=repeats,
            sequence_algorithm=algorithm,
            frame_mode=frame_mode,
            materialize_barcodes=materialize_barcodes,
        )
        for algorithm in _normalize_morse_algorithm_portfolio(algorithms)
    )


def select_morse_sequence_algorithm(
    complex_: FilteredComplex,
    *,
    algorithms: Iterable[str] | str | None = DEFAULT_MORSE_ALGORITHM_PORTFOLIO,
    repeats: int = 1,
    frame_mode: str = DEFAULT_MORSE_FRAME_MODE,
    selection_metric: str = "morse_seconds",
    materialize_barcodes: bool = False,
) -> PersistenceBenchmark:
    """Return the best measured Morse strategy for one complex."""

    candidates = benchmark_sequence_algorithms(
        complex_,
        algorithms=algorithms,
        repeats=repeats,
        frame_mode=frame_mode,
        materialize_barcodes=materialize_barcodes,
    )
    return min(candidates, key=lambda benchmark: _benchmark_selection_key(benchmark, selection_metric))


def _make_persistence_benchmark(
    complex_: FilteredComplex,
    *,
    num_critical_simplices: int,
    sequence_algorithm: str,
    frame_mode: str,
    frame_metrics: dict[str, object],
    reducer_metrics: dict[str, object],
    sequence_seconds: float,
    reference_seconds: float,
    morse_reduction_seconds: float,
    morse_seconds: float,
    standard_seconds: float,
    repeats: int,
    finite_interval_count: int,
    essential_interval_count: int,
    barcodes_materialized: bool,
    validation_mode: str,
    finite_barcode: tuple[tuple[int, float, float], ...],
    essential_barcode: tuple[tuple[int, float], ...],
) -> PersistenceBenchmark:
    inverse_metrics = reducer_metrics["inverse_store"]
    return PersistenceBenchmark(
        num_simplices=complex_.size,
        num_levels=complex_.num_levels,
        num_critical_simplices=num_critical_simplices,
        sequence_algorithm=sequence_algorithm,
        frame_mode=frame_mode,
        reference_final_live_nonempty_annotations=int(
            frame_metrics["final_live_nonempty_annotations"]
        ),
        reference_final_live_total_annotation_size=int(
            frame_metrics["final_live_total_annotation_size"]
        ),
        reference_peak_live_nonempty_annotations=int(
            frame_metrics["peak_live_nonempty_annotations"]
        ),
        reference_peak_live_total_annotation_size=int(
            frame_metrics["peak_live_total_annotation_size"]
        ),
        reference_released_annotations=int(frame_metrics["released_annotations"]),
        reference_released_total_annotation_size=int(
            frame_metrics["released_total_annotation_size"]
        ),
        reducer_working_set_size=int(reducer_metrics["working_set_size"]),
        reducer_initial_nonempty_annotations=int(
            inverse_metrics["initial_nonempty_annotations"]  # type: ignore[index]
        ),
        reducer_initial_total_annotation_size=int(
            inverse_metrics["initial_total_annotation_size"]  # type: ignore[index]
        ),
        reducer_initial_max_annotation_size=int(
            inverse_metrics["initial_max_annotation_size"]  # type: ignore[index]
        ),
        reducer_initial_inverse_list_entries=int(
            inverse_metrics["initial_inverse_list_entries"]  # type: ignore[index]
        ),
        reducer_boundary_plan_face_scans=int(reducer_metrics["boundary_plan_face_scans"]),
        reducer_boundary_annotation_candidate_criticals=int(
            reducer_metrics["boundary_annotation_candidate_criticals"]
        ),
        reducer_boundary_annotation_zero_skipped_criticals=int(
            reducer_metrics["boundary_annotation_zero_skipped_criticals"]
        ),
        reducer_boundary_annotation_zero_skipped_faces=int(
            reducer_metrics["boundary_annotation_zero_skipped_faces"]
        ),
        reducer_boundary_annotation_xors=int(reducer_metrics["boundary_annotation_xors"]),
        reducer_boundary_annotation_total_input_size=int(
            reducer_metrics["boundary_annotation_total_input_size"]
        ),
        reducer_boundary_annotation_total_output_size=int(
            reducer_metrics["boundary_annotation_total_output_size"]
        ),
        reducer_boundary_annotation_max_size=int(
            reducer_metrics["boundary_annotation_max_size"]
        ),
        reducer_boundary_annotation_max_output_size=int(
            reducer_metrics["boundary_annotation_max_output_size"]
        ),
        reducer_pivot_eliminations=int(reducer_metrics["pivot_eliminations"]),
        reducer_remove_candidate_scans=int(
            inverse_metrics["remove_candidate_scans"]  # type: ignore[index]
        ),
        reducer_remove_applied=int(
            inverse_metrics["remove_applied"]  # type: ignore[index]
        ),
        reducer_remove_total_annotation_size=int(
            inverse_metrics["remove_total_annotation_size"]  # type: ignore[index]
        ),
        reducer_remove_max_annotation_size=int(
            inverse_metrics["remove_max_annotation_size"]  # type: ignore[index]
        ),
        reducer_xor_candidate_scans=int(
            inverse_metrics["xor_candidate_scans"]  # type: ignore[index]
        ),
        reducer_xor_applied=int(
            inverse_metrics["xor_applied"]  # type: ignore[index]
        ),
        reducer_xor_changed_labels=int(
            inverse_metrics["xor_changed_labels"]  # type: ignore[index]
        ),
        reducer_xor_total_input_size=int(
            inverse_metrics["xor_total_input_size"]  # type: ignore[index]
        ),
        reducer_xor_total_output_size=int(
            inverse_metrics["xor_total_output_size"]  # type: ignore[index]
        ),
        reducer_xor_max_input_size=int(
            inverse_metrics["xor_max_input_size"]  # type: ignore[index]
        ),
        reducer_xor_max_output_size=int(
            inverse_metrics["xor_max_output_size"]  # type: ignore[index]
        ),
        reducer_xor_inserted_labels=int(
            inverse_metrics["xor_inserted_labels"]  # type: ignore[index]
        ),
        reducer_xor_removed_labels=int(
            inverse_metrics["xor_removed_labels"]  # type: ignore[index]
        ),
        reducer_inverse_list_appends=int(
            inverse_metrics["inverse_list_appends"]  # type: ignore[index]
        ),
        sequence_seconds=sequence_seconds,
        reference_seconds=reference_seconds,
        morse_reduction_seconds=morse_reduction_seconds,
        reducer_setup_seconds=1.0e-9 * int(reducer_metrics["reducer_setup_nanoseconds"]),
        reducer_compute_seconds=1.0e-9 * int(reducer_metrics["reducer_compute_nanoseconds"]),
        morse_seconds=morse_seconds,
        standard_seconds=standard_seconds,
        repeats=repeats,
        finite_interval_count=finite_interval_count,
        essential_interval_count=essential_interval_count,
        barcodes_materialized=barcodes_materialized,
        validation_mode=validation_mode,
        finite_barcode=finite_barcode,
        essential_barcode=essential_barcode,
    )


def _make_morse_reference_profile(payload: dict[str, object]) -> MorseReferenceProfile:
    metrics = payload["metrics"]  # type: ignore[assignment]
    frame_metrics = payload["frame_metrics"]  # type: ignore[assignment]
    inverse_metrics = metrics["inverse_store"]  # type: ignore[index]
    profile_nanoseconds = int(payload.get("profile_nanoseconds", 0))
    profile_seconds = (
        float(payload["profile_seconds"])
        if "profile_seconds" in payload
        else 1.0e-9 * profile_nanoseconds
    )
    return MorseReferenceProfile(
        num_simplices=int(payload["num_simplices"]),
        num_levels=int(payload["num_levels"]),
        num_critical_simplices=int(payload["num_critical_simplices"]),
        num_regular_pairs=int(payload["num_regular_pairs"]),
        sequence_algorithm=str(payload["sequence_algorithm"]),
        frame_mode=str(payload["frame_mode"]),
        reference_final_live_nonempty_annotations=int(
            frame_metrics["final_live_nonempty_annotations"]  # type: ignore[index]
        ),
        reference_final_live_total_annotation_size=int(
            frame_metrics["final_live_total_annotation_size"]  # type: ignore[index]
        ),
        reference_peak_live_nonempty_annotations=int(
            frame_metrics["peak_live_nonempty_annotations"]  # type: ignore[index]
        ),
        reference_peak_live_total_annotation_size=int(
            frame_metrics["peak_live_total_annotation_size"]  # type: ignore[index]
        ),
        reference_released_annotations=int(frame_metrics["released_annotations"]),  # type: ignore[index]
        reference_released_total_annotation_size=int(
            frame_metrics["released_total_annotation_size"]  # type: ignore[index]
        ),
        reducer_working_set_size=int(metrics["working_set_size"]),  # type: ignore[index]
        reducer_initial_nonempty_annotations=int(
            inverse_metrics["initial_nonempty_annotations"]  # type: ignore[index]
        ),
        reducer_initial_total_annotation_size=int(
            inverse_metrics["initial_total_annotation_size"]  # type: ignore[index]
        ),
        reducer_initial_max_annotation_size=int(
            inverse_metrics["initial_max_annotation_size"]  # type: ignore[index]
        ),
        reducer_initial_inverse_list_entries=int(
            inverse_metrics["initial_inverse_list_entries"]  # type: ignore[index]
        ),
        reducer_boundary_plan_face_scans=int(metrics["boundary_plan_face_scans"]),  # type: ignore[index]
        reducer_boundary_annotation_candidate_criticals=int(
            metrics["boundary_annotation_candidate_criticals"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_zero_skipped_criticals=int(
            metrics["boundary_annotation_zero_skipped_criticals"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_zero_skipped_faces=int(
            metrics["boundary_annotation_zero_skipped_faces"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_xors=int(metrics["boundary_annotation_xors"]),  # type: ignore[index]
        reducer_boundary_annotation_total_input_size=int(
            metrics["boundary_annotation_total_input_size"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_total_output_size=int(
            metrics["boundary_annotation_total_output_size"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_max_size=int(
            metrics["boundary_annotation_max_size"]  # type: ignore[index]
        ),
        reducer_boundary_annotation_max_output_size=int(
            metrics["boundary_annotation_max_output_size"]  # type: ignore[index]
        ),
        estimated_reducer_work=int(payload["estimated_reducer_work"]),
        profile_seconds=profile_seconds,
        metrics=metrics,  # type: ignore[arg-type]
    )


def _profile_selection_key(
    profile: MorseReferenceProfile,
    selection_metric: str,
) -> tuple[float, float, float]:
    metric = selection_metric.lower().replace("_", "-")
    critical_ratio = profile.critical_ratio
    estimated_work = float(profile.estimated_reducer_work)
    boundary_work = float(profile.boundary_annotation_work)
    if metric in {
        "estimated-reducer-work",
        "profile-work",
        "reducer-work",
        "annotation-work",
        "work",
    }:
        return (estimated_work, critical_ratio, profile.profile_seconds)
    if metric in {
        "profile-total-work",
        "total-profile-work",
        "time-adjusted-work",
        "hybrid-work",
    }:
        total_work = estimated_work + 1.0e9 * profile.profile_seconds
        return (total_work, estimated_work, critical_ratio)
    if metric in {"critical-ratio", "critical-count", "critical"}:
        return (critical_ratio, estimated_work, profile.profile_seconds)
    if metric in {"profile-seconds", "time", "frame-seconds"}:
        return (profile.profile_seconds, estimated_work, critical_ratio)
    if metric in {"working-set-size", "working-set"}:
        return (float(profile.reducer_working_set_size), estimated_work, critical_ratio)
    if metric in {"initial-annotation-size", "initial-labels", "initial-volume"}:
        return (float(profile.reducer_initial_total_annotation_size), estimated_work, critical_ratio)
    if metric in {"boundary-annotation-work", "boundary-work"}:
        return (boundary_work, estimated_work, critical_ratio)
    supported = ", ".join(MORSE_PROFILE_SELECTION_METRICS)
    raise ValueError(f"Unknown profile selection metric {selection_metric!r}. Supported: {supported}.")


def _profile_morse_reference_frame_python(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    references: Sequence[Annotation],
    *,
    algorithm: str,
    profile_seconds: float,
) -> dict[str, object]:
    present = [False] * complex_.size
    boundary_candidates: list[tuple[int, int]] = []
    zero_boundary_critical_ids: list[int] = []
    boundary_plan_face_scans = 0
    zero_boundary_skipped_faces = 0

    for critical_id, sigma in enumerate(sequence.critical_simplices):
        present[sigma] = True
        boundary = complex_.boundary(sigma)
        boundary_plan_face_scans += len(boundary)
        may_have_nonzero_boundary = False
        for face in boundary:
            present[face] = True
            if references[face]:
                may_have_nonzero_boundary = True
        if may_have_nonzero_boundary:
            boundary_candidates.append((critical_id, sigma))
        else:
            zero_boundary_critical_ids.append(critical_id)
            zero_boundary_skipped_faces += len(boundary)

    working_set = tuple(simplex for simplex, is_present in enumerate(present) if is_present)
    working_set_index = {simplex: index for index, simplex in enumerate(working_set)}
    annotations = tuple(references[simplex] for simplex in working_set)

    inverse_metrics = _empty_reducer_metrics()["inverse_store"]  # type: ignore[index,assignment]
    for annotation in annotations:
        if not annotation:
            continue
        inverse_metrics["initial_nonempty_annotations"] += 1  # type: ignore[index,operator]
        inverse_metrics["initial_total_annotation_size"] += len(annotation)  # type: ignore[index,operator]
        inverse_metrics["initial_inverse_list_entries"] += len(annotation)  # type: ignore[index,operator]
        inverse_metrics["initial_max_annotation_size"] = max(  # type: ignore[index]
            inverse_metrics["initial_max_annotation_size"],  # type: ignore[index]
            len(annotation),
        )

    boundary_annotation_xors = 0
    boundary_annotation_total_input_size = 0
    boundary_annotation_total_output_size = 0
    boundary_annotation_max_size = 0
    boundary_annotation_max_output_size = 0
    for _, sigma in boundary_candidates:
        boundary_annotation: Annotation = ()
        for face in complex_.boundary(sigma):
            face_annotation = annotations[working_set_index[face]]
            boundary_annotation_xors += 1
            boundary_annotation_total_input_size += len(face_annotation)
            boundary_annotation = _xor_sorted(boundary_annotation, face_annotation)
        boundary_annotation_max_size = max(boundary_annotation_max_size, len(boundary_annotation))
        boundary_annotation_total_output_size += len(boundary_annotation)
        boundary_annotation_max_output_size = max(
            boundary_annotation_max_output_size,
            len(boundary_annotation),
        )

    full_reference_nonempty = sum(1 for annotation in references if annotation)
    full_reference_total = sum(len(annotation) for annotation in references)
    frame_metrics = {
        "remaining_cofaces_nanoseconds": 0,
        "sequence_total_nanoseconds": 0,
        "sequence_core_nanoseconds": 0,
        "reference_update_nanoseconds": 0,
        "reduction_plan_nanoseconds": 0,
        "release_cleanup_nanoseconds": 0,
        "working_set_pack_nanoseconds": 0,
        "local_index_nanoseconds": 0,
        "sequence_init_nanoseconds": 0,
        "sequence_candidate_seed_nanoseconds": 0,
        "sequence_candidate_loop_nanoseconds": 0,
        "sequence_emit_nanoseconds": 0,
        "sequence_callback_nanoseconds": 0,
        "sequence_replay_nanoseconds": 0,
        "sequence_candidate_pushes": 0,
        "sequence_candidate_pops": 0,
        "sequence_stale_candidate_skips": 0,
        "sequence_level_mismatch_skips": 0,
        "sequence_regular_pairs": 0,
        "sequence_criticals": 0,
        "final_live_nonempty_annotations": full_reference_nonempty,
        "final_live_total_annotation_size": full_reference_total,
        "peak_live_nonempty_annotations": full_reference_nonempty,
        "peak_live_total_annotation_size": full_reference_total,
        "released_annotations": 0,
        "released_total_annotation_size": 0,
    }
    metrics = _empty_reducer_metrics()
    metrics["working_set_size"] = len(working_set)
    metrics["critical_count"] = len(sequence.critical_simplices)
    metrics["boundary_plan_face_scans"] = boundary_plan_face_scans
    metrics["boundary_annotation_candidate_criticals"] = len(boundary_candidates)
    metrics["boundary_annotation_zero_skipped_criticals"] = len(zero_boundary_critical_ids)
    metrics["boundary_annotation_zero_skipped_faces"] = zero_boundary_skipped_faces
    metrics["boundary_annotation_xors"] = boundary_annotation_xors
    metrics["boundary_annotation_total_input_size"] = boundary_annotation_total_input_size
    metrics["boundary_annotation_total_output_size"] = boundary_annotation_total_output_size
    metrics["boundary_annotation_max_size"] = boundary_annotation_max_size
    metrics["boundary_annotation_max_output_size"] = boundary_annotation_max_output_size
    metrics["inverse_store"] = inverse_metrics
    estimated_work = (
        boundary_plan_face_scans
        + int(inverse_metrics["initial_total_annotation_size"])  # type: ignore[index]
        + boundary_annotation_xors
        + boundary_annotation_total_input_size
        + boundary_annotation_total_output_size
    )
    return {
        "num_simplices": complex_.size,
        "num_levels": complex_.num_levels,
        "num_critical_simplices": len(sequence.critical_simplices),
        "num_regular_pairs": sum(1 for step in sequence.steps if step.type == REGULAR_PAIR),
        "sequence_algorithm": algorithm,
        "frame_mode": FUSED_FRAME,
        "frame_metrics": frame_metrics,
        "metrics": metrics,
        "estimated_reducer_work": estimated_work,
        "profile_seconds": profile_seconds,
    }


def _require_cpp_backend() -> None:
    if _morse_core is None:
        raise RuntimeError("The C++ backend is not available. Build `_morse_core` first.")


def _empty_reducer_metrics() -> dict[str, object]:
    return {
        "working_set_size": 0,
        "critical_count": 0,
        "reducer_setup_nanoseconds": 0,
        "reducer_compute_nanoseconds": 0,
        "boundary_plan_face_scans": 0,
        "boundary_annotation_candidate_criticals": 0,
        "boundary_annotation_zero_skipped_criticals": 0,
        "boundary_annotation_zero_skipped_faces": 0,
        "boundary_annotation_xors": 0,
        "boundary_annotation_total_input_size": 0,
        "boundary_annotation_total_output_size": 0,
        "boundary_annotation_max_size": 0,
        "boundary_annotation_max_output_size": 0,
        "pivot_eliminations": 0,
        "finite_pairs": 0,
        "essential_intervals": 0,
        "inverse_store": {
            "initial_nonempty_annotations": 0,
            "initial_total_annotation_size": 0,
            "initial_max_annotation_size": 0,
            "initial_inverse_list_entries": 0,
            "remove_candidate_scans": 0,
            "remove_applied": 0,
            "remove_total_annotation_size": 0,
            "remove_max_annotation_size": 0,
            "xor_candidate_scans": 0,
            "xor_applied": 0,
            "xor_changed_labels": 0,
            "xor_total_input_size": 0,
            "xor_total_output_size": 0,
            "xor_max_input_size": 0,
            "xor_max_output_size": 0,
            "xor_inserted_labels": 0,
            "xor_removed_labels": 0,
            "inverse_list_appends": 0,
        },
    }


def _empty_frame_metrics() -> dict[str, object]:
    return {
        "remaining_cofaces_nanoseconds": 0,
        "sequence_total_nanoseconds": 0,
        "sequence_core_nanoseconds": 0,
        "reference_update_nanoseconds": 0,
        "reduction_plan_nanoseconds": 0,
        "release_cleanup_nanoseconds": 0,
        "working_set_pack_nanoseconds": 0,
        "local_index_nanoseconds": 0,
        "sequence_init_nanoseconds": 0,
        "sequence_candidate_seed_nanoseconds": 0,
        "sequence_candidate_loop_nanoseconds": 0,
        "sequence_emit_nanoseconds": 0,
        "sequence_callback_nanoseconds": 0,
        "sequence_replay_nanoseconds": 0,
        "sequence_candidate_pushes": 0,
        "sequence_candidate_pops": 0,
        "sequence_stale_candidate_skips": 0,
        "sequence_level_mismatch_skips": 0,
        "sequence_regular_pairs": 0,
        "sequence_criticals": 0,
        "final_live_nonempty_annotations": 0,
        "final_live_total_annotation_size": 0,
        "peak_live_nonempty_annotations": 0,
        "peak_live_total_annotation_size": 0,
        "released_annotations": 0,
        "released_total_annotation_size": 0,
    }


def _compute_saturated_morse_reference_frame_python(
    complex_: FilteredComplex,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> MorseReferenceFrame:
    n = complex_.size
    steps: list[MorseStep] = []
    critical_simplices: list[int] = []
    critical_index_of_simplex = [-1] * n
    paired_with: list[int | None] = [None] * n
    references: list[Annotation] = [() for _ in range(n)]
    inserted = [False] * n
    missing_count = [len(complex_.boundary(simplex)) for simplex in range(n)]
    unique_missing = [None] * n
    remaining_by_level = [0] * complex_.num_levels

    for simplex in range(n):
        if missing_count[simplex] == 1:
            unique_missing[simplex] = complex_.boundary(simplex)[0]
        remaining_by_level[complex_.level(simplex)] += 1

    def refresh_unique_missing(simplex: int) -> None:
        unique_missing[simplex] = None
        if missing_count[simplex] != 1:
            return
        for face in complex_.boundary(simplex):
            if not inserted[face]:
                unique_missing[simplex] = face
                return
        raise RuntimeError("Missing-face count is inconsistent.")

    def insert_simplex(simplex: int) -> None:
        if inserted[simplex]:
            raise RuntimeError("Tried to insert a simplex twice.")
        inserted[simplex] = True
        remaining_by_level[complex_.level(simplex)] -= 1
        for coface in complex_.coboundary(simplex):
            if inserted[coface]:
                continue
            if missing_count[coface] == 0:
                raise RuntimeError("Missing-face count underflow.")
            missing_count[coface] -= 1
            refresh_unique_missing(coface)

    for level in range(complex_.num_levels):
        bucket = complex_.simplices_of_level(level)
        while remaining_by_level[level] > 0:
            inserted_pair = False
            for tau in bucket:
                if inserted[tau] or complex_.dimension(tau) == 0 or missing_count[tau] != 1:
                    continue
                sigma = unique_missing[tau]
                if sigma is not None and not inserted[sigma] and complex_.level(sigma) == level:
                    steps.append(MorseStep(REGULAR_PAIR, sigma, tau, level))
                    paired_with[sigma] = tau
                    paired_with[tau] = sigma
                    references[tau] = ()
                    lower_reference: Annotation = ()
                    for face in complex_.boundary(tau):
                        if face != sigma:
                            lower_reference = _xor_sorted(lower_reference, references[face])
                    references[sigma] = lower_reference
                    insert_simplex(sigma)
                    insert_simplex(tau)
                    inserted_pair = True
                    break

            if inserted_pair:
                continue

            fillable = None
            for simplex in bucket:
                if not inserted[simplex] and missing_count[simplex] == 0:
                    fillable = simplex
                    break

            if fillable is None:
                raise RuntimeError("No valid F-sequence step found.")

            critical_id = len(critical_simplices)
            steps.append(MorseStep(CRITICAL, fillable, None, level))
            critical_simplices.append(fillable)
            critical_index_of_simplex[fillable] = critical_id
            references[fillable] = (critical_id,)
            insert_simplex(fillable)

    return MorseReferenceFrame(
        sequence=MorseSequence(
            steps=tuple(steps),
            critical_simplices=tuple(critical_simplices),
            critical_index_of_simplex=tuple(critical_index_of_simplex),
            paired_with=tuple(paired_with),
            algorithm=algorithm,
        ),
        _references=tuple(references),
    )


def _normalize_morse_sequence_algorithm(algorithm: str) -> str:
    normalized = str(algorithm).lower().replace("_", "-")
    aliases = {
        "f-sequence": SATURATED_SEQUENCE,
        "saturated-f-sequence": SATURATED_SEQUENCE,
        "paper-max": F_MAX_SEQUENCE,
        "max-s-f": F_MAX_SEQUENCE,
        "max-sf": F_MAX_SEQUENCE,
        "maximal-f-sequence": F_MAX_SEQUENCE,
        "paper-min": F_MIN_SEQUENCE,
        "min-s-f": F_MIN_SEQUENCE,
        "min-sf": F_MIN_SEQUENCE,
        "minimal-f-sequence": F_MIN_SEQUENCE,
        "plateau": PLATEAU_GREEDY_SEQUENCE,
        "plateau-greedy-f-sequence": PLATEAU_GREEDY_SEQUENCE,
        "coreduction": SAME_LEVEL_REDUCTION_SEQUENCE,
        "same-level-coreduction": SAME_LEVEL_REDUCTION_SEQUENCE,
        "coreduction-f-sequence": SAME_LEVEL_REDUCTION_SEQUENCE,
        "reduction-reverse": SAME_LEVEL_REDUCTION_SEQUENCE,
        "same-level-reduction-reverse": SAME_LEVEL_REDUCTION_SEQUENCE,
        "flooding": FLOODING_MINMAX_SEQUENCE,
        "flooding-maximal": FLOODING_MAX_SEQUENCE,
        "maximal-flooding": FLOODING_MAX_SEQUENCE,
        "flooding-minimal": FLOODING_MIN_SEQUENCE,
        "minimal-flooding": FLOODING_MIN_SEQUENCE,
        "flooding-minmax": FLOODING_MINMAX_SEQUENCE,
        "flooding-min-max": FLOODING_MINMAX_SEQUENCE,
        "min-max": FLOODING_MINMAX_SEQUENCE,
        "minmax": FLOODING_MINMAX_SEQUENCE,
        "flooding-maxmin": FLOODING_MAXMIN_SEQUENCE,
        "flooding-max-min": FLOODING_MAXMIN_SEQUENCE,
        "max-min": FLOODING_MAXMIN_SEQUENCE,
        "maxmin": FLOODING_MAXMIN_SEQUENCE,
    }
    normalized = aliases.get(normalized, normalized)
    if normalized in MORSE_SEQUENCE_ALGORITHMS:
        return normalized
    if normalized in RESERVED_MORSE_SEQUENCE_ALGORITHMS:
        raise NotImplementedError(
            f"Morse sequence algorithm {algorithm!r} is reserved for a future implementation."
        )
    supported = ", ".join(MORSE_SEQUENCE_ALGORITHMS)
    reserved = ", ".join(RESERVED_MORSE_SEQUENCE_ALGORITHMS)
    raise ValueError(
        f"Unknown Morse sequence algorithm {algorithm!r}. "
        f"Implemented: {supported}. Reserved for future work: {reserved}."
    )


def _normalize_morse_algorithm_portfolio(
    algorithms: Iterable[str] | str | None,
) -> tuple[str, ...]:
    if algorithms is None:
        values: Iterable[str] = DEFAULT_MORSE_ALGORITHM_PORTFOLIO
    elif isinstance(algorithms, str):
        normalized = algorithms.lower().replace("_", "-")
        if normalized in {"all", "portfolio", AUTO_MORSE_SEQUENCE_ALGORITHM}:
            values = DEFAULT_MORSE_ALGORITHM_PORTFOLIO
        else:
            values = (algorithms,)
    else:
        values = algorithms

    normalized_values: list[str] = []
    seen: set[str] = set()
    for algorithm in values:
        normalized_algorithm = _normalize_morse_sequence_algorithm(algorithm)
        if normalized_algorithm in seen:
            continue
        normalized_values.append(normalized_algorithm)
        seen.add(normalized_algorithm)

    if not normalized_values:
        raise ValueError("At least one Morse sequence algorithm is required.")
    return tuple(normalized_values)


def _normalize_morse_frame_mode(frame_mode: str) -> str:
    normalized = str(frame_mode).lower().replace("_", "-")
    if normalized in MORSE_FRAME_MODES:
        return normalized
    supported = ", ".join(MORSE_FRAME_MODES)
    raise ValueError(f"Unknown Morse frame mode {frame_mode!r}. Supported modes: {supported}.")


def _normalize_morse_selection_mode(selection_mode: str) -> str:
    normalized = str(selection_mode).lower().replace("_", "-")
    if normalized in MORSE_ALGORITHM_SELECTION_MODES:
        return normalized
    supported = ", ".join(MORSE_ALGORITHM_SELECTION_MODES)
    raise ValueError(f"Unknown Morse selection mode {selection_mode!r}. Supported modes: {supported}.")


def _core_compute_morse_sequence(cpp_complex: object, algorithm: str) -> object:
    try:
        return _morse_core.compute_morse_sequence(cpp_complex, algorithm)
    except TypeError:
        if algorithm != SATURATED_SEQUENCE:
            raise
        return _morse_core.compute_morse_sequence(cpp_complex)


def _compute_cpp_sequence(
    complex_: FilteredComplex,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> object | None:
    if _morse_core is None or complex_._cpp is None:
        return None
    return _core_compute_morse_sequence(complex_._cpp, algorithm)


def _compute_cpp_reference_frame(
    complex_: FilteredComplex,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> object | None:
    if _morse_core is None or complex_._cpp is None:
        return None
    builder = getattr(_morse_core, "compute_morse_sequence_and_reference_map_object", None)
    if builder is None:
        return None
    return builder(complex_._cpp, algorithm)


def _compute_cpp_coreference_frame(
    complex_: FilteredComplex,
    algorithm: str = COREDUCTION_SEQUENCE,
) -> object | None:
    if _morse_core is None or complex_._cpp is None:
        return None
    builder = getattr(_morse_core, "compute_morse_sequence_and_coreference_map_object", None)
    if builder is None:
        return None
    return builder(complex_._cpp, algorithm)


def _cpp_sequence_for(
    complex_: FilteredComplex,
    sequence: MorseSequence | None,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> object | None:
    if _morse_core is None or complex_._cpp is None:
        return None
    if sequence is None:
        return _compute_cpp_sequence(
            complex_,
            _normalize_morse_sequence_algorithm(algorithm),
        )
    if (
        sequence._cpp_sequence is not None
        and sequence._complex_signature == _complex_signature(complex_)
    ):
        return sequence._cpp_sequence
    return None


def _sequence_from_cpp_sequence(
    complex_: FilteredComplex,
    cpp_sequence: object,
    *,
    algorithm: str = DEFAULT_MORSE_SEQUENCE_ALGORITHM,
) -> MorseSequence:
    signature = _complex_signature(complex_)
    return MorseSequence(
        steps=tuple(
            MorseStep(
                type=str(kind),
                sigma=int(sigma),
                tau=None if tau is None else int(tau),
                level=int(level),
            )
            for kind, sigma, tau, level in cpp_sequence.steps()
        ),
        critical_simplices=tuple(int(simplex) for simplex in cpp_sequence.critical_simplices),
        critical_index_of_simplex=tuple(
            int(index) for index in cpp_sequence.critical_index_of_simplex
        ),
        paired_with=tuple(
            None if simplex is None else int(simplex)
            for simplex in cpp_sequence.paired_with(complex_._cpp)
        ),
        algorithm=algorithm,
        _cpp_sequence=cpp_sequence,
        _complex_signature=signature,
    )


def _references_from_cpp_annotations(annotations: Sequence[Sequence[int]]) -> tuple[Annotation, ...]:
    return tuple(tuple(int(label) for label in annotation) for annotation in annotations)


def _field_references_from_cpp_annotations(
    annotations: Sequence[Sequence[Sequence[int]]],
) -> tuple[FieldAnnotation, ...]:
    return tuple(
        tuple((int(label), int(coefficient)) for label, coefficient in annotation)
        for annotation in annotations
    )


def _diagram_from_cpp_dict(diagram: dict[str, object]) -> PersistenceDiagram:
    return PersistenceDiagram(
        finite_pairs=tuple(
            PersistencePair(
                birth=int(birth),
                death=int(death),
                dimension=int(dimension),
                birth_value=float(birth_value),
                death_value=float(death_value),
            )
            for birth, death, dimension, birth_value, death_value in diagram["finite_pairs"]  # type: ignore[index]
        ),
        essential=tuple(
            EssentialInterval(
                birth=int(birth),
                dimension=int(dimension),
                birth_value=float(birth_value),
            )
            for birth, dimension, birth_value in diagram["essential"]  # type: ignore[index]
        ),
    )


def _compute_cpp_payload(complex_: FilteredComplex) -> dict[str, object] | None:
    if _morse_core is None:
        return None
    if not getattr(complex_, "_simplicial_payload_supported", True):
        return None
    simplices = [
        (list(record.vertices), record.filtration)
        for record in complex_.simplices()
    ]
    return _morse_core.analyze(simplices)


def _payload_from_sequence_or_compute(
    complex_: FilteredComplex, sequence: MorseSequence | None
) -> dict[str, object] | None:
    signature = _complex_signature(complex_)
    if sequence is not None:
        if sequence._cpp_payload is not None and sequence._complex_signature == signature:
            return sequence._cpp_payload
        return None
    return _compute_cpp_payload(complex_)


def _complex_signature(complex_: FilteredComplex) -> tuple[tuple[Simplex, float], ...]:
    return tuple((record.vertices, record.filtration) for record in complex_.simplices())


def _sequence_from_cpp_payload(
    payload: dict[str, object], signature: tuple[tuple[Simplex, float], ...]
) -> MorseSequence:
    return MorseSequence(
        steps=tuple(
            MorseStep(
                type=str(kind),
                sigma=int(sigma),
                tau=None if tau is None else int(tau),
                level=int(level),
            )
            for kind, sigma, tau, level in payload["steps"]  # type: ignore[index]
        ),
        critical_simplices=tuple(int(simplex) for simplex in payload["critical_simplices"]),  # type: ignore[index]
        critical_index_of_simplex=tuple(
            int(index) for index in payload["critical_index_of_simplex"]  # type: ignore[index]
        ),
        paired_with=tuple(
            None if simplex is None else int(simplex)
            for simplex in payload["paired_with"]  # type: ignore[index]
        ),
        algorithm=DEFAULT_MORSE_SEQUENCE_ALGORITHM,
        _cpp_payload=payload,
        _complex_signature=signature,
    )


def _references_from_cpp_payload(payload: dict[str, object]) -> tuple[Annotation, ...]:
    return tuple(
        tuple(int(label) for label in annotation)
        for annotation in payload["references"]  # type: ignore[index]
    )


def _diagram_from_cpp_payload(payload: dict[str, object], key: str) -> PersistenceDiagram:
    return _diagram_from_cpp_dict(payload[key])  # type: ignore[index]


def _canonical_simplex(vertices: Sequence[int]) -> Simplex:
    simplex = tuple(sorted(int(vertex) for vertex in vertices))
    if not simplex:
        raise ValueError("A simplex must contain at least one vertex.")
    if len(set(simplex)) != len(simplex):
        raise ValueError("A simplex cannot contain duplicate vertices.")
    return simplex


def _all_nonempty_faces(facet: Simplex) -> Iterator[Simplex]:
    for size in range(1, len(facet) + 1):
        for simplex in combinations(facet, size):
            yield simplex


def _simplex_filtration_value(spec: FiltrationSpec, simplex: Simplex) -> float:
    if callable(spec):
        return float(spec(simplex))
    if isinstance(spec, Mapping):
        try:
            return float(spec[simplex])
        except KeyError as exc:
            raise KeyError(f"Missing filtration value for simplex {simplex}.") from exc
    return float(spec)


def _vertex_filtration_value(spec: VertexFiltrationSpec, vertex: int) -> float:
    if callable(spec):
        return float(spec(vertex))
    if isinstance(spec, Mapping):
        try:
            return float(spec[vertex])
        except KeyError as exc:
            raise KeyError(f"Missing filtration value for vertex {vertex}.") from exc
    if isinstance(spec, Sequence) and not isinstance(spec, (str, bytes)):
        return float(spec[vertex])
    return float(spec)


def _edge_filtration_value(spec: EdgeFiltrationSpec, edge: tuple[int, int]) -> float:
    canonical = (min(edge), max(edge))
    if callable(spec):
        return float(spec(canonical))
    if isinstance(spec, Mapping):
        try:
            return float(spec[canonical])
        except KeyError as exc:
            raise KeyError(f"Missing filtration value for edge {canonical}.") from exc
    return float(spec)


def _codimension_one_faces(simplex: Simplex) -> Iterator[Simplex]:
    if len(simplex) == 1:
        return
    for removed in range(len(simplex)):
        yield simplex[:removed] + simplex[removed + 1 :]


def _validate_prime_modulus(modulus: int) -> int:
    if isinstance(modulus, bool):
        raise TypeError("Coefficient modulus must be an integer prime.")
    try:
        value = modulus.__index__()
    except AttributeError as exc:
        raise TypeError("Coefficient modulus must be an integer prime.") from exc

    if value < 2:
        raise ValueError("Coefficient modulus must be a prime integer at least 2.")
    if value == 2:
        return value
    if value % 2 == 0:
        raise ValueError("Coefficient modulus must be prime; composite moduli are not supported.")

    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            raise ValueError("Coefficient modulus must be prime; composite moduli are not supported.")
        divisor += 2
    return value


def _oriented_boundary_column(
    complex_: FilteredComplex,
    simplex_id: int,
    order_index: Mapping[int, int],
    modulus: int,
) -> dict[int, int]:
    column: dict[int, int] = {}
    for removed_index, face in enumerate(complex_.boundary(simplex_id)):
        coefficient = _boundary_coefficient_for(complex_, simplex_id, removed_index, modulus)
        row = order_index[face]
        column[row] = coefficient
    return column


def _boundary_coefficient(removed_index: int, modulus: int) -> int:
    return 1 if removed_index % 2 == 0 else modulus - 1


def _boundary_coefficient_for(
    complex_: FilteredComplex,
    simplex_id: int,
    boundary_index: int,
    modulus: int,
) -> int:
    coefficient_method = getattr(complex_, "boundary_coefficient", None)
    if coefficient_method is None:
        return _boundary_coefficient(boundary_index, modulus)
    return int(coefficient_method(simplex_id, boundary_index, modulus)) % modulus


def _boundary_coefficient_of_face(
    complex_: FilteredComplex,
    *,
    coface: int,
    face: int,
    modulus: int,
) -> int:
    for removed_index, candidate in enumerate(complex_.boundary(coface)):
        if candidate == face:
            return _boundary_coefficient_for(complex_, coface, removed_index, modulus)
    raise RuntimeError("Expected a codimension-one face/coface incidence.")


def _field_annotation_from_labels(annotation: Annotation) -> FieldAnnotation:
    return tuple((int(label), 1) for label in annotation)


def _field_annotation_from_dict(annotation: Mapping[int, int]) -> FieldAnnotation:
    return tuple(
        (int(label), int(coefficient))
        for label, coefficient in sorted(annotation.items())
        if coefficient
    )


def _field_annotation_to_dict(
    annotation: Mapping[int, int] | Sequence[tuple[int, int]] | Sequence[int],
    modulus: int,
) -> dict[int, int]:
    result: dict[int, int] = {}
    if isinstance(annotation, Mapping):
        items = annotation.items()
    else:
        items = []
        for item in annotation:
            if isinstance(item, (tuple, list)) and len(item) == 2:
                label, coefficient = item
            else:
                label, coefficient = item, 1
            items.append((label, coefficient))

    for label, coefficient in items:
        value = int(coefficient) % modulus
        if value:
            result[int(label)] = (result.get(int(label), 0) + value) % modulus
            if result[int(label)] == 0:
                result.pop(int(label), None)
    return result


def _add_oriented_boundary_references(
    target: dict[int, int],
    complex_: FilteredComplex,
    simplex_id: int,
    references: Sequence[Mapping[int, int]],
    modulus: int,
) -> None:
    for removed_index, face in enumerate(complex_.boundary(simplex_id)):
        _modp_add_scaled(
            target,
            references[face],
            _boundary_coefficient_for(complex_, simplex_id, removed_index, modulus),
            modulus,
        )


def _add_oriented_coboundary_coreferences(
    target: dict[int, int],
    complex_: FilteredComplex,
    simplex_id: int,
    coreferences: Sequence[Mapping[int, int]],
    modulus: int,
) -> None:
    for coface in complex_.coboundary(simplex_id):
        _modp_add_scaled(
            target,
            coreferences[coface],
            _boundary_coefficient_of_face(
                complex_,
                coface=coface,
                face=simplex_id,
                modulus=modulus,
            ),
            modulus,
        )


def _modp_low(column: Mapping[int, int]) -> int | None:
    if not column:
        return None
    return max(column)


def _modp_min(column: Mapping[int, int]) -> int | None:
    if not column:
        return None
    return min(column)


def _modp_inverse(value: int, modulus: int) -> int:
    return pow(value % modulus, modulus - 2, modulus)


def _modp_scale_in_place(column: dict[int, int], scale: int, modulus: int) -> None:
    scale %= modulus
    if scale == 0:
        column.clear()
        return
    for row in list(column):
        value = (column[row] * scale) % modulus
        if value:
            column[row] = value
        else:
            column.pop(row, None)


def _modp_add_scaled(
    target: dict[int, int],
    source: Mapping[int, int],
    scale: int,
    modulus: int,
) -> None:
    scale %= modulus
    if scale == 0:
        return
    for row, coefficient in source.items():
        value = (target.get(row, 0) + scale * coefficient) % modulus
        if value:
            target[row] = value
        else:
            target.pop(row, None)


def _xor_sorted(lhs: Annotation, rhs: Annotation) -> Annotation:
    result: list[int] = []
    i = j = 0
    while i < len(lhs) or j < len(rhs):
        if j == len(rhs) or (i < len(lhs) and lhs[i] < rhs[j]):
            result.append(lhs[i])
            i += 1
        elif i == len(lhs) or rhs[j] < lhs[i]:
            result.append(rhs[j])
            j += 1
        else:
            i += 1
            j += 1
    return tuple(result)


def _remove_label_from_all(annotations: list[Annotation], label: int) -> None:
    for simplex_id, annotation in enumerate(annotations):
        if label in annotation:
            annotations[simplex_id] = tuple(value for value in annotation if value != label)


def _xor_into_all_containing(annotations: list[Annotation], label: int, value: Annotation) -> None:
    for simplex_id, annotation in enumerate(annotations):
        if label in annotation:
            annotations[simplex_id] = _xor_sorted(annotation, value)


__all__ = [
    "__version__",
    "AUTO_MORSE_SEQUENCE_ALGORITHM",
    "CRITICAL",
    "COREDUCTION_SEQUENCE",
    "DEFAULT_MORSE_ALGORITHM_PORTFOLIO",
    "DEFAULT_MORSE_SEQUENCE_ALGORITHM",
    "DEFAULT_MORSE_FRAME_MODE",
    "FUSED_FRAME",
    "MORSE_ALGORITHM_SELECTION_METRICS",
    "MORSE_ALGORITHM_SELECTION_MODES",
    "MORSE_FRAME_MODES",
    "MORSE_PROFILE_SELECTION_METRICS",
    "MORSE_SEQUENCE_ALGORITHMS",
    "PLATEAU_GREEDY_SEQUENCE",
    "PROFILE_SELECTION_MODE",
    "REGULAR_PAIR",
    "RESERVED_MORSE_SEQUENCE_ALGORITHMS",
    "SATURATED_SEQUENCE",
    "SEPARATE_FRAME",
    "SAME_LEVEL_REDUCTION_SEQUENCE",
    "BENCHMARK_SELECTION_MODE",
    "AdaptivePersistenceResult",
    "Annotation",
    "CppCubicalGrid2DComplex",
    "CppFilteredComplex",
    "CppMorseCoreferenceFrame",
    "CppMorseReferenceFrame",
    "CppMorseSequence",
    "CppReferenceMap",
    "CppSimplexTreeBuilder",
    "CoreductionDirectionBenchmark",
    "CubicalCellRecord",
    "CubicalGrid2DComplex",
    "EdgeFiltrationSpec",
    "EssentialInterval",
    "EssentialIntervalSimplices",
    "FieldAnnotation",
    "FiltrationSpec",
    "FilteredComplex",
    "FilteredComplexBuilder",
    "F_MAX_SEQUENCE",
    "F_MIN_SEQUENCE",
    "FLOODING_MAXMIN_SEQUENCE",
    "FLOODING_MAX_SEQUENCE",
    "FLOODING_MINMAX_SEQUENCE",
    "FLOODING_MIN_SEQUENCE",
    "PersistenceBenchmark",
    "MorseComplex",
    "MorseCoreferenceFrame",
    "MorseReferenceFrame",
    "MorseReferenceProfile",
    "MorseSequence",
    "MorseStep",
    "MorseStepSimplices",
    "PersistenceDiagram",
    "PersistencePair",
    "PersistencePairSimplices",
    "Simplex",
    "SimplexTreeBuilder",
    "SameLevelReductionDirectionBenchmark",
    "SimplexRecord",
    "VertexFiltrationSpec",
    "assert_matches_gudhi_cubical",
    "assert_matches_gudhi",
    "assert_matches_standard",
    "benchmark_coreduction_directions",
    "benchmark_same_level_reduction_directions",
    "benchmark_persistence",
    "benchmark_sequence_algorithms",
    "annotation_as_simplices",
    "compute_coreference_complex",
    "compute_coreference_map_modp",
    "compute_morse_complex",
    "compute_morse_coreference_persistence",
    "compute_morse_coreference_persistence_modp",
    "compute_morse_persistence_modp",
    "compute_morse_persistence",
    "compute_morse_sequence_and_coreference_map",
    "compute_morse_sequence_and_reference_map",
    "compute_persistence_portfolio",
    "compute_persistence_profile_portfolio",
    "compute_persistence_adaptive",
    "compute_reference_complex",
    "compute_morse_sequence",
    "compute_coreference_map",
    "compute_reference_map_modp",
    "compute_reference_map",
    "compute_standard_persistence",
    "compute_standard_persistence_modp",
    "coreference_map_as_simplices",
    "cpp_backend_active",
    "cpp_backend_available",
    "cpp_compute_coreference_map",
    "cpp_compute_coreference_map_modp",
    "cpp_compute_coreference_map_object",
    "cpp_compute_morse_coreference_persistence",
    "cpp_compute_morse_coreference_persistence_modp",
    "cpp_compute_morse_persistence",
    "cpp_compute_morse_persistence_modp",
    "cpp_compute_morse_sequence",
    "cpp_compute_morse_sequence_and_coreference_map_object",
    "cpp_compute_morse_sequence_and_reference_map_object",
    "cpp_compute_reference_map",
    "cpp_compute_reference_map_modp",
    "cpp_compute_reference_map_object",
    "cpp_compute_standard_persistence",
    "cpp_compute_standard_persistence_modp",
    "cpp_cubical_grid_2d_from_vertex_values",
    "cpp_filtered_complex_from_simplices",
    "cpp_profile_morse_reference_frame",
    "cpp_reference_map_to_tuple",
    "cpp_reduce_morse_coreference_persistence",
    "cpp_reduce_morse_persistence",
    "cpp_reduce_morse_persistence_with_metrics",
    "gudhi_cubical_barcode",
    "gudhi_barcode",
    "morse_sequence_as_simplices",
    "profile_morse_reference_frame",
    "profile_morse_sequence_algorithms",
    "reference_map_as_simplices",
    "select_morse_sequence_algorithm",
    "select_morse_sequence_profile",
]
