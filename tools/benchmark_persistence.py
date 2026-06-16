#!/usr/bin/env python3
"""Run reproducible Morse-vs-standard persistence benchmark families."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path
from time import perf_counter
from typing import Iterable, TextIO
from urllib.error import URLError
from urllib.request import urlretrieve


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import morseframes as mp  # noqa: E402
from calibrate_profile_gate import shape_bin_from_values  # noqa: E402


BASE_FAMILIES = ("lower-star", "plateau", "rips")
TERRAIN_FAMILIES = ("terrain",)
IMAGE_FAMILIES = ("image-grid",)
CAM_FAMILIES = ("cam-s4-rips",)
ROADMAP_FAMILIES = ("roadmap-rips",)
FAMILIES = BASE_FAMILIES + TERRAIN_FAMILIES + IMAGE_FAMILIES + CAM_FAMILIES + ROADMAP_FAMILIES
DEFAULT_ROADMAP_CACHE = ROOT.parent / "work" / "roadmap-data"
PROFILE_CANDIDATE_GATES = ("all", "family-aware", "shape-aware", "tabulated")
BENCHMARK_PROFILE_SELECTION_METRICS = mp.MORSE_PROFILE_SELECTION_METRICS + (
    "adaptive",
    "family-adaptive",
    "adaptive_structured",
    "structured-adaptive",
    "grid-adaptive",
)


@dataclass(frozen=True)
class BenchmarkRow:
    family: str
    name: str
    seed: int
    size: int
    num_simplices: int
    num_levels: int
    num_critical_simplices: int
    sequence_algorithm: str
    frame_mode: str
    validation_mode: str
    barcodes_materialized: bool
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
    finite_intervals: int
    essential_intervals: int
    repeats: int
    sequence_seconds: float
    reference_seconds: float
    morse_reduction_seconds: float
    reducer_setup_seconds: float
    reducer_compute_seconds: float
    morse_seconds: float
    standard_seconds: float
    gudhi_cam_seconds: float | None
    gudhi_cam_speedup: float | None
    perseus_seconds: float | None
    perseus_speedup: float | None
    perseus_matches: bool | None
    speedup: float
    cpp_backend: bool


@dataclass(frozen=True)
class PerseusTiming:
    seconds: float
    finite_barcode: tuple[tuple[int, float, float], ...]
    essential_barcode: tuple[tuple[int, float], ...]
    matches_standard: bool


_EXTERNAL_TIMING_UNSET = object()


@dataclass(frozen=True)
class BenchmarkPreset:
    families: tuple[str, ...]
    sizes: tuple[int, ...]
    seeds: tuple[int, ...]
    repeats: int
    plateau_levels: int = 3


@dataclass(frozen=True)
class RoadmapDataset:
    label: str
    filename: str
    url: str
    max_dimension: int = 2


@dataclass(frozen=True)
class BenchmarkSummary:
    family: str
    sequence_algorithm: str
    frame_mode: str
    validation_mode: str
    size: int
    cases: int
    avg_num_simplices: float
    avg_num_levels: float
    avg_num_critical_simplices: float
    avg_critical_ratio: float
    avg_reference_final_live_nonempty_annotations: float
    avg_reference_final_live_total_annotation_size: float
    avg_reference_peak_live_nonempty_annotations: float
    avg_reference_peak_live_total_annotation_size: float
    avg_reference_released_annotations: float
    avg_reference_released_total_annotation_size: float
    avg_reducer_working_set_size: float
    avg_reducer_initial_nonempty_annotations: float
    avg_reducer_initial_total_annotation_size: float
    avg_reducer_initial_max_annotation_size: float
    avg_reducer_initial_inverse_list_entries: float
    avg_reducer_boundary_plan_face_scans: float
    avg_reducer_boundary_annotation_candidate_criticals: float
    avg_reducer_boundary_annotation_zero_skipped_criticals: float
    avg_reducer_boundary_annotation_zero_skipped_faces: float
    avg_reducer_boundary_annotation_xors: float
    avg_reducer_boundary_annotation_total_output_size: float
    avg_reducer_pivot_eliminations: float
    avg_reducer_remove_candidate_scans: float
    avg_reducer_xor_candidate_scans: float
    avg_reducer_xor_applied: float
    avg_reducer_xor_total_input_size: float
    avg_reducer_xor_total_output_size: float
    avg_reducer_xor_max_input_size: float
    avg_reducer_xor_max_output_size: float
    avg_reducer_xor_inserted_labels: float
    avg_reducer_xor_removed_labels: float
    avg_sequence_share: float
    avg_reference_share: float
    avg_reduction_share: float
    avg_reducer_setup_seconds: float
    avg_reducer_compute_seconds: float
    avg_morse_seconds: float
    avg_standard_seconds: float
    avg_gudhi_cam_seconds: float | None
    avg_gudhi_cam_speedup: float | None
    avg_perseus_seconds: float | None
    avg_perseus_speedup: float | None
    avg_speedup: float
    min_speedup: float
    max_speedup: float


@dataclass(frozen=True)
class ComplexShapeProfile:
    max_dimension: int
    num_vertices: int
    num_edges: int
    num_triangles: int
    num_higher_simplices: int
    simplex_vertex_ratio: float
    vertex_sqrt: float
    edge_density: float
    triangle_density: float
    avg_boundary_size: float
    max_boundary_size: int
    avg_coboundary_size: float
    max_coboundary_size: int
    largest_level_size: int
    avg_level_size: float
    largest_level_ratio: float
    singleton_level_ratio: float
    level_count_ratio: float
    level_concentration: float
    level_entropy_ratio: float


@dataclass(frozen=True)
class ProfileVsMeasuredRow:
    family: str
    name: str
    seed: int
    size: int
    num_simplices: int
    num_levels: int
    shape_max_dimension: int
    shape_num_vertices: int
    shape_num_edges: int
    shape_num_triangles: int
    shape_num_higher_simplices: int
    shape_simplex_vertex_ratio: float
    shape_vertex_sqrt: float
    shape_edge_density: float
    shape_triangle_density: float
    shape_avg_boundary_size: float
    shape_avg_coboundary_size: float
    shape_max_coboundary_size: int
    shape_largest_level_size: int
    shape_avg_level_size: float
    shape_largest_level_ratio: float
    shape_singleton_level_ratio: float
    shape_level_count_ratio: float
    shape_level_concentration: float
    shape_level_entropy_ratio: float
    candidate_count: int
    profile_candidate_count: int
    profile_candidate_gate: str
    frame_mode: str
    validation_mode: str
    profile_selection_metric: str
    effective_profile_selection_metric: str
    measured_selection_metric: str
    profile_selected_algorithm: str
    measured_best_algorithm: str
    profile_matches_measured: bool
    profile_selected_measured_rank: int
    measured_best_profile_rank: int
    measured_best_morse_seconds: float
    profile_selected_morse_seconds: float
    profile_penalty_seconds: float
    profile_penalty_ratio: float
    profile_penalty_percent: float
    measured_best_critical_ratio: float
    profile_selected_critical_ratio: float
    measured_best_profile_score: float | None
    profile_selected_profile_score: float
    measured_best_estimated_reducer_work: int | None
    profile_selected_estimated_reducer_work: int
    measured_best_profile_seconds: float | None
    profile_selected_profile_seconds: float
    measured_best_sequence_seconds: float
    profile_selected_sequence_seconds: float
    measured_best_reducer_compute_seconds: float
    profile_selected_reducer_compute_seconds: float
    measured_best_reducer_work: float
    profile_selected_reducer_work: float
    total_profile_seconds: float
    total_measured_morse_seconds: float
    standard_seconds: float
    gudhi_cam_seconds: float | None
    measured_best_gudhi_cam_speedup: float | None
    profile_selected_gudhi_cam_speedup: float | None
    perseus_seconds: float | None
    measured_best_perseus_speedup: float | None
    profile_selected_perseus_speedup: float | None
    perseus_matches: bool | None
    cpp_backend: bool


PRESETS = {
    "smoke": BenchmarkPreset(BASE_FAMILIES, (5,), (0,), 1),
    "grid": BenchmarkPreset(BASE_FAMILIES, (8, 12, 16, 20, 24), (0, 1, 2, 3, 4), 5),
    "scale": BenchmarkPreset(BASE_FAMILIES, (16, 24, 32, 40, 48), (0, 1, 2), 3),
    "terrain": BenchmarkPreset(TERRAIN_FAMILIES, (6, 8, 10, 12), (0, 1, 2), 3, 6),
    "image-grid": BenchmarkPreset(IMAGE_FAMILIES, (6, 8, 10, 12), (0, 1, 2), 3, 4),
    "cam": BenchmarkPreset(CAM_FAMILIES, (16, 24, 32), (0, 1, 2), 3),
    "roadmap": BenchmarkPreset(ROADMAP_FAMILIES, (50, 100, 104), (0,), 3),
}

DEFAULT_PRESET = BenchmarkPreset(BASE_FAMILIES, (6, 8, 10), (0, 1, 2), 3)

ROADMAP_DATASETS = {
    50: RoadmapDataset(
        label="random50-16d",
        filename="random_point_cloud_50_16_.txt_distmat.txt",
        url=(
            "https://raw.githubusercontent.com/n-otter/PH-roadmap/master/"
            "data_sets/roadmap_datasets_distmat/random_point_cloud_50_16_.txt_distmat.txt"
        ),
    ),
    100: RoadmapDataset(
        label="random100-4d",
        filename="random_point_cloud_100_4_.txt_distmat.txt",
        url=(
            "https://raw.githubusercontent.com/n-otter/PH-roadmap/master/"
            "data_sets/roadmap_datasets_distmat/random_point_cloud_100_4_.txt_distmat.txt"
        ),
    ),
    104: RoadmapDataset(
        label="senate104",
        filename="senate104_edge_list.txt_0.68902_distmat.txt",
        url=(
            "https://raw.githubusercontent.com/n-otter/PH-roadmap/master/"
            "data_sets/roadmap_datasets_distmat/senate104_edge_list.txt_0.68902_distmat.txt"
        ),
    ),
}


def make_lower_star(seed: int, n: int) -> mp.FilteredComplex:
    rng = random.Random(seed)
    facets = [[vertex] for vertex in range(n)]

    if n >= 2:
        for _ in range(max(2, n)):
            facet_size = rng.randint(2, min(4, n))
            facets.append(rng.sample(range(n), facet_size))

    vertex_order = list(range(n))
    rng.shuffle(vertex_order)
    vertex_values = {vertex: float(value) for value, vertex in enumerate(vertex_order)}
    return mp.FilteredComplex.from_lower_star(
        facets,
        vertex_values,
        dimension_offset=0.0,
    )


def facet_closure(facets: Iterable[Iterable[int]]) -> list[tuple[int, ...]]:
    simplices: set[tuple[int, ...]] = set()
    for facet in facets:
        vertices = tuple(sorted(set(facet)))
        for dimension in range(1, len(vertices) + 1):
            simplices.update(combinations(vertices, dimension))
    return sorted(simplices, key=lambda simplex: (len(simplex), simplex))


def make_plateau(seed: int, n: int, *, levels: int = 3) -> mp.FilteredComplex:
    if levels <= 0:
        raise ValueError("levels must be positive")

    rng = random.Random(seed)
    facets = [[vertex] for vertex in range(n)]

    if n >= 2:
        for _ in range(max(2, n)):
            facet_size = rng.randint(2, min(4, n))
            facets.append(rng.sample(range(n), facet_size))

    filtration: dict[tuple[int, ...], float] = {}
    for simplex in facet_closure(facets):
        if len(simplex) == 1:
            lower_bound = 0
        else:
            lower_bound = max(
                int(filtration[face])
                for face in combinations(simplex, len(simplex) - 1)
            )
        filtration[simplex] = float(max(lower_bound, rng.randrange(levels)))

    return mp.FilteredComplex.from_simplices(
        (simplex, value) for simplex, value in filtration.items()
    )


def make_terrain(seed: int, n: int, *, levels: int = 6) -> mp.FilteredComplex:
    """Triangulated 2D terrain with quantized vertex values.

    The size parameter is the grid side length. Quantization intentionally creates
    same-level plateaus while keeping the filtration tied to a spatial scalar field.
    """

    if n < 2:
        raise ValueError("terrain size must be at least 2")
    if levels <= 0:
        raise ValueError("levels must be positive")

    rng = random.Random(seed)
    bumps = [
        (
            rng.random(),
            rng.random(),
            rng.uniform(-0.8, 0.8),
            rng.uniform(0.08, 0.22),
        )
        for _ in range(5)
    ]

    def vertex_id(row: int, col: int) -> int:
        return row * n + col

    raw_values: dict[int, float] = {}
    for row in range(n):
        x = row / float(n - 1)
        for col in range(n):
            y = col / float(n - 1)
            value = (
                0.35 * x
                + 0.20 * y
                + 0.35 * math.sin(2.0 * math.pi * (x + 0.17 * seed))
                + 0.25 * math.cos(2.0 * math.pi * (1.7 * y - 0.11 * seed))
                + 0.15 * math.sin(2.0 * math.pi * (x + y))
            )
            for center_x, center_y, weight, sigma in bumps:
                squared_distance = (x - center_x) ** 2 + (y - center_y) ** 2
                value += weight * math.exp(-squared_distance / (2.0 * sigma * sigma))
            value += rng.uniform(-0.03, 0.03)
            raw_values[vertex_id(row, col)] = value

    min_value = min(raw_values.values())
    max_value = max(raw_values.values())
    span = max_value - min_value
    vertex_values: dict[int, float] = {}
    for vertex, value in raw_values.items():
        normalized = 0.0 if span == 0.0 else (value - min_value) / span
        quantized = int(round(normalized * float(levels - 1)))
        vertex_values[vertex] = float(quantized)

    facets: list[tuple[int, int, int]] = []
    for row in range(n - 1):
        for col in range(n - 1):
            v00 = vertex_id(row, col)
            v10 = vertex_id(row + 1, col)
            v01 = vertex_id(row, col + 1)
            v11 = vertex_id(row + 1, col + 1)
            if (row + col + seed) % 2 == 0:
                facets.append((v00, v10, v11))
                facets.append((v00, v11, v01))
            else:
                facets.append((v00, v10, v01))
                facets.append((v10, v11, v01))

    return mp.FilteredComplex.from_lower_star(
        facets,
        vertex_values,
        dimension_offset=0.0,
    )


def make_image_grid(seed: int, n: int, *, levels: int = 4) -> mp.FilteredComplex:
    """Image-like triangulated grid with deliberately strong plateaus."""

    if n < 2:
        raise ValueError("image-grid size must be at least 2")
    if levels <= 0:
        raise ValueError("levels must be positive")

    rng = random.Random(seed)
    coarse_side = max(2, min(5, n // 2))
    coarse_values = [
        [rng.randrange(levels) for _col in range(coarse_side)]
        for _row in range(coarse_side)
    ]
    disk_center = (rng.uniform(0.25, 0.75), rng.uniform(0.25, 0.75))
    disk_radius = rng.uniform(0.16, 0.30)
    stripe_offset = rng.uniform(-0.25, 0.25)
    stripe_width = rng.uniform(0.08, 0.16)
    rect_x0 = rng.uniform(0.05, 0.35)
    rect_y0 = rng.uniform(0.05, 0.35)
    rect_x1 = min(0.95, rect_x0 + rng.uniform(0.25, 0.45))
    rect_y1 = min(0.95, rect_y0 + rng.uniform(0.25, 0.45))

    def vertex_id(row: int, col: int) -> int:
        return row * n + col

    vertex_values: dict[int, float] = {}
    for row in range(n):
        x = row / float(n - 1)
        coarse_row = min(coarse_side - 1, int(x * coarse_side))
        for col in range(n):
            y = col / float(n - 1)
            coarse_col = min(coarse_side - 1, int(y * coarse_side))
            value = coarse_values[coarse_row][coarse_col]

            if rect_x0 <= x <= rect_x1 and rect_y0 <= y <= rect_y1:
                value = 0
            squared_distance = (x - disk_center[0]) ** 2 + (y - disk_center[1]) ** 2
            if squared_distance <= disk_radius * disk_radius:
                value = levels - 1
            if abs((x - y) - stripe_offset) <= stripe_width:
                value = min(levels - 1, value + 1)

            vertex_values[vertex_id(row, col)] = float(value)

    facets: list[tuple[int, int, int]] = []
    for row in range(n - 1):
        for col in range(n - 1):
            v00 = vertex_id(row, col)
            v10 = vertex_id(row + 1, col)
            v01 = vertex_id(row, col + 1)
            v11 = vertex_id(row + 1, col + 1)
            if (row + col + seed) % 2 == 0:
                facets.append((v00, v10, v11))
                facets.append((v00, v11, v01))
            else:
                facets.append((v00, v10, v01))
                facets.append((v10, v11, v01))

    return mp.FilteredComplex.from_lower_star(
        facets,
        vertex_values,
        dimension_offset=0.0,
    )


def make_rips(seed: int, n: int) -> mp.FilteredComplex:
    rng = random.Random(seed)
    distances = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            value = float(rng.randint(1, 6))
            distances[i][j] = value
            distances[j][i] = value

    return mp.FilteredComplex.from_rips_distance_matrix(
        distances,
        max_edge_length=3.0,
        max_dimension=2,
    )


def make_cam_s4_rips(seed: int, n: int) -> mp.FilteredComplex:
    # CAM paper S4 profile: points on unit S^4 in R^5, Rips threshold 0.715,
    # and dimension cap 5. We downscale n for local prototype benchmarks.
    rng = random.Random(seed)
    points = []
    for _ in range(n):
        point = [rng.gauss(0.0, 1.0) for _ in range(5)]
        norm = math.sqrt(sum(coordinate * coordinate for coordinate in point))
        points.append([coordinate / norm for coordinate in point])

    distances = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            value = math.sqrt(
                sum((points[i][coordinate] - points[j][coordinate]) ** 2 for coordinate in range(5))
            )
            distances[i][j] = value
            distances[j][i] = value

    return mp.FilteredComplex.from_rips_distance_matrix(
        distances,
        max_edge_length=0.715,
        max_dimension=5,
    )


def read_distance_matrix(path: Path) -> list[list[float]]:
    distances = [
        [float(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip()
    ]
    if not distances:
        raise ValueError(f"Distance matrix file is empty: {path}")
    n = len(distances)
    if any(len(row) != n for row in distances):
        raise ValueError(f"Distance matrix must be square: {path}")
    return distances


def roadmap_dataset_path(
    size: int,
    *,
    cache_dir: Path = DEFAULT_ROADMAP_CACHE,
    download: bool = False,
) -> Path:
    try:
        dataset = ROADMAP_DATASETS[size]
    except KeyError as exc:
        supported = ", ".join(str(key) for key in sorted(ROADMAP_DATASETS))
        raise ValueError(f"Unsupported Roadmap dataset size {size}; supported sizes: {supported}") from exc

    path = cache_dir / dataset.filename
    if path.exists():
        return path
    if not download:
        raise FileNotFoundError(
            f"Missing Roadmap dataset {dataset.label} at {path}. "
            "Re-run with --download-roadmap-data to download it."
        )

    cache_dir.mkdir(parents=True, exist_ok=True)
    try:
        urlretrieve(dataset.url, path)
    except URLError as exc:
        raise RuntimeError(f"Could not download Roadmap dataset {dataset.label}.") from exc
    return path


def make_roadmap_rips(
    size: int,
    *,
    cache_dir: Path = DEFAULT_ROADMAP_CACHE,
    download: bool = False,
) -> mp.FilteredComplex:
    dataset = ROADMAP_DATASETS[size]
    distances = read_distance_matrix(
        roadmap_dataset_path(size, cache_dir=cache_dir, download=download)
    )
    return mp.FilteredComplex.from_rips_distance_matrix(
        distances,
        max_edge_length=math.inf,
        max_dimension=dataset.max_dimension,
    )


def time_gudhi_cam(complex_: mp.FilteredComplex, *, repeats: int = 3) -> float | None:
    try:
        import gudhi  # type: ignore[import-not-found]
    except ImportError:
        return None

    best = math.inf
    simplices = [
        (list(complex_.vertices(simplex_id)), complex_.filtration(simplex_id))
        for simplex_id in complex_.filtration_order
    ]
    for _ in range(repeats):
        simplex_tree = gudhi.SimplexTree()
        for simplex, filtration in simplices:
            simplex_tree.insert(simplex, filtration=filtration)

        started = perf_counter()
        simplex_tree.persistence(min_persistence=0, persistence_dim_max=True)
        best = min(best, perf_counter() - started)

    return best


def _resolve_perseus_executable(executable: str | Path) -> str | None:
    value = str(executable)
    path = Path(value)
    if path.is_file() and os.access(path, os.X_OK):
        return str(path)
    return shutil.which(value)


def write_perseus_nmfsimtop(
    complex_: mp.FilteredComplex,
    path: Path,
) -> dict[int, float]:
    """Write a Perseus nmfsimtop file and return integer birth -> filtration value."""

    level_to_birth = {
        value: index + 1
        for index, value in enumerate(complex_.level_values)
    }
    birth_to_level = {
        birth: value
        for value, birth in level_to_birth.items()
    }
    vertices = sorted(
        {
            vertex
            for simplex_id in complex_.filtration_order
            for vertex in complex_.vertices(simplex_id)
        }
    )
    vertex_coordinates = {vertex: index + 1 for index, vertex in enumerate(vertices)}

    lines = ["1"]
    for simplex_id in complex_.filtration_order:
        simplex = complex_.vertices(simplex_id)
        coordinates = " ".join(str(vertex_coordinates[vertex]) for vertex in simplex)
        birth = level_to_birth[complex_.filtration(simplex_id)]
        lines.append(f"{len(simplex) - 1} {coordinates} {birth}")

    path.write_text("\n".join(lines) + "\n")
    return birth_to_level


def _perseus_output_path(prefix: Path, dimension: int) -> Path:
    return Path(f"{prefix}_{dimension}.txt")


def _read_perseus_birth(
    birth_index: int,
    birth_to_level: dict[int, float],
    *,
    path: Path,
) -> float:
    try:
        return birth_to_level[birth_index]
    except KeyError as exc:
        supported = ", ".join(str(index) for index in sorted(birth_to_level))
        raise ValueError(
            f"Perseus output {path} uses unknown birth time {birth_index}; "
            f"expected one of: {supported}"
        ) from exc


def read_perseus_barcode(
    prefix: Path,
    *,
    birth_to_level: dict[int, float],
    max_dimension: int,
) -> tuple[tuple[tuple[int, float, float], ...], tuple[tuple[int, float], ...]]:
    finite: list[tuple[int, float, float]] = []
    essential: list[tuple[int, float]] = []

    for dimension in range(max_dimension + 1):
        path = _perseus_output_path(prefix, dimension)
        if not path.exists():
            continue
        for line_number, line in enumerate(path.read_text().splitlines(), start=1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) != 2:
                raise ValueError(
                    f"Expected two columns in Perseus output {path}:{line_number}, "
                    f"got {len(parts)}."
                )
            birth_index = int(parts[0])
            death_index = int(parts[1])
            birth = _read_perseus_birth(birth_index, birth_to_level, path=path)
            if death_index == -1:
                essential.append((dimension, birth))
            else:
                death = _read_perseus_birth(death_index, birth_to_level, path=path)
                if birth < death:
                    finite.append((dimension, birth, death))

    return tuple(sorted(finite)), tuple(sorted(essential))


def _remove_perseus_outputs(prefix: Path) -> None:
    for path in prefix.parent.glob(f"{prefix.name}_*.txt"):
        path.unlink()


def time_perseus_persistence(
    complex_: mp.FilteredComplex,
    *,
    repeats: int = 3,
    executable: str | Path = "perseus",
) -> PerseusTiming | None:
    resolved = _resolve_perseus_executable(executable)
    if resolved is None:
        return None
    if repeats <= 0:
        raise ValueError("repeats must be positive")

    max_dimension = max(record.dimension for record in complex_.simplices())
    standard = mp.compute_standard_persistence(complex_)
    standard_finite = standard.finite_barcode()
    standard_essential = standard.essential_barcode()

    with tempfile.TemporaryDirectory(prefix="morse-perseus-") as tmpdir:
        tmp = Path(tmpdir)
        input_path = tmp / "input_nmfsimtop.txt"
        output_prefix = tmp / "perseus"
        birth_to_level = write_perseus_nmfsimtop(complex_, input_path)

        best = math.inf
        finite: tuple[tuple[int, float, float], ...] = ()
        essential: tuple[tuple[int, float], ...] = ()
        for _ in range(repeats):
            _remove_perseus_outputs(output_prefix)
            started = perf_counter()
            completed = subprocess.run(
                [resolved, "nmfsimtop", str(input_path), str(output_prefix)],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            best = min(best, perf_counter() - started)
            if completed.returncode != 0:
                message = completed.stderr.strip() or completed.stdout.strip()
                raise RuntimeError(f"Perseus failed with exit code {completed.returncode}: {message}")
            finite, essential = read_perseus_barcode(
                output_prefix,
                birth_to_level=birth_to_level,
                max_dimension=max_dimension,
            )

    matches_standard = finite == standard_finite and essential == standard_essential
    if not matches_standard:
        raise AssertionError(
            f"Perseus barcode {finite!r}, {essential!r} "
            f"!= standard {standard_finite!r}, {standard_essential!r}"
        )

    return PerseusTiming(
        seconds=best,
        finite_barcode=finite,
        essential_barcode=essential,
        matches_standard=matches_standard,
    )


def sequence_algorithm_options(sequence_algorithm: str | Iterable[str]) -> tuple[str, ...]:
    if isinstance(sequence_algorithm, str):
        normalized = sequence_algorithm.lower().replace("_", "-")
        if normalized in {"all", "portfolio", mp.AUTO_MORSE_SEQUENCE_ALGORITHM}:
            return mp.DEFAULT_MORSE_ALGORITHM_PORTFOLIO
        return (sequence_algorithm,)

    values = tuple(sequence_algorithm)
    if not values:
        raise ValueError("At least one sequence algorithm is required.")
    return values


def profile_complex_shape(complex_: mp.FilteredComplex) -> ComplexShapeProfile:
    records = tuple(complex_.simplices())
    num_simplices = len(records)
    dimension_counts = Counter(record.dimension for record in records)
    num_vertices = dimension_counts.get(0, 0)
    num_edges = dimension_counts.get(1, 0)
    num_triangles = dimension_counts.get(2, 0)
    max_dimension = max(dimension_counts, default=0)
    num_higher_simplices = sum(
        count for dimension, count in dimension_counts.items() if dimension > 2
    )

    boundary_sizes = [len(record.boundary) for record in records]
    coboundary_sizes = [len(record.coboundary) for record in records]
    level_sizes = [len(complex_.simplices_of_level(level)) for level in range(complex_.num_levels)]
    largest_level_size = max(level_sizes, default=0)
    level_probabilities = [
        float(size) / float(num_simplices)
        for size in level_sizes
        if num_simplices and size
    ]
    level_entropy = -sum(probability * math.log(probability) for probability in level_probabilities)
    max_level_entropy = math.log(len(level_probabilities)) if len(level_probabilities) > 1 else 0.0

    edge_denominator = math.comb(num_vertices, 2) if num_vertices >= 2 else 0
    triangle_denominator = math.comb(num_vertices, 3) if num_vertices >= 3 else 0
    return ComplexShapeProfile(
        max_dimension=max_dimension,
        num_vertices=num_vertices,
        num_edges=num_edges,
        num_triangles=num_triangles,
        num_higher_simplices=num_higher_simplices,
        simplex_vertex_ratio=(
            float(num_simplices) / float(num_vertices) if num_vertices else 0.0
        ),
        vertex_sqrt=math.sqrt(float(num_vertices)) if num_vertices else 0.0,
        edge_density=(float(num_edges) / float(edge_denominator)) if edge_denominator else 0.0,
        triangle_density=(
            float(num_triangles) / float(triangle_denominator)
            if triangle_denominator
            else 0.0
        ),
        avg_boundary_size=(
            sum(boundary_sizes) / float(num_simplices) if num_simplices else 0.0
        ),
        max_boundary_size=max(boundary_sizes, default=0),
        avg_coboundary_size=(
            sum(coboundary_sizes) / float(num_simplices) if num_simplices else 0.0
        ),
        max_coboundary_size=max(coboundary_sizes, default=0),
        largest_level_size=largest_level_size,
        avg_level_size=(
            float(num_simplices) / float(len(level_sizes)) if level_sizes else 0.0
        ),
        largest_level_ratio=(
            float(largest_level_size) / float(num_simplices) if num_simplices else 0.0
        ),
        singleton_level_ratio=(
            sum(1 for size in level_sizes if size == 1) / float(len(level_sizes))
            if level_sizes
            else 0.0
        ),
        level_count_ratio=(
            float(complex_.num_levels) / float(num_simplices) if num_simplices else 0.0
        ),
        level_concentration=sum(probability * probability for probability in level_probabilities),
        level_entropy_ratio=(
            level_entropy / max_level_entropy if max_level_entropy > 0.0 else 0.0
        ),
    )


def shape_bin_from_profile(shape: ComplexShapeProfile) -> str:
    return shape_bin_from_values(
        max_dimension=shape.max_dimension,
        edge_density=shape.edge_density,
        triangle_density=shape.triangle_density,
        largest_level_ratio=shape.largest_level_ratio,
        max_coboundary=shape.max_coboundary_size,
    )


def read_profile_candidate_table(path: Path) -> dict[str, tuple[str, ...]]:
    table: dict[str, tuple[str, ...]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or ())
        missing = sorted({"group", "proposed_candidates"} - fieldnames)
        if missing:
            raise ValueError(
                f"Profile candidate table {path} is missing columns: {', '.join(missing)}"
            )
        for row in reader:
            group = row.get("group", "").strip()
            candidates = tuple(
                candidate.strip()
                for candidate in row.get("proposed_candidates", "").split(";")
                if candidate.strip()
            )
            if group and candidates:
                table[group] = candidates
    return table


def profile_candidate_table_keys(
    *,
    family: str,
    size: int,
    shape: ComplexShapeProfile | None,
) -> tuple[str, ...]:
    keys: list[str] = []
    if shape is not None:
        bin_name = shape_bin_from_profile(shape)
        keys.append(f"{family}:{bin_name}")
        keys.append(bin_name)
    keys.extend((f"{family}:n={size}", family, "all"))
    return tuple(keys)


def profile_candidate_algorithm_options(
    *,
    family: str,
    size: int,
    algorithms: tuple[str, ...],
    gate: str,
    shape: ComplexShapeProfile | None = None,
    candidate_table: dict[str, tuple[str, ...]] | None = None,
) -> tuple[str, ...]:
    """Return the algorithms profiled by an experimental selector gate.

    The measured comparison still runs the full ``algorithms`` tuple.  This gate only
    reduces the candidates seen by the cheap profile selector, so the resulting penalty
    remains measured against the full portfolio.
    """

    normalized_gate = gate.lower().replace("_", "-")
    if normalized_gate == "all":
        return algorithms
    if normalized_gate not in {
        "family-aware",
        "family",
        "adaptive-family",
        "shape-aware",
        "shape",
        "complex-shape",
        "tabulated",
        "table",
        "calibrated",
        "learned",
    }:
        supported = ", ".join(PROFILE_CANDIDATE_GATES)
        raise ValueError(f"Unknown profile candidate gate {gate!r}. Supported: {supported}.")

    if normalized_gate in {"tabulated", "table", "calibrated", "learned"}:
        if candidate_table is None:
            raise ValueError(
                "The tabulated profile candidate gate requires a profile candidate table."
            )
        preferred = algorithms
        for key in profile_candidate_table_keys(family=family, size=size, shape=shape):
            if key in candidate_table:
                preferred = candidate_table[key]
                break
    elif normalized_gate in {"shape-aware", "shape", "complex-shape"} and shape is not None:
        shape_size = (
            shape.num_vertices
            + shape.num_edges
            + shape.num_triangles
            + shape.num_higher_simplices
        )
        if (
            shape.max_dimension == 2
            and shape.edge_density >= 0.9
            and shape.triangle_density >= 0.9
        ):
            preferred = (mp.F_MAX_SEQUENCE, mp.COREDUCTION_SEQUENCE)
        elif shape_size <= 250:
            preferred = (
                mp.F_MAX_SEQUENCE,
                mp.F_MIN_SEQUENCE,
                mp.SATURATED_SEQUENCE,
                mp.COREDUCTION_SEQUENCE,
            )
        elif shape.max_dimension >= 3:
            preferred = (mp.F_MAX_SEQUENCE, mp.F_MIN_SEQUENCE, mp.SATURATED_SEQUENCE)
        elif shape.num_edges >= 5000 or shape.max_coboundary_size >= 100:
            preferred = (mp.F_MAX_SEQUENCE, mp.COREDUCTION_SEQUENCE)
        elif (
            shape.edge_density >= 0.45
            or shape.triangle_density >= 0.08
            or shape.largest_level_ratio >= 0.5
        ):
            preferred = (mp.F_MAX_SEQUENCE, mp.F_MIN_SEQUENCE, mp.COREDUCTION_SEQUENCE)
        else:
            preferred = (mp.F_MAX_SEQUENCE, mp.F_MIN_SEQUENCE, mp.SATURATED_SEQUENCE)
    elif family in {"lower-star", "plateau", "terrain", "image-grid"}:
        preferred = (mp.F_MAX_SEQUENCE, mp.F_MIN_SEQUENCE, mp.SATURATED_SEQUENCE)
    elif family == "roadmap-rips":
        preferred = (mp.F_MAX_SEQUENCE, mp.COREDUCTION_SEQUENCE)
    elif family in {"rips", "cam-s4-rips"}:
        preferred = (mp.F_MAX_SEQUENCE, mp.F_MIN_SEQUENCE, mp.COREDUCTION_SEQUENCE)
    else:
        preferred = algorithms

    available = set(algorithms)
    gated = tuple(algorithm for algorithm in preferred if algorithm in available)
    return gated if gated else algorithms


def make_benchmark_complex(
    family: str,
    seed: int,
    size: int,
    *,
    plateau_levels: int = 3,
    roadmap_cache: Path = DEFAULT_ROADMAP_CACHE,
    download_roadmap_data: bool = False,
) -> tuple[mp.FilteredComplex, str]:
    if family == "lower-star":
        complex_ = make_lower_star(seed, size)
        name = f"{family}-n{size}-seed{seed}"
    elif family == "plateau":
        complex_ = make_plateau(seed, size, levels=plateau_levels)
        name = f"{family}-n{size}-L{plateau_levels}-seed{seed}"
    elif family == "terrain":
        complex_ = make_terrain(seed, size, levels=plateau_levels)
        name = f"{family}-g{size}-L{plateau_levels}-seed{seed}"
    elif family == "image-grid":
        complex_ = make_image_grid(seed, size, levels=plateau_levels)
        name = f"{family}-g{size}-L{plateau_levels}-seed{seed}"
    elif family == "rips":
        complex_ = make_rips(seed, size)
        name = f"{family}-n{size}-seed{seed}"
    elif family == "cam-s4-rips":
        complex_ = make_cam_s4_rips(seed, size)
        name = f"{family}-n{size}-seed{seed}"
    elif family == "roadmap-rips":
        complex_ = make_roadmap_rips(
            size,
            cache_dir=roadmap_cache,
            download=download_roadmap_data,
        )
        name = f"{family}-{ROADMAP_DATASETS[size].label}"
    else:
        raise ValueError(f"Unknown benchmark family: {family}")
    return complex_, name


def benchmark_complex_case(
    family: str,
    name: str,
    seed: int,
    size: int,
    complex_: mp.FilteredComplex,
    repeats: int,
    *,
    time_gudhi: bool = False,
    time_perseus: bool = False,
    perseus_executable: str | Path = "perseus",
    sequence_algorithm: str = mp.DEFAULT_MORSE_SEQUENCE_ALGORITHM,
    frame_mode: str = mp.DEFAULT_MORSE_FRAME_MODE,
    validation_mode: str = "core",
    gudhi_cam_seconds_override: float | None | object = _EXTERNAL_TIMING_UNSET,
    perseus_result_override: PerseusTiming | None | object = _EXTERNAL_TIMING_UNSET,
) -> BenchmarkRow:
    result = mp.benchmark_persistence(
        complex_,
        repeats=repeats,
        sequence_algorithm=sequence_algorithm,
        frame_mode=frame_mode,
        materialize_barcodes=validation_mode == "materialized",
    )
    if gudhi_cam_seconds_override is _EXTERNAL_TIMING_UNSET:
        gudhi_cam_seconds = time_gudhi_cam(complex_, repeats=repeats) if time_gudhi else None
    else:
        gudhi_cam_seconds = gudhi_cam_seconds_override
    gudhi_cam_speedup = (
        None
        if gudhi_cam_seconds is None or result.morse_seconds == 0.0
        else gudhi_cam_seconds / result.morse_seconds
    )
    if perseus_result_override is _EXTERNAL_TIMING_UNSET:
        perseus_result = (
            time_perseus_persistence(
                complex_,
                repeats=repeats,
                executable=perseus_executable,
            )
            if time_perseus
            else None
        )
    else:
        perseus_result = perseus_result_override
    perseus_seconds = None if perseus_result is None else perseus_result.seconds
    perseus_speedup = (
        None
        if perseus_seconds is None or result.morse_seconds == 0.0
        else perseus_seconds / result.morse_seconds
    )
    return BenchmarkRow(
        family=family,
        name=name,
        seed=seed,
        size=size,
        num_simplices=result.num_simplices,
        num_levels=result.num_levels,
        num_critical_simplices=result.num_critical_simplices,
        sequence_algorithm=result.sequence_algorithm,
        frame_mode=result.frame_mode,
        validation_mode=result.validation_mode,
        barcodes_materialized=result.barcodes_materialized,
        reference_final_live_nonempty_annotations=(
            result.reference_final_live_nonempty_annotations
        ),
        reference_final_live_total_annotation_size=(
            result.reference_final_live_total_annotation_size
        ),
        reference_peak_live_nonempty_annotations=(
            result.reference_peak_live_nonempty_annotations
        ),
        reference_peak_live_total_annotation_size=(
            result.reference_peak_live_total_annotation_size
        ),
        reference_released_annotations=result.reference_released_annotations,
        reference_released_total_annotation_size=(
            result.reference_released_total_annotation_size
        ),
        reducer_working_set_size=result.reducer_working_set_size,
        reducer_initial_nonempty_annotations=result.reducer_initial_nonempty_annotations,
        reducer_initial_total_annotation_size=result.reducer_initial_total_annotation_size,
        reducer_initial_max_annotation_size=result.reducer_initial_max_annotation_size,
        reducer_initial_inverse_list_entries=result.reducer_initial_inverse_list_entries,
        reducer_boundary_plan_face_scans=result.reducer_boundary_plan_face_scans,
        reducer_boundary_annotation_candidate_criticals=(
            result.reducer_boundary_annotation_candidate_criticals
        ),
        reducer_boundary_annotation_zero_skipped_criticals=(
            result.reducer_boundary_annotation_zero_skipped_criticals
        ),
        reducer_boundary_annotation_zero_skipped_faces=(
            result.reducer_boundary_annotation_zero_skipped_faces
        ),
        reducer_boundary_annotation_xors=result.reducer_boundary_annotation_xors,
        reducer_boundary_annotation_total_input_size=(
            result.reducer_boundary_annotation_total_input_size
        ),
        reducer_boundary_annotation_total_output_size=(
            result.reducer_boundary_annotation_total_output_size
        ),
        reducer_boundary_annotation_max_size=result.reducer_boundary_annotation_max_size,
        reducer_boundary_annotation_max_output_size=(
            result.reducer_boundary_annotation_max_output_size
        ),
        reducer_pivot_eliminations=result.reducer_pivot_eliminations,
        reducer_remove_candidate_scans=result.reducer_remove_candidate_scans,
        reducer_remove_applied=result.reducer_remove_applied,
        reducer_remove_total_annotation_size=result.reducer_remove_total_annotation_size,
        reducer_remove_max_annotation_size=result.reducer_remove_max_annotation_size,
        reducer_xor_candidate_scans=result.reducer_xor_candidate_scans,
        reducer_xor_applied=result.reducer_xor_applied,
        reducer_xor_changed_labels=result.reducer_xor_changed_labels,
        reducer_xor_total_input_size=result.reducer_xor_total_input_size,
        reducer_xor_total_output_size=result.reducer_xor_total_output_size,
        reducer_xor_max_input_size=result.reducer_xor_max_input_size,
        reducer_xor_max_output_size=result.reducer_xor_max_output_size,
        reducer_xor_inserted_labels=result.reducer_xor_inserted_labels,
        reducer_xor_removed_labels=result.reducer_xor_removed_labels,
        reducer_inverse_list_appends=result.reducer_inverse_list_appends,
        finite_intervals=result.finite_interval_count,
        essential_intervals=result.essential_interval_count,
        repeats=result.repeats,
        sequence_seconds=result.sequence_seconds,
        reference_seconds=result.reference_seconds,
        morse_reduction_seconds=result.morse_reduction_seconds,
        reducer_setup_seconds=result.reducer_setup_seconds,
        reducer_compute_seconds=result.reducer_compute_seconds,
        morse_seconds=result.morse_seconds,
        standard_seconds=result.standard_seconds,
        gudhi_cam_seconds=gudhi_cam_seconds,
        gudhi_cam_speedup=gudhi_cam_speedup,
        perseus_seconds=perseus_seconds,
        perseus_speedup=perseus_speedup,
        perseus_matches=None if perseus_result is None else perseus_result.matches_standard,
        speedup=result.speedup,
        cpp_backend=complex_.cpp_backend_active(),
    )


def benchmark_case(
    family: str,
    seed: int,
    size: int,
    repeats: int,
    *,
    plateau_levels: int = 3,
    time_gudhi: bool = False,
    time_perseus: bool = False,
    perseus_executable: str | Path = "perseus",
    sequence_algorithm: str = mp.DEFAULT_MORSE_SEQUENCE_ALGORITHM,
    frame_mode: str = mp.DEFAULT_MORSE_FRAME_MODE,
    validation_mode: str = "core",
    roadmap_cache: Path = DEFAULT_ROADMAP_CACHE,
    download_roadmap_data: bool = False,
) -> BenchmarkRow:
    complex_, name = make_benchmark_complex(
        family,
        seed,
        size,
        plateau_levels=plateau_levels,
        roadmap_cache=roadmap_cache,
        download_roadmap_data=download_roadmap_data,
    )
    return benchmark_complex_case(
        family,
        name,
        seed,
        size,
        complex_,
        repeats,
        time_gudhi=time_gudhi,
        time_perseus=time_perseus,
        perseus_executable=perseus_executable,
        sequence_algorithm=sequence_algorithm,
        frame_mode=frame_mode,
        validation_mode=validation_mode,
    )


def run_benchmarks(
    *,
    families: Iterable[str],
    sizes: Iterable[int],
    seeds: Iterable[int],
    repeats: int,
    plateau_levels: int = 3,
    time_gudhi: bool = False,
    time_perseus: bool = False,
    perseus_executable: str | Path = "perseus",
    sequence_algorithm: str | Iterable[str] = mp.DEFAULT_MORSE_SEQUENCE_ALGORITHM,
    frame_mode: str = mp.DEFAULT_MORSE_FRAME_MODE,
    validation_mode: str = "core",
    roadmap_cache: Path = DEFAULT_ROADMAP_CACHE,
    download_roadmap_data: bool = False,
) -> list[BenchmarkRow]:
    rows: list[BenchmarkRow] = []
    algorithms = sequence_algorithm_options(sequence_algorithm)
    for family in families:
        for size in sizes:
            for seed in seeds:
                complex_, name = make_benchmark_complex(
                    family,
                    seed,
                    size,
                    plateau_levels=plateau_levels,
                    roadmap_cache=roadmap_cache,
                    download_roadmap_data=download_roadmap_data,
                )
                gudhi_cam_seconds = (
                    time_gudhi_cam(complex_, repeats=repeats)
                    if time_gudhi
                    else None
                )
                perseus_result = (
                    time_perseus_persistence(
                        complex_,
                        repeats=repeats,
                        executable=perseus_executable,
                    )
                    if time_perseus
                    else None
                )
                for algorithm in algorithms:
                    rows.append(
                        benchmark_complex_case(
                            family,
                            name,
                            seed,
                            size,
                            complex_,
                            repeats,
                            time_gudhi=time_gudhi,
                            time_perseus=time_perseus,
                            perseus_executable=perseus_executable,
                            sequence_algorithm=algorithm,
                            frame_mode=frame_mode,
                            validation_mode=validation_mode,
                            gudhi_cam_seconds_override=gudhi_cam_seconds,
                            perseus_result_override=perseus_result,
                        )
                    )
    return rows


def profile_selection_key(
    profile: mp.MorseReferenceProfile,
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
    supported = ", ".join(mp.MORSE_PROFILE_SELECTION_METRICS)
    raise ValueError(f"Unknown profile selection metric {selection_metric!r}. Supported: {supported}.")


def resolve_profile_selection_metric(family: str, selection_metric: str) -> str:
    metric = selection_metric.lower().replace("_", "-")
    if metric in {"adaptive", "family-adaptive", "critical-except-roadmap"}:
        if family == "roadmap-rips":
            return "profile_total_work"
        return "critical_ratio"
    if metric in {
        "adaptive-structured",
        "structured-adaptive",
        "grid-adaptive",
        "critical-except-structured",
    }:
        if family in {"roadmap-rips", "image-grid"}:
            return "profile_total_work"
        return "critical_ratio"
    return selection_metric


def profile_selection_metric_options(
    profile_selection_metric: str | Iterable[str],
) -> tuple[str, ...]:
    if isinstance(profile_selection_metric, str):
        values = tuple(
            metric.strip()
            for metric in profile_selection_metric.replace(",", " ").split()
            if metric.strip()
        )
        return values if values else (profile_selection_metric,)

    values = tuple(profile_selection_metric)
    if not values:
        raise ValueError("At least one profile selection metric is required.")
    return values


def benchmark_row_reducer_work(row: BenchmarkRow) -> float:
    return float(
        row.reducer_initial_total_annotation_size
        + row.reducer_boundary_annotation_total_input_size
        + row.reducer_boundary_annotation_total_output_size
        + row.reducer_xor_total_input_size
        + row.reducer_xor_total_output_size
        + row.reducer_remove_total_annotation_size
    )


def benchmark_row_critical_ratio(row: BenchmarkRow) -> float:
    if row.num_simplices == 0:
        return 0.0
    return float(row.num_critical_simplices) / float(row.num_simplices)


def benchmark_row_selection_key(
    row: BenchmarkRow,
    selection_metric: str,
) -> tuple[float, float, float]:
    metric = selection_metric.lower().replace("_", "-")
    critical_ratio = benchmark_row_critical_ratio(row)
    reducer_work = benchmark_row_reducer_work(row)
    if metric in {"morse-seconds", "morse-seconds-best", "time"}:
        return (row.morse_seconds, reducer_work, critical_ratio)
    if metric in {"critical-ratio", "critical-count", "critical"}:
        return (critical_ratio, row.morse_seconds, reducer_work)
    if metric in {"reducer-work", "annotation-work", "work"}:
        return (reducer_work, row.morse_seconds, critical_ratio)
    supported = ", ".join(mp.MORSE_ALGORITHM_SELECTION_METRICS)
    raise ValueError(f"Unknown measured selection metric {selection_metric!r}. Supported: {supported}.")


def _profile_penalty_ratio(selected_seconds: float, best_seconds: float) -> float:
    if best_seconds == 0.0:
        return 1.0 if selected_seconds == 0.0 else math.inf
    return selected_seconds / best_seconds


def run_profile_vs_measured(
    *,
    families: Iterable[str],
    sizes: Iterable[int],
    seeds: Iterable[int],
    repeats: int,
    plateau_levels: int = 3,
    profile_repeats: int = 1,
    time_gudhi: bool = False,
    time_perseus: bool = False,
    perseus_executable: str | Path = "perseus",
    sequence_algorithm: str | Iterable[str] = "portfolio",
    frame_mode: str = mp.DEFAULT_MORSE_FRAME_MODE,
    validation_mode: str = "core",
    profile_selection_metric: str | Iterable[str] = "estimated_reducer_work",
    profile_candidate_gate: str = "all",
    profile_candidate_table: Path | None = None,
    measured_selection_metric: str = "morse_seconds",
    roadmap_cache: Path = DEFAULT_ROADMAP_CACHE,
    download_roadmap_data: bool = False,
) -> list[ProfileVsMeasuredRow]:
    if profile_repeats < 1:
        raise ValueError("profile_repeats must be positive")

    rows: list[ProfileVsMeasuredRow] = []
    algorithms = sequence_algorithm_options(sequence_algorithm)
    profile_selection_metrics = profile_selection_metric_options(profile_selection_metric)
    candidate_table = (
        read_profile_candidate_table(profile_candidate_table)
        if profile_candidate_table is not None
        else None
    )
    for family in families:
        for size in sizes:
            for seed in seeds:
                complex_, name = make_benchmark_complex(
                    family,
                    seed,
                    size,
                    plateau_levels=plateau_levels,
                    roadmap_cache=roadmap_cache,
                    download_roadmap_data=download_roadmap_data,
                )
                shape = profile_complex_shape(complex_)
                profile_algorithms = profile_candidate_algorithm_options(
                    family=family,
                    size=size,
                    algorithms=algorithms,
                    gate=profile_candidate_gate,
                    shape=shape,
                    candidate_table=candidate_table,
                )
                profiles = mp.profile_morse_sequence_algorithms(
                    complex_,
                    algorithms=profile_algorithms,
                    repeats=profile_repeats,
                )
                profile_by_algorithm = {
                    profile.sequence_algorithm: profile
                    for profile in profiles
                }

                gudhi_cam_seconds = (
                    time_gudhi_cam(complex_, repeats=repeats)
                    if time_gudhi
                    else None
                )
                perseus_result = (
                    time_perseus_persistence(
                        complex_,
                        repeats=repeats,
                        executable=perseus_executable,
                    )
                    if time_perseus
                    else None
                )

                measured_rows = [
                    benchmark_complex_case(
                        family,
                        name,
                        seed,
                        size,
                        complex_,
                        repeats,
                        time_gudhi=time_gudhi,
                        time_perseus=time_perseus,
                        perseus_executable=perseus_executable,
                        sequence_algorithm=algorithm,
                        frame_mode=frame_mode,
                        validation_mode=validation_mode,
                        gudhi_cam_seconds_override=gudhi_cam_seconds,
                        perseus_result_override=perseus_result,
                    )
                    for algorithm in algorithms
                ]
                measured_by_algorithm = {
                    row.sequence_algorithm: row
                    for row in measured_rows
                }
                measured_best = min(
                    measured_rows,
                    key=lambda row: benchmark_row_selection_key(row, measured_selection_metric),
                )

                measured_ranking = sorted(
                    measured_rows,
                    key=lambda row: benchmark_row_selection_key(row, measured_selection_metric),
                )
                for requested_profile_selection_metric in profile_selection_metrics:
                    effective_profile_selection_metric = resolve_profile_selection_metric(
                        family,
                        requested_profile_selection_metric,
                    )
                    profile_selected = min(
                        profiles,
                        key=lambda profile: profile_selection_key(
                            profile,
                            effective_profile_selection_metric,
                        ),
                    )
                    try:
                        profile_selected_row = measured_by_algorithm[
                            profile_selected.sequence_algorithm
                        ]
                    except KeyError as exc:
                        raise RuntimeError(
                            "Profile and measured benchmark algorithms did not align."
                        ) from exc

                    profile_ranking = sorted(
                        profiles,
                        key=lambda profile: profile_selection_key(
                            profile,
                            effective_profile_selection_metric,
                        ),
                    )
                    profile_selected_rank = (
                        next(
                            index
                            for index, row in enumerate(measured_ranking, start=1)
                            if row.sequence_algorithm == profile_selected.sequence_algorithm
                        )
                    )
                    measured_best_profile_rank = (
                        next(
                            (
                                index
                                for index, profile in enumerate(profile_ranking, start=1)
                                if profile.sequence_algorithm == measured_best.sequence_algorithm
                            ),
                            0,
                        )
                    )
                    measured_best_profile = profile_by_algorithm.get(
                        measured_best.sequence_algorithm
                    )

                    ratio = _profile_penalty_ratio(
                        profile_selected_row.morse_seconds,
                        measured_best.morse_seconds,
                    )
                    penalty_percent = math.inf if math.isinf(ratio) else 100.0 * (ratio - 1.0)

                    rows.append(
                        ProfileVsMeasuredRow(
                        family=family,
                        name=name,
                        seed=seed,
                        size=size,
                        num_simplices=measured_best.num_simplices,
                        num_levels=measured_best.num_levels,
                        shape_max_dimension=shape.max_dimension,
                        shape_num_vertices=shape.num_vertices,
                        shape_num_edges=shape.num_edges,
                        shape_num_triangles=shape.num_triangles,
                        shape_num_higher_simplices=shape.num_higher_simplices,
                        shape_simplex_vertex_ratio=shape.simplex_vertex_ratio,
                        shape_vertex_sqrt=shape.vertex_sqrt,
                        shape_edge_density=shape.edge_density,
                        shape_triangle_density=shape.triangle_density,
                        shape_avg_boundary_size=shape.avg_boundary_size,
                        shape_avg_coboundary_size=shape.avg_coboundary_size,
                        shape_max_coboundary_size=shape.max_coboundary_size,
                        shape_largest_level_size=shape.largest_level_size,
                        shape_avg_level_size=shape.avg_level_size,
                        shape_largest_level_ratio=shape.largest_level_ratio,
                        shape_singleton_level_ratio=shape.singleton_level_ratio,
                        shape_level_count_ratio=shape.level_count_ratio,
                        shape_level_concentration=shape.level_concentration,
                        shape_level_entropy_ratio=shape.level_entropy_ratio,
                        candidate_count=len(measured_rows),
                        profile_candidate_count=len(profiles),
                        profile_candidate_gate=profile_candidate_gate,
                        frame_mode=frame_mode,
                        validation_mode=validation_mode,
                        profile_selection_metric=requested_profile_selection_metric,
                        effective_profile_selection_metric=effective_profile_selection_metric,
                        measured_selection_metric=measured_selection_metric,
                        profile_selected_algorithm=profile_selected.sequence_algorithm,
                        measured_best_algorithm=measured_best.sequence_algorithm,
                        profile_matches_measured=(
                            profile_selected.sequence_algorithm
                            == measured_best.sequence_algorithm
                        ),
                        profile_selected_measured_rank=profile_selected_rank,
                        measured_best_profile_rank=measured_best_profile_rank,
                        measured_best_morse_seconds=measured_best.morse_seconds,
                        profile_selected_morse_seconds=profile_selected_row.morse_seconds,
                        profile_penalty_seconds=(
                            profile_selected_row.morse_seconds - measured_best.morse_seconds
                        ),
                        profile_penalty_ratio=ratio,
                        profile_penalty_percent=penalty_percent,
                        measured_best_critical_ratio=benchmark_row_critical_ratio(measured_best),
                        profile_selected_critical_ratio=benchmark_row_critical_ratio(
                            profile_selected_row
                        ),
                        measured_best_profile_score=(
                            None
                            if measured_best_profile is None
                            else profile_selection_key(
                                measured_best_profile,
                                effective_profile_selection_metric,
                            )[0]
                        ),
                        profile_selected_profile_score=profile_selection_key(
                            profile_selected,
                            effective_profile_selection_metric,
                        )[0],
                        measured_best_estimated_reducer_work=(
                            None
                            if measured_best_profile is None
                            else measured_best_profile.estimated_reducer_work
                        ),
                        profile_selected_estimated_reducer_work=(
                            profile_selected.estimated_reducer_work
                        ),
                        measured_best_profile_seconds=(
                            None
                            if measured_best_profile is None
                            else measured_best_profile.profile_seconds
                        ),
                        profile_selected_profile_seconds=profile_selected.profile_seconds,
                        measured_best_sequence_seconds=measured_best.sequence_seconds,
                        profile_selected_sequence_seconds=(
                            profile_selected_row.sequence_seconds
                        ),
                        measured_best_reducer_compute_seconds=(
                            measured_best.reducer_compute_seconds
                        ),
                        profile_selected_reducer_compute_seconds=(
                            profile_selected_row.reducer_compute_seconds
                        ),
                        measured_best_reducer_work=benchmark_row_reducer_work(
                            measured_best
                        ),
                        profile_selected_reducer_work=benchmark_row_reducer_work(
                            profile_selected_row
                        ),
                        total_profile_seconds=sum(
                            profile.profile_seconds for profile in profiles
                        ),
                        total_measured_morse_seconds=sum(
                            row.morse_seconds for row in measured_rows
                        ),
                        standard_seconds=measured_best.standard_seconds,
                        gudhi_cam_seconds=measured_best.gudhi_cam_seconds,
                        measured_best_gudhi_cam_speedup=measured_best.gudhi_cam_speedup,
                        profile_selected_gudhi_cam_speedup=(
                            profile_selected_row.gudhi_cam_speedup
                        ),
                        perseus_seconds=measured_best.perseus_seconds,
                        measured_best_perseus_speedup=measured_best.perseus_speedup,
                        profile_selected_perseus_speedup=(
                            profile_selected_row.perseus_speedup
                        ),
                        perseus_matches=measured_best.perseus_matches,
                        cpp_backend=complex_.cpp_backend_active(),
                    )
                )
    return rows


def summarize_rows(rows: list[BenchmarkRow]) -> list[BenchmarkSummary]:
    grouped: dict[tuple[str, str, str, str, int], list[BenchmarkRow]] = {}
    for row in rows:
        grouped.setdefault(
            (row.family, row.sequence_algorithm, row.frame_mode, row.validation_mode, row.size),
            [],
        ).append(row)

    summaries: list[BenchmarkSummary] = []
    for (family, sequence_algorithm, frame_mode, validation_mode, size), group in sorted(
        grouped.items()
    ):
        cases = len(group)

        def average(attribute: str) -> float:
            return sum(float(getattr(row, attribute)) for row in group) / cases

        def average_optional(attribute: str) -> float | None:
            values = [getattr(row, attribute) for row in group]
            numeric = [float(value) for value in values if value is not None]
            if not numeric:
                return None
            return sum(numeric) / len(numeric)

        avg_num_simplices = average("num_simplices")
        avg_morse_seconds = average("morse_seconds")
        if avg_num_simplices == 0.0:
            avg_critical_ratio = 0.0
        else:
            avg_critical_ratio = average("num_critical_simplices") / avg_num_simplices

        if avg_morse_seconds == 0.0:
            avg_sequence_share = 0.0
            avg_reference_share = 0.0
            avg_reduction_share = 0.0
        else:
            avg_sequence_share = average("sequence_seconds") / avg_morse_seconds
            avg_reference_share = average("reference_seconds") / avg_morse_seconds
            avg_reduction_share = average("morse_reduction_seconds") / avg_morse_seconds

        speedups = [row.speedup for row in group]
        summaries.append(
            BenchmarkSummary(
                family=family,
                sequence_algorithm=sequence_algorithm,
                frame_mode=frame_mode,
                validation_mode=validation_mode,
                size=size,
                cases=cases,
                avg_num_simplices=avg_num_simplices,
                avg_num_levels=average("num_levels"),
                avg_num_critical_simplices=average("num_critical_simplices"),
                avg_critical_ratio=avg_critical_ratio,
                avg_reference_final_live_nonempty_annotations=average(
                    "reference_final_live_nonempty_annotations"
                ),
                avg_reference_final_live_total_annotation_size=average(
                    "reference_final_live_total_annotation_size"
                ),
                avg_reference_peak_live_nonempty_annotations=average(
                    "reference_peak_live_nonempty_annotations"
                ),
                avg_reference_peak_live_total_annotation_size=average(
                    "reference_peak_live_total_annotation_size"
                ),
                avg_reference_released_annotations=average("reference_released_annotations"),
                avg_reference_released_total_annotation_size=average(
                    "reference_released_total_annotation_size"
                ),
                avg_reducer_working_set_size=average("reducer_working_set_size"),
                avg_reducer_initial_nonempty_annotations=average(
                    "reducer_initial_nonempty_annotations"
                ),
                avg_reducer_initial_total_annotation_size=average(
                    "reducer_initial_total_annotation_size"
                ),
                avg_reducer_initial_max_annotation_size=average(
                    "reducer_initial_max_annotation_size"
                ),
                avg_reducer_initial_inverse_list_entries=average(
                    "reducer_initial_inverse_list_entries"
                ),
                avg_reducer_boundary_plan_face_scans=average(
                    "reducer_boundary_plan_face_scans"
                ),
                avg_reducer_boundary_annotation_candidate_criticals=average(
                    "reducer_boundary_annotation_candidate_criticals"
                ),
                avg_reducer_boundary_annotation_zero_skipped_criticals=average(
                    "reducer_boundary_annotation_zero_skipped_criticals"
                ),
                avg_reducer_boundary_annotation_zero_skipped_faces=average(
                    "reducer_boundary_annotation_zero_skipped_faces"
                ),
                avg_reducer_boundary_annotation_xors=average("reducer_boundary_annotation_xors"),
                avg_reducer_boundary_annotation_total_output_size=average(
                    "reducer_boundary_annotation_total_output_size"
                ),
                avg_reducer_pivot_eliminations=average("reducer_pivot_eliminations"),
                avg_reducer_remove_candidate_scans=average("reducer_remove_candidate_scans"),
                avg_reducer_xor_candidate_scans=average("reducer_xor_candidate_scans"),
                avg_reducer_xor_applied=average("reducer_xor_applied"),
                avg_reducer_xor_total_input_size=average("reducer_xor_total_input_size"),
                avg_reducer_xor_total_output_size=average("reducer_xor_total_output_size"),
                avg_reducer_xor_max_input_size=average("reducer_xor_max_input_size"),
                avg_reducer_xor_max_output_size=average("reducer_xor_max_output_size"),
                avg_reducer_xor_inserted_labels=average("reducer_xor_inserted_labels"),
                avg_reducer_xor_removed_labels=average("reducer_xor_removed_labels"),
                avg_sequence_share=avg_sequence_share,
                avg_reference_share=avg_reference_share,
                avg_reduction_share=avg_reduction_share,
                avg_reducer_setup_seconds=average("reducer_setup_seconds"),
                avg_reducer_compute_seconds=average("reducer_compute_seconds"),
                avg_morse_seconds=avg_morse_seconds,
                avg_standard_seconds=average("standard_seconds"),
                avg_gudhi_cam_seconds=average_optional("gudhi_cam_seconds"),
                avg_gudhi_cam_speedup=average_optional("gudhi_cam_speedup"),
                avg_perseus_seconds=average_optional("perseus_seconds"),
                avg_perseus_speedup=average_optional("perseus_speedup"),
                avg_speedup=average("speedup"),
                min_speedup=min(speedups),
                max_speedup=max(speedups),
            )
        )

    return summaries


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return "-"
    return f"{value:.6g}"


def _format_summary_optional_us(value: float | None) -> str:
    if value is None:
        return "-"
    return f"{1.0e6 * value:.2f}"


def write_table(rows: list[BenchmarkRow], output: TextIO) -> None:
    headers = [
        "family",
        "seq_alg",
        "frame",
        "valid",
        "seed",
        "size",
        "simplices",
        "critical",
        "ref_peak",
        "ref_rel",
        "workset",
        "init_non",
        "init_lbl",
        "plan_face",
        "cand",
        "skip0",
        "skipface",
        "bnd_xor",
        "rm_scan",
        "xor_scan",
        "xor_apply",
        "seq_s",
        "ref_s",
        "reduce_s",
        "setup_s",
        "loop_s",
        "morse_s",
        "standard_s",
        "cam_s",
        "cam/morse",
        "perseus_s",
        "perseus/morse",
        "perseus_ok",
        "speedup",
        "cpp",
    ]
    output.write(" ".join(f"{header:>12}" for header in headers) + "\n")
    for row in rows:
        output.write(
            f"{row.family:>12} "
            f"{row.sequence_algorithm:>12} "
            f"{row.frame_mode:>12} "
            f"{row.validation_mode:>12} "
            f"{row.seed:12d} "
            f"{row.size:12d} "
            f"{row.num_simplices:12d} "
            f"{row.num_critical_simplices:12d} "
            f"{row.reference_peak_live_total_annotation_size:12d} "
            f"{row.reference_released_total_annotation_size:12d} "
            f"{row.reducer_working_set_size:12d} "
            f"{row.reducer_initial_nonempty_annotations:12d} "
            f"{row.reducer_initial_total_annotation_size:12d} "
            f"{row.reducer_boundary_plan_face_scans:12d} "
            f"{row.reducer_boundary_annotation_candidate_criticals:12d} "
            f"{row.reducer_boundary_annotation_zero_skipped_criticals:12d} "
            f"{row.reducer_boundary_annotation_zero_skipped_faces:12d} "
            f"{row.reducer_boundary_annotation_xors:12d} "
            f"{row.reducer_remove_candidate_scans:12d} "
            f"{row.reducer_xor_candidate_scans:12d} "
            f"{row.reducer_xor_applied:12d} "
            f"{row.sequence_seconds:12.6g} "
            f"{row.reference_seconds:12.6g} "
            f"{row.morse_reduction_seconds:12.6g} "
            f"{row.reducer_setup_seconds:12.6g} "
            f"{row.reducer_compute_seconds:12.6g} "
            f"{row.morse_seconds:12.6g} "
            f"{row.standard_seconds:12.6g} "
            f"{_format_optional_float(row.gudhi_cam_seconds):>12} "
            f"{_format_optional_float(row.gudhi_cam_speedup):>12} "
            f"{_format_optional_float(row.perseus_seconds):>12} "
            f"{_format_optional_float(row.perseus_speedup):>12} "
            f"{('-' if row.perseus_matches is None else str(row.perseus_matches)):>12} "
            f"{row.speedup:12.6g} "
            f"{str(row.cpp_backend):>12}\n"
        )


def write_summary(rows: list[BenchmarkRow], output: TextIO) -> None:
    summaries = summarize_rows(rows)
    headers = [
        "family",
        "seq_alg",
        "frame",
        "valid",
        "size",
        "cases",
        "simp",
        "crit",
        "crit%",
        "speed",
        "min",
        "max",
        "seq%",
        "ref%",
        "red%",
        "refpeak",
        "refrel",
        "morse_us",
        "setup_us",
        "loop_us",
        "std_us",
        "cam_us",
        "cam/morse",
        "pers_us",
        "pers/morse",
        "workset",
        "init_non",
        "init_lbl",
        "init_max",
        "planface",
        "cand",
        "skip0",
        "skipface",
        "bnd_xor",
        "xor_scan",
        "bnd_out",
        "pivots",
        "xor_in",
        "xor_out",
        "xor_maxin",
        "xor_maxout",
        "xor_ins",
        "xor_rem",
    ]
    output.write(" ".join(f"{header:>10}" for header in headers) + "\n")
    for row in summaries:
        output.write(
            f"{row.family:>10} "
            f"{row.sequence_algorithm:>10} "
            f"{row.frame_mode:>10} "
            f"{row.validation_mode:>10} "
            f"{row.size:10d} "
            f"{row.cases:10d} "
            f"{row.avg_num_simplices:10.1f} "
            f"{row.avg_num_critical_simplices:10.1f} "
            f"{100.0 * row.avg_critical_ratio:10.1f} "
            f"{row.avg_speedup:10.3f} "
            f"{row.min_speedup:10.3f} "
            f"{row.max_speedup:10.3f} "
            f"{100.0 * row.avg_sequence_share:10.1f} "
            f"{100.0 * row.avg_reference_share:10.1f} "
            f"{100.0 * row.avg_reduction_share:10.1f} "
            f"{row.avg_reference_peak_live_total_annotation_size:10.1f} "
            f"{row.avg_reference_released_total_annotation_size:10.1f} "
            f"{1.0e6 * row.avg_morse_seconds:10.2f} "
            f"{1.0e6 * row.avg_reducer_setup_seconds:10.2f} "
            f"{1.0e6 * row.avg_reducer_compute_seconds:10.2f} "
            f"{1.0e6 * row.avg_standard_seconds:10.2f} "
            f"{_format_summary_optional_us(row.avg_gudhi_cam_seconds):>10} "
            f"{_format_optional_float(row.avg_gudhi_cam_speedup):>10} "
            f"{_format_summary_optional_us(row.avg_perseus_seconds):>10} "
            f"{_format_optional_float(row.avg_perseus_speedup):>10} "
            f"{row.avg_reducer_working_set_size:10.1f} "
            f"{row.avg_reducer_initial_nonempty_annotations:10.1f} "
            f"{row.avg_reducer_initial_total_annotation_size:10.1f} "
            f"{row.avg_reducer_initial_max_annotation_size:10.1f} "
            f"{row.avg_reducer_boundary_plan_face_scans:10.1f} "
            f"{row.avg_reducer_boundary_annotation_candidate_criticals:10.1f} "
            f"{row.avg_reducer_boundary_annotation_zero_skipped_criticals:10.1f} "
            f"{row.avg_reducer_boundary_annotation_zero_skipped_faces:10.1f} "
            f"{row.avg_reducer_boundary_annotation_xors:10.1f} "
            f"{row.avg_reducer_xor_candidate_scans:10.1f} "
            f"{row.avg_reducer_boundary_annotation_total_output_size:10.1f} "
            f"{row.avg_reducer_pivot_eliminations:10.1f} "
            f"{row.avg_reducer_xor_total_input_size:10.1f} "
            f"{row.avg_reducer_xor_total_output_size:10.1f} "
            f"{row.avg_reducer_xor_max_input_size:10.1f} "
            f"{row.avg_reducer_xor_max_output_size:10.1f} "
            f"{row.avg_reducer_xor_inserted_labels:10.1f} "
            f"{row.avg_reducer_xor_removed_labels:10.1f}\n"
        )


def write_csv(rows: list[BenchmarkRow], output: TextIO) -> None:
    writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0]).keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(asdict(row))


def write_json(rows: list[BenchmarkRow], output: TextIO) -> None:
    json.dump([asdict(row) for row in rows], output, indent=2)
    output.write("\n")


def write_rows(rows: list[BenchmarkRow], *, output_format: str, output_path: Path | None) -> None:
    if not rows:
        raise ValueError("No benchmark rows were generated.")

    if output_path is None:
        output = sys.stdout
        close_output = False
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output = output_path.open("w", newline="")
        close_output = True

    try:
        if output_format == "table":
            write_table(rows, output)
        elif output_format == "summary":
            write_summary(rows, output)
        elif output_format == "csv":
            write_csv(rows, output)
        elif output_format == "json":
            write_json(rows, output)
        else:
            raise ValueError(f"Unknown output format: {output_format}")
    finally:
        if close_output:
            output.close()


def write_profile_vs_measured_table(rows: list[ProfileVsMeasuredRow], output: TextIO) -> None:
    headers = [
        "family",
        "seed",
        "size",
        "metric",
        "eff_metric",
        "simplices",
        "cands",
        "prof_cands",
        "profile",
        "measured",
        "match",
        "rank",
        "penalty%",
        "profile_us",
        "best_us",
        "prof_score",
        "best_score",
        "sel_prof_us",
        "best_prof_us",
        "sel_seq_us",
        "best_seq_us",
        "sel_work",
        "best_work",
        "cam_s",
        "pers_s",
    ]
    output.write(" ".join(f"{header:>12}" for header in headers) + "\n")
    for row in rows:
        output.write(
            f"{row.family:>12} "
            f"{row.seed:12d} "
            f"{row.size:12d} "
            f"{row.profile_selection_metric:>12} "
            f"{row.effective_profile_selection_metric:>12} "
            f"{row.num_simplices:12d} "
            f"{row.candidate_count:12d} "
            f"{row.profile_candidate_count:12d} "
            f"{row.profile_selected_algorithm:>12} "
            f"{row.measured_best_algorithm:>12} "
            f"{str(row.profile_matches_measured):>12} "
            f"{row.profile_selected_measured_rank:12d} "
            f"{row.profile_penalty_percent:12.3f} "
            f"{1.0e6 * row.profile_selected_morse_seconds:12.2f} "
            f"{1.0e6 * row.measured_best_morse_seconds:12.2f} "
            f"{row.profile_selected_profile_score:12.6g} "
            f"{_format_optional_float(row.measured_best_profile_score):>12} "
            f"{1.0e6 * row.profile_selected_profile_seconds:12.2f} "
            f"{_format_summary_optional_us(row.measured_best_profile_seconds):>12} "
            f"{1.0e6 * row.profile_selected_sequence_seconds:12.2f} "
            f"{1.0e6 * row.measured_best_sequence_seconds:12.2f} "
            f"{row.profile_selected_reducer_work:12.0f} "
            f"{row.measured_best_reducer_work:12.0f} "
            f"{_format_optional_float(row.gudhi_cam_seconds):>12} "
            f"{_format_optional_float(row.perseus_seconds):>12}\n"
        )


def write_profile_vs_measured_summary(rows: list[ProfileVsMeasuredRow], output: TextIO) -> None:
    grouped: dict[tuple[str, str, int], list[ProfileVsMeasuredRow]] = {}
    for row in rows:
        grouped.setdefault((row.profile_selection_metric, row.family, row.size), []).append(row)

    headers = [
        "metric",
        "eff_metric",
        "family",
        "size",
        "cases",
        "prof_cands",
        "match%",
        "avg_pen%",
        "max_pen%",
        "avg_rank",
        "profile_us",
        "best_us",
        "prof_over%",
        "seq_gap%",
        "work_gap%",
        "prof_sec",
        "meas_sec",
    ]
    output.write(" ".join(f"{header:>12}" for header in headers) + "\n")
    for (metric, family, size), group in sorted(grouped.items()):
        cases = len(group)
        effective_metric = "/".join(
            sorted({row.effective_profile_selection_metric for row in group})
        )
        matches = sum(1 for row in group if row.profile_matches_measured)
        finite_penalties = [
            row.profile_penalty_percent
            for row in group
            if not math.isinf(row.profile_penalty_percent)
        ]
        avg_penalty = (
            sum(finite_penalties) / len(finite_penalties)
            if finite_penalties
            else math.inf
        )
        max_penalty = max((row.profile_penalty_percent for row in group), default=0.0)
        profile_overhead_percent = (
            100.0
            * sum(row.total_profile_seconds for row in group)
            / sum(row.total_measured_morse_seconds for row in group)
            if sum(row.total_measured_morse_seconds for row in group) != 0.0
            else math.inf
        )
        sequence_gap_percent = (
            100.0
            * sum(
                row.profile_selected_sequence_seconds - row.measured_best_sequence_seconds
                for row in group
            )
            / sum(row.measured_best_morse_seconds for row in group)
            if sum(row.measured_best_morse_seconds for row in group) != 0.0
            else math.inf
        )
        work_gap_percent = (
            100.0
            * sum(row.profile_selected_reducer_work - row.measured_best_reducer_work for row in group)
            / sum(row.measured_best_reducer_work for row in group)
            if sum(row.measured_best_reducer_work for row in group) != 0.0
            else math.inf
        )
        output.write(
            f"{metric:>12} "
            f"{effective_metric:>12} "
            f"{family:>12} "
            f"{size:12d} "
            f"{cases:12d} "
            f"{sum(row.profile_candidate_count for row in group) / cases:12.2f} "
            f"{100.0 * matches / cases:12.1f} "
            f"{avg_penalty:12.3f} "
            f"{max_penalty:12.3f} "
            f"{sum(row.profile_selected_measured_rank for row in group) / cases:12.2f} "
            f"{1.0e6 * sum(row.profile_selected_morse_seconds for row in group) / cases:12.2f} "
            f"{1.0e6 * sum(row.measured_best_morse_seconds for row in group) / cases:12.2f} "
            f"{profile_overhead_percent:12.3f} "
            f"{sequence_gap_percent:12.3f} "
            f"{work_gap_percent:12.3f} "
            f"{sum(row.total_profile_seconds for row in group) / cases:12.6g} "
            f"{sum(row.total_measured_morse_seconds for row in group) / cases:12.6g}\n"
        )


def write_profile_vs_measured_csv(
    rows: list[ProfileVsMeasuredRow],
    output: TextIO,
) -> None:
    writer = csv.DictWriter(output, fieldnames=list(asdict(rows[0]).keys()))
    writer.writeheader()
    for row in rows:
        writer.writerow(asdict(row))


def write_profile_vs_measured_json(
    rows: list[ProfileVsMeasuredRow],
    output: TextIO,
) -> None:
    json.dump([asdict(row) for row in rows], output, indent=2)
    output.write("\n")


def write_profile_vs_measured_rows(
    rows: list[ProfileVsMeasuredRow],
    *,
    output_format: str,
    output_path: Path | None,
) -> None:
    if not rows:
        raise ValueError("No profile-vs-measured rows were generated.")

    if output_path is None:
        output = sys.stdout
        close_output = False
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output = output_path.open("w", newline="")
        close_output = True

    try:
        if output_format == "table":
            write_profile_vs_measured_table(rows, output)
        elif output_format == "summary":
            write_profile_vs_measured_summary(rows, output)
        elif output_format == "csv":
            write_profile_vs_measured_csv(rows, output)
        elif output_format == "json":
            write_profile_vs_measured_json(rows, output)
        else:
            raise ValueError(f"Unknown output format: {output_format}")
    finally:
        if close_output:
            output.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=["benchmark", "profile-vs-measured"],
        default="benchmark",
        help=(
            "Run the ordinary benchmark table, or evaluate cheap profile selection "
            "against measured portfolio timings."
        ),
    )
    parser.add_argument("--preset", choices=sorted(PRESETS))
    parser.add_argument(
        "--families",
        nargs="+",
        choices=list(FAMILIES),
    )
    parser.add_argument("--sizes", nargs="+", type=int)
    parser.add_argument("--seeds", nargs="+", type=int)
    parser.add_argument("--repeats", type=int)
    parser.add_argument(
        "--profile-repeats",
        type=int,
        default=1,
        help="Number of cheap profile repetitions for --mode profile-vs-measured.",
    )
    parser.add_argument("--plateau-levels", type=int)
    parser.add_argument(
        "--sequence-algorithm",
        default=mp.DEFAULT_MORSE_SEQUENCE_ALGORITHM,
        help=(
            "Morse sequence constructor to benchmark. "
            f"Implemented: {', '.join(mp.MORSE_SEQUENCE_ALGORITHMS)}. "
            "Use 'portfolio' or 'all' to run the default strategy portfolio."
        ),
    )
    parser.add_argument(
        "--frame-mode",
        choices=list(mp.MORSE_FRAME_MODES),
        default=mp.DEFAULT_MORSE_FRAME_MODE,
        help="Build the Morse sequence/reference frame in one pass or as separate passes.",
    )
    parser.add_argument(
        "--validation-mode",
        choices=["core", "materialized"],
        default="core",
        help=(
            "Use compact C++ barcode signatures for validation, or materialize all intervals "
            "in Python for debugging."
        ),
    )
    parser.add_argument("--roadmap-cache", type=Path, default=DEFAULT_ROADMAP_CACHE)
    parser.add_argument(
        "--download-roadmap-data",
        action="store_true",
        help="Download missing PH-roadmap benchmark datasets into --roadmap-cache.",
    )
    parser.add_argument(
        "--time-gudhi-cam",
        action="store_true",
        help="Also time GUDHI persistence as a current CAM implementation comparison.",
    )
    parser.add_argument(
        "--time-perseus",
        action="store_true",
        help="Also time Perseus on an exported nmfsimtop complex, if the executable is available.",
    )
    parser.add_argument(
        "--perseus-executable",
        default="perseus",
        help="Path or command name for the Perseus executable used with --time-perseus.",
    )
    parser.add_argument(
        "--profile-selection-metric",
        default="estimated_reducer_work",
        help=(
            "Cheap profile metric for --mode profile-vs-measured. "
            f"Supported: {', '.join(BENCHMARK_PROFILE_SELECTION_METRICS)}."
        ),
    )
    parser.add_argument(
        "--profile-selection-metrics",
        nargs="+",
        help=(
            "Evaluate several cheap profile metrics in one --mode profile-vs-measured "
            "run, reusing the same profiles and measured portfolio timings for each "
            "complex."
        ),
    )
    parser.add_argument(
        "--profile-candidate-gate",
        choices=list(PROFILE_CANDIDATE_GATES),
        default="all",
        help=(
            "Candidate subset profiled by --mode profile-vs-measured. "
            "'all' profiles every measured strategy; 'family-aware' profiles a small "
            "experimental subset chosen from the benchmark family; 'tabulated' reads "
            "candidate sets from --profile-candidate-table."
        ),
    )
    parser.add_argument(
        "--profile-candidate-table",
        type=Path,
        help=(
            "CSV produced by calibrate_profile_gate.py for "
            "--profile-candidate-gate tabulated."
        ),
    )
    parser.add_argument(
        "--measured-selection-metric",
        default="morse_seconds",
        help=(
            "Measured best-strategy metric for --mode profile-vs-measured. "
            f"Supported: {', '.join(mp.MORSE_ALGORITHM_SELECTION_METRICS)}."
        ),
    )
    parser.add_argument("--format", choices=["table", "summary", "csv", "json"], default="table")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    preset = PRESETS[args.preset] if args.preset is not None else DEFAULT_PRESET
    families = args.families if args.families is not None else preset.families
    sizes = args.sizes if args.sizes is not None else preset.sizes
    seeds = args.seeds if args.seeds is not None else preset.seeds
    repeats = args.repeats if args.repeats is not None else preset.repeats
    plateau_levels = (
        args.plateau_levels if args.plateau_levels is not None else preset.plateau_levels
    )
    if args.mode == "profile-vs-measured":
        profile_rows = run_profile_vs_measured(
            families=families,
            sizes=sizes,
            seeds=seeds,
            repeats=repeats,
            plateau_levels=plateau_levels,
            profile_repeats=args.profile_repeats,
            time_gudhi=args.time_gudhi_cam,
            time_perseus=args.time_perseus,
            perseus_executable=args.perseus_executable,
            sequence_algorithm=args.sequence_algorithm,
            frame_mode=args.frame_mode,
            validation_mode=args.validation_mode,
            profile_selection_metric=(
                args.profile_selection_metrics
                if args.profile_selection_metrics is not None
                else args.profile_selection_metric
            ),
            profile_candidate_gate=args.profile_candidate_gate,
            profile_candidate_table=args.profile_candidate_table,
            measured_selection_metric=args.measured_selection_metric,
            roadmap_cache=args.roadmap_cache,
            download_roadmap_data=args.download_roadmap_data,
        )
        write_profile_vs_measured_rows(
            profile_rows,
            output_format=args.format,
            output_path=args.output,
        )
    else:
        rows = run_benchmarks(
            families=families,
            sizes=sizes,
            seeds=seeds,
            repeats=repeats,
            plateau_levels=plateau_levels,
            time_gudhi=args.time_gudhi_cam,
            time_perseus=args.time_perseus,
            perseus_executable=args.perseus_executable,
            sequence_algorithm=args.sequence_algorithm,
            frame_mode=args.frame_mode,
            validation_mode=args.validation_mode,
            roadmap_cache=args.roadmap_cache,
            download_roadmap_data=args.download_roadmap_data,
        )
        write_rows(rows, output_format=args.format, output_path=args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
