"""Plotly helpers for inspecting low-dimensional Morse persistence examples."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Iterable, Sequence

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from . import (
    FilteredComplex,
    MorseSequence,
    PersistenceDiagram,
    Simplex,
    compute_coreference_map,
    compute_morse_complex,
    compute_morse_persistence,
    compute_morse_sequence_and_reference_map,
)


Point2D = tuple[float, float]
Point3D = tuple[float, float, float]
SimplexLink = tuple[Simplex, Simplex]
SimplexSupports = dict[Simplex, tuple[Simplex, ...]]


@dataclass(frozen=True)
class TriangulatedScalarField:
    complex: FilteredComplex
    coordinates: dict[int, Point2D]
    values: dict[int, float]
    triangles: tuple[Simplex, ...]
    original_values: dict[int, float] | None = None
    expected_minima: tuple[int, ...] = ()
    expected_maxima: tuple[int, ...] = ()


@dataclass(frozen=True)
class TraceHeightBinding:
    trace_index: int
    trace_type: str
    noisy_z: tuple[float | None, ...]
    original_z: tuple[float | None, ...]
    noisy_intensity: tuple[float, ...] | None = None
    original_intensity: tuple[float, ...] | None = None


def build_noisy_sine_square(
    *,
    size: int = 18,
    noise: float = 0.035,
    seed: int = 7,
) -> TriangulatedScalarField:
    """Build a lower-star triangulated square from sin(x)sin(y) plus noise.

    The domain is ``[0, 2*pi]^2``. The noiseless function has two maxima and
    two minima in the interior of this domain.
    """

    if size < 2:
        raise ValueError("size must be at least 2.")

    rng = random.Random(seed)

    def vertex_id(i: int, j: int) -> int:
        return i * (size + 1) + j

    def nearest_vertex(x: float, y: float) -> int:
        i = min(size, max(0, int(round((x / (2.0 * math.pi)) * size))))
        j = min(size, max(0, int(round((y / (2.0 * math.pi)) * size))))
        return vertex_id(i, j)

    coordinates: dict[int, Point2D] = {}
    values: dict[int, float] = {}
    original_values: dict[int, float] = {}
    for i in range(size + 1):
        x = 2.0 * math.pi * float(i) / float(size)
        for j in range(size + 1):
            y = 2.0 * math.pi * float(j) / float(size)
            vertex = vertex_id(i, j)
            original_value = math.sin(x) * math.sin(y)
            coordinates[vertex] = (x, y)
            original_values[vertex] = original_value
            values[vertex] = original_value + noise * rng.gauss(0.0, 1.0)

    triangles: list[Simplex] = []
    for i in range(size):
        for j in range(size):
            v00 = vertex_id(i, j)
            v10 = vertex_id(i + 1, j)
            v01 = vertex_id(i, j + 1)
            v11 = vertex_id(i + 1, j + 1)
            if (i + j) % 2 == 0:
                triangles.append((v00, v10, v11))
                triangles.append((v00, v11, v01))
            else:
                triangles.append((v00, v10, v01))
                triangles.append((v10, v11, v01))

    def lower_star(simplex: Sequence[int]) -> float:
        return max(values[vertex] for vertex in simplex)

    complex_ = FilteredComplex.from_facets(triangles, simplex_filtration=lower_star)
    return TriangulatedScalarField(
        complex=complex_,
        coordinates=coordinates,
        values=values,
        triangles=tuple(triangles),
        original_values=original_values,
        expected_minima=(
            nearest_vertex(0.5 * math.pi, 1.5 * math.pi),
            nearest_vertex(1.5 * math.pi, 0.5 * math.pi),
        ),
        expected_maxima=(
            nearest_vertex(0.5 * math.pi, 0.5 * math.pi),
            nearest_vertex(1.5 * math.pi, 1.5 * math.pi),
        ),
    )


def simplex_center(
    simplex: Sequence[int],
    coordinates: dict[int, Point2D],
    values: dict[int, float],
    *,
    lift: float = 0.0,
) -> Point3D:
    xs = [coordinates[vertex][0] for vertex in simplex]
    ys = [coordinates[vertex][1] for vertex in simplex]
    zs = [values[vertex] for vertex in simplex]
    return (
        sum(xs) / len(xs),
        sum(ys) / len(ys),
        sum(zs) / len(zs) + lift,
    )


def simplex_centers(
    simplexes: Iterable[Simplex],
    coordinates: dict[int, Point2D],
    values: dict[int, float],
    *,
    lift: float = 0.0,
) -> tuple[list[float], list[float], list[float]]:
    x: list[float] = []
    y: list[float] = []
    z: list[float] = []
    for simplex in simplexes:
        cx, cy, cz = simplex_center(simplex, coordinates, values, lift=lift)
        x.append(cx)
        y.append(cy)
        z.append(cz)
    return x, y, z


def simplex_links_to_trace_points(
    links: Iterable[SimplexLink],
    coordinates: dict[int, Point2D],
    values: dict[int, float],
    *,
    lift: float = 0.0,
) -> tuple[list[float | None], list[float | None], list[float | None]]:
    x: list[float | None] = []
    y: list[float | None] = []
    z: list[float | None] = []
    for source, target in links:
        sx, sy, sz = simplex_center(source, coordinates, values, lift=lift)
        tx, ty, tz = simplex_center(target, coordinates, values, lift=lift)
        x.extend([sx, tx, None])
        y.extend([sy, ty, None])
        z.extend([sz, tz, None])
    return x, y, z


def edge_trace_points(
    edges: Iterable[Simplex],
    coordinates: dict[int, Point2D],
    values: dict[int, float],
    *,
    lift: float = 0.0,
) -> tuple[list[float | None], list[float | None], list[float | None]]:
    x: list[float | None] = []
    y: list[float | None] = []
    z: list[float | None] = []
    for edge in edges:
        if len(edge) != 2:
            continue
        v0, v1 = edge
        x0, y0 = coordinates[v0]
        x1, y1 = coordinates[v1]
        z0 = values[v0] + lift
        z1 = values[v1] + lift
        x.extend([x0, x1, None])
        y.extend([y0, y1, None])
        z.extend([z0, z1, None])
    return x, y, z


def edge_group_trace_points(
    edge_groups: Iterable[Iterable[Simplex]],
    coordinates: dict[int, Point2D],
    values: dict[int, float],
    *,
    lift: float = 0.0,
) -> tuple[list[float | None], list[float | None], list[float | None]]:
    x: list[float | None] = []
    y: list[float | None] = []
    z: list[float | None] = []
    for edges in edge_groups:
        x_edges, y_edges, z_edges = edge_trace_points(
            edges,
            coordinates,
            values,
            lift=lift,
        )
        x.extend(x_edges)
        y.extend(y_edges)
        z.extend(z_edges)
    return x, y, z


def triangulation_edges(complex_: FilteredComplex) -> tuple[Simplex, ...]:
    return tuple(record.vertices for record in complex_.simplices() if record.dimension == 1)


def _xor_sorted(lhs: Sequence[int], rhs: Sequence[int]) -> tuple[int, ...]:
    labels = set(lhs)
    for label in rhs:
        if label in labels:
            labels.remove(label)
        else:
            labels.add(label)
    return tuple(sorted(labels))


def _critical_ids_for_support(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    critical_simplexes: Iterable[Sequence[int]] | None,
    *,
    dimension: int,
) -> tuple[int, ...]:
    if critical_simplexes is None:
        return tuple(
            critical_id
            for critical_id, simplex_id in enumerate(sequence.critical_simplices)
            if complex_.dimension(simplex_id) == dimension
        )

    critical_ids: list[int] = []
    for simplex in critical_simplexes:
        simplex_id = complex_.require_simplex_id(simplex)
        if complex_.dimension(simplex_id) != dimension:
            continue
        critical_id = sequence.critical_index(simplex_id)
        if critical_id < 0:
            raise ValueError(f"Simplex {complex_.vertices(simplex_id)} is not critical.")
        critical_ids.append(critical_id)
    return tuple(dict.fromkeys(critical_ids))


def _annotation_supports_by_critical(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    annotations: Sequence[Sequence[int]],
    critical_simplexes: Iterable[Sequence[int]] | None,
    *,
    dimension: int,
) -> SimplexSupports:
    """Collect same-dimensional simplex supports for selected critical cells."""

    critical_ids = _critical_ids_for_support(
        complex_,
        sequence,
        critical_simplexes,
        dimension=dimension,
    )
    selected = set(critical_ids)
    supports: dict[int, list[Simplex]] = {critical_id: [] for critical_id in critical_ids}
    for simplex_id, annotation in enumerate(annotations):
        if complex_.dimension(simplex_id) != dimension:
            continue
        simplex = complex_.vertices(simplex_id)
        for critical_id in annotation:
            if critical_id in selected:
                supports[critical_id].append(simplex)

    return {
        complex_.vertices(sequence.critical_simplices[critical_id]): tuple(supports[critical_id])
        for critical_id in critical_ids
    }


def reference_cocycle_edges(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    references: Sequence[Sequence[int]],
    critical_simplexes: Iterable[Sequence[int]] | None = None,
    *,
    dimension: int = 1,
) -> SimplexSupports:
    """Return reference cocycle supports as original same-dimensional simplexes.

    For the surface demo with ``dimension=1``, this matches the private
    ``referenceMap.cocycles(crit1)`` convention: each critical edge is mapped
    to the mesh edges whose reference annotation contains that critical edge.
    """

    return _annotation_supports_by_critical(
        complex_,
        sequence,
        references,
        critical_simplexes,
        dimension=dimension,
    )


def coreference_cycle_edges(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    coreferences: Sequence[Sequence[int]],
    critical_simplexes: Iterable[Sequence[int]] | None = None,
    *,
    dimension: int = 1,
) -> SimplexSupports:
    """Return coreference cycle supports as original same-dimensional simplexes.

    For the surface demo with ``dimension=1``, this matches the private
    ``coreferenceMap.cycles(crit1)`` convention.
    """

    return _annotation_supports_by_critical(
        complex_,
        sequence,
        coreferences,
        critical_simplexes,
        dimension=dimension,
    )


def morse_boundary_links(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    references: Sequence[Sequence[int]],
) -> tuple[SimplexLink, ...]:
    """Return links from critical 2-simplices to critical 1-simplices."""

    morse_complex = compute_morse_complex(complex_, sequence, references)
    links: list[SimplexLink] = []
    for critical_id, simplex_id in enumerate(sequence.critical_simplices):
        if complex_.dimension(simplex_id) != 2:
            continue
        triangle = complex_.vertices(simplex_id)
        for boundary_critical_id in morse_complex.boundary(critical_id):
            boundary_simplex_id = sequence.critical_simplices[boundary_critical_id]
            if complex_.dimension(boundary_simplex_id) == 1:
                links.append((complex_.vertices(boundary_simplex_id), triangle))
    return tuple(links)


def _persistent_critical_simplexes_by_dimension(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    diagram: PersistenceDiagram,
    *,
    threshold: float,
    dimension: int,
) -> tuple[Simplex, ...]:
    selected_ids: set[int] = set()
    for pair in diagram.finite_pairs:
        persistence = pair.death_value - pair.birth_value
        if persistence < threshold - 1e-12:
            continue
        if complex_.dimension(pair.birth) == dimension:
            selected_ids.add(pair.birth)
        if complex_.dimension(pair.death) == dimension:
            selected_ids.add(pair.death)
    for interval in diagram.essential:
        if interval.dimension == dimension:
            selected_ids.add(interval.birth)

    return tuple(
        complex_.vertices(simplex_id)
        for simplex_id in sequence.critical_simplices
        if simplex_id in selected_ids
    )


def morse_coboundary_links(
    complex_: FilteredComplex,
    sequence: MorseSequence,
    coreferences: Sequence[Sequence[int]],
) -> tuple[SimplexLink, ...]:
    """Return links from critical 1-simplices to critical 2-simplices."""

    links: list[SimplexLink] = []
    for critical_id, simplex_id in enumerate(sequence.critical_simplices):
        if complex_.dimension(simplex_id) != 1:
            continue
        edge = complex_.vertices(simplex_id)
        coboundary_annotation: tuple[int, ...] = ()
        for coface in complex_.coboundary(simplex_id):
            if complex_.dimension(coface) == 2:
                coboundary_annotation = _xor_sorted(
                    coboundary_annotation,
                    coreferences[coface],
                )
        for coboundary_critical_id in coboundary_annotation:
            coboundary_simplex_id = sequence.critical_simplices[coboundary_critical_id]
            if complex_.dimension(coboundary_simplex_id) == 2:
                links.append((edge, complex_.vertices(coboundary_simplex_id)))
    return tuple(links)


def thresholds_from_diagram(
    diagram: PersistenceDiagram,
    *,
    steps: int = 12,
) -> tuple[float, ...]:
    persistences = [pair.death_value - pair.birth_value for pair in diagram.finite_pairs]
    max_persistence = max(persistences, default=0.0)
    if max_persistence <= 0.0:
        return (0.0,)
    count = max(2, steps)
    return tuple(max_persistence * float(index) / float(count - 1) for index in range(count))


def _sparse_slider_label(index: int, count: int, value: float) -> str:
    if count <= 6:
        return f"{value:.3g}"
    stride = max(1, math.ceil((count - 1) / 4))
    if index in {0, count - 1} or index % stride == 0:
        return f"{value:.3g}"
    return ""


def _surface_vertex_values(
    field: TriangulatedScalarField,
    values: dict[int, float],
) -> tuple[float, ...]:
    return tuple(values[vertex] for vertex in sorted(field.coordinates))


def _surface_trace(
    field: TriangulatedScalarField,
    *,
    opacity: float,
    values: dict[int, float] | None = None,
) -> go.Mesh3d:
    surface_values = values or field.values
    vertices = sorted(field.coordinates)
    vertex_to_mesh = {vertex: index for index, vertex in enumerate(vertices)}
    return go.Mesh3d(
        x=[field.coordinates[vertex][0] for vertex in vertices],
        y=[field.coordinates[vertex][1] for vertex in vertices],
        z=[surface_values[vertex] for vertex in vertices],
        i=[vertex_to_mesh[triangle[0]] for triangle in field.triangles],
        j=[vertex_to_mesh[triangle[1]] for triangle in field.triangles],
        k=[vertex_to_mesh[triangle[2]] for triangle in field.triangles],
        intensity=[surface_values[vertex] for vertex in vertices],
        colorscale="RdBu",
        reversescale=True,
        opacity=opacity,
        showscale=True,
        colorbar=dict(title="f"),
        name="triangulated surface",
        hovertemplate="x=%{x:.3f}<br>y=%{y:.3f}<br>f=%{z:.3f}<extra></extra>",
    )


def _critical_triangle_trace(
    simplexes: Sequence[Simplex],
    field: TriangulatedScalarField,
    *,
    lift: float,
    values: dict[int, float] | None = None,
) -> go.Mesh3d:
    height_values = values or field.values
    triangles = [simplex for simplex in simplexes if len(simplex) == 3]
    x: list[float] = []
    y: list[float] = []
    z: list[float] = []
    i: list[int] = []
    j: list[int] = []
    k: list[int] = []
    for triangle in triangles:
        offset = len(x)
        for vertex in triangle:
            vx, vy = field.coordinates[vertex]
            x.append(vx)
            y.append(vy)
            z.append(height_values[vertex] + lift)
        i.append(offset)
        j.append(offset + 1)
        k.append(offset + 2)
    return go.Mesh3d(
        x=x,
        y=y,
        z=z,
        i=i,
        j=j,
        k=k,
        color="#9467bd",
        opacity=0.55,
        name="critical 2",
        hoverinfo="skip",
        showlegend=True,
    )


def _add_critical_traces(
    fig: go.Figure,
    field: TriangulatedScalarField,
    critical_simplexes: Sequence[Simplex],
    *,
    lift: float,
    values: dict[int, float] | None = None,
    original_values: dict[int, float] | None = None,
) -> tuple[TraceHeightBinding, ...]:
    height_values = values or field.values
    bindings: list[TraceHeightBinding] = []
    critical_vertices = [simplex for simplex in critical_simplexes if len(simplex) == 1]
    critical_edges = [simplex for simplex in critical_simplexes if len(simplex) == 2]
    critical_triangles = [simplex for simplex in critical_simplexes if len(simplex) == 3]

    x0, y0, z0 = simplex_centers(critical_vertices, field.coordinates, height_values, lift=lift)
    fig.add_trace(
        go.Scatter3d(
            x=x0,
            y=y0,
            z=z0,
            mode="markers",
            marker=dict(size=6, color="#1f77b4", symbol="circle"),
            text=[str(simplex) for simplex in critical_vertices],
            hovertemplate="critical vertex %{text}<extra></extra>",
            name="critical 0",
        ),
        row=1,
        col=1,
    )
    if original_values is not None:
        _, _, original_z0 = simplex_centers(
            critical_vertices,
            field.coordinates,
            original_values,
            lift=lift,
        )
        bindings.append(
            TraceHeightBinding(
                trace_index=len(fig.data) - 1,
                trace_type="scatter3d",
                noisy_z=tuple(z0),
                original_z=tuple(original_z0),
            )
        )

    x1, y1, z1 = edge_trace_points(
        critical_edges,
        field.coordinates,
        height_values,
        lift=lift,
    )
    fig.add_trace(
        go.Scatter3d(
            x=x1,
            y=y1,
            z=z1,
            mode="lines",
            line=dict(width=7, color="#ff7f0e"),
            hoverinfo="skip",
            name="critical 1",
        ),
        row=1,
        col=1,
    )
    if original_values is not None:
        _, _, original_z1 = edge_trace_points(
            critical_edges,
            field.coordinates,
            original_values,
            lift=lift,
        )
        bindings.append(
            TraceHeightBinding(
                trace_index=len(fig.data) - 1,
                trace_type="scatter3d",
                noisy_z=tuple(z1),
                original_z=tuple(original_z1),
            )
        )

    if critical_triangles:
        fig.add_trace(
            _critical_triangle_trace(
                critical_triangles,
                field,
                lift=lift,
                values=height_values,
            ),
            row=1,
            col=1,
        )
        triangle_mesh_index = len(fig.data) - 1
        if original_values is not None:
            original_triangle_z_list: list[float | None] = []
            for triangle in critical_triangles:
                for vertex in triangle:
                    original_triangle_z_list.append(original_values[vertex] + lift)
            noisy_triangle_z_list: list[float | None] = []
            for triangle in critical_triangles:
                for vertex in triangle:
                    noisy_triangle_z_list.append(height_values[vertex] + lift)
            bindings.append(
                TraceHeightBinding(
                    trace_index=triangle_mesh_index,
                    trace_type="mesh3d",
                    noisy_z=tuple(noisy_triangle_z_list),
                    original_z=tuple(original_triangle_z_list),
                )
            )
        x2, y2, z2 = simplex_centers(
            critical_triangles,
            field.coordinates,
            height_values,
            lift=2.0 * lift,
        )
        fig.add_trace(
            go.Scatter3d(
                x=x2,
                y=y2,
                z=z2,
                mode="markers",
                marker=dict(size=7, color="#9467bd", symbol="diamond"),
                text=[str(simplex) for simplex in critical_triangles],
                hovertemplate="critical triangle %{text}<extra></extra>",
                name="critical 2 centers",
                showlegend=False,
            ),
            row=1,
            col=1,
        )
        if original_values is not None:
            _, _, original_z2 = simplex_centers(
                critical_triangles,
                field.coordinates,
                original_values,
                lift=2.0 * lift,
            )
            bindings.append(
                TraceHeightBinding(
                    trace_index=len(fig.data) - 1,
                    trace_type="scatter3d",
                    noisy_z=tuple(z2),
                    original_z=tuple(original_z2),
                )
            )
    return tuple(bindings)


def _add_morse_link_traces(
    fig: go.Figure,
    field: TriangulatedScalarField,
    boundary_links: Sequence[SimplexLink],
    coboundary_links: Sequence[SimplexLink],
    *,
    lift: float,
) -> None:
    xb, yb, zb = simplex_links_to_trace_points(
        boundary_links,
        field.coordinates,
        field.values,
        lift=2.5 * lift,
    )
    fig.add_trace(
        go.Scatter3d(
            x=xb,
            y=yb,
            z=zb,
            mode="lines",
            line=dict(width=5, color="rgba(20,20,20,0.85)"),
            hoverinfo="skip",
            name="Morse boundary 2-to-1",
        ),
        row=1,
        col=1,
    )
    xc, yc, zc = simplex_links_to_trace_points(
        coboundary_links,
        field.coordinates,
        field.values,
        lift=3.2 * lift,
    )
    fig.add_trace(
        go.Scatter3d(
            x=xc,
            y=yc,
            z=zc,
            mode="lines",
            line=dict(width=3, color="rgba(23,190,207,0.9)"),
            hoverinfo="skip",
            name="Morse coboundary 1-to-2",
        ),
        row=1,
        col=1,
    )


def _support_groups_for_selected(
    supports: SimplexSupports,
    selected_critical_simplexes: Sequence[Simplex],
) -> tuple[tuple[Simplex, ...], ...]:
    return tuple(
        supports[critical_simplex]
        for critical_simplex in selected_critical_simplexes
        if critical_simplex in supports
    )


def _add_cycle_cocycle_threshold_traces(
    fig: go.Figure,
    field: TriangulatedScalarField,
    sequence: MorseSequence,
    diagram: PersistenceDiagram,
    thresholds: Sequence[float],
    cycle_supports: SimplexSupports,
    cocycle_supports: SimplexSupports,
    *,
    lift: float,
    show_coreference_cycles: bool,
    values: dict[int, float] | None = None,
    original_values: dict[int, float] | None = None,
) -> tuple[tuple[tuple[int, ...], ...], tuple[TraceHeightBinding, ...]]:
    height_values = values or field.values
    trace_groups: list[tuple[int, ...]] = []
    bindings: list[TraceHeightBinding] = []
    for threshold_index, threshold in enumerate(thresholds):
        selected = _persistent_critical_simplexes_by_dimension(
            field.complex,
            sequence,
            diagram,
            threshold=threshold,
            dimension=1,
        )
        visible = threshold_index == 0
        threshold_trace_indices: list[int] = []

        if show_coreference_cycles:
            cycle_groups = _support_groups_for_selected(cycle_supports, selected)
            x_cycle, y_cycle, z_cycle = edge_group_trace_points(
                cycle_groups,
                field.coordinates,
                height_values,
                lift=2.6 * lift,
            )
            fig.add_trace(
                go.Scatter3d(
                    x=x_cycle,
                    y=y_cycle,
                    z=z_cycle,
                    mode="lines",
                    line=dict(width=5, color="rgba(255,127,14,0.92)"),
                    hoverinfo="skip",
                    visible=visible,
                    name="1-cycles",
                    showlegend=threshold_index == 0,
                ),
                row=1,
                col=1,
            )
            threshold_trace_indices.append(len(fig.data) - 1)
            if original_values is not None:
                _, _, original_z_cycle = edge_group_trace_points(
                    cycle_groups,
                    field.coordinates,
                    original_values,
                    lift=2.6 * lift,
                )
                bindings.append(
                    TraceHeightBinding(
                        trace_index=len(fig.data) - 1,
                        trace_type="scatter3d",
                        noisy_z=tuple(z_cycle),
                        original_z=tuple(original_z_cycle),
                    )
                )

        cocycle_groups = _support_groups_for_selected(cocycle_supports, selected)
        x_cocycle, y_cocycle, z_cocycle = edge_group_trace_points(
            cocycle_groups,
            field.coordinates,
            height_values,
            lift=3.2 * lift,
        )
        fig.add_trace(
            go.Scatter3d(
                x=x_cocycle,
                y=y_cocycle,
                z=z_cocycle,
                mode="lines",
                line=dict(width=4, color="rgba(241,196,15,0.98)"),
                hoverinfo="skip",
                visible=visible,
                name="1-cocycles",
                showlegend=threshold_index == 0,
            ),
            row=1,
            col=1,
        )
        threshold_trace_indices.append(len(fig.data) - 1)
        if original_values is not None:
            _, _, original_z_cocycle = edge_group_trace_points(
                cocycle_groups,
                field.coordinates,
                original_values,
                lift=3.2 * lift,
            )
            bindings.append(
                TraceHeightBinding(
                    trace_index=len(fig.data) - 1,
                    trace_type="scatter3d",
                    noisy_z=tuple(z_cocycle),
                    original_z=tuple(original_z_cocycle),
                )
            )
        trace_groups.append(tuple(threshold_trace_indices))
    return tuple(trace_groups), tuple(bindings)


def _interval_payload(
    diagram: PersistenceDiagram,
    complex_: FilteredComplex,
) -> list[dict[str, object]]:
    payload: list[dict[str, object]] = []
    for pair in diagram.finite_pairs:
        payload.append(
            {
                "pair": pair,
                "persistence": pair.death_value - pair.birth_value,
                "birth_simplex": complex_.vertices(pair.birth),
                "death_simplex": complex_.vertices(pair.death),
            }
        )
    return payload


def _add_persistence_traces(
    fig: go.Figure,
    field: TriangulatedScalarField,
    diagram: PersistenceDiagram,
    thresholds: Sequence[float],
    *,
    lift: float,
    values: dict[int, float] | None = None,
    original_values: dict[int, float] | None = None,
) -> tuple[list[int], list[int], tuple[TraceHeightBinding, ...]]:
    height_values = values or field.values
    payload = _interval_payload(diagram, field.complex)
    axis_values = [
        abs(value)
        for value in height_values.values()
    ] + [
        abs(value)
        for value in (original_values or {}).values()
    ] + [
        abs(pair.birth_value)
        for pair in diagram.finite_pairs
    ] + [
        abs(pair.death_value)
        for pair in diagram.finite_pairs
    ] + [1.0]
    axis_radius = max(axis_values) + 0.1
    fig.add_trace(
        go.Scatter(
            x=[-axis_radius, axis_radius],
            y=[-axis_radius, axis_radius],
            mode="lines",
            line=dict(color="rgba(80,80,80,0.45)", width=1, dash="dash"),
            hoverinfo="skip",
            showlegend=False,
            name="diagonal",
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Scatter(
            x=[item["pair"].birth_value for item in payload],
            y=[item["pair"].death_value for item in payload],
            mode="markers",
            marker=dict(size=7, color="rgba(80,80,80,0.22)"),
            hoverinfo="skip",
            name="all finite intervals",
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    persistence_trace_indices: list[int] = []
    surface_highlight_indices: list[int] = []
    bindings: list[TraceHeightBinding] = []
    for threshold_index, threshold in enumerate(thresholds):
        selected = [
            item
            for item in payload
            if float(item["persistence"]) >= threshold - 1e-12
        ]
        visible = threshold_index == 0
        fig.add_trace(
            go.Scatter(
                x=[item["pair"].birth_value for item in selected],
                y=[item["pair"].death_value for item in selected],
                mode="markers",
                marker=dict(
                    size=11,
                    color=[item["pair"].dimension for item in selected],
                    colorscale="Viridis",
                    cmin=0,
                    cmax=2,
                    line=dict(color="white", width=1),
                ),
                customdata=[
                    [
                        item["pair"].dimension,
                        item["persistence"],
                        item["birth_simplex"],
                        item["death_simplex"],
                    ]
                    for item in selected
                ],
                hovertemplate=(
                    "dim %{customdata[0]}<br>"
                    "birth=%{x:.4f}<br>"
                    "death=%{y:.4f}<br>"
                    "persistence=%{customdata[1]:.4f}<br>"
                    "birth simplex=%{customdata[2]}<br>"
                    "death simplex=%{customdata[3]}<extra></extra>"
                ),
                visible=visible,
                name=f"persistence >= {threshold:.3g}",
                showlegend=False,
            ),
            row=1,
            col=2,
        )
        persistence_trace_indices.append(len(fig.data) - 1)

        birth_simplexes = [item["birth_simplex"] for item in selected]
        xb, yb, zb = simplex_centers(
            birth_simplexes,
            field.coordinates,
            height_values,
            lift=4.0 * lift,
        )
        fig.add_trace(
            go.Scatter3d(
                x=xb,
                y=yb,
                z=zb,
                mode="markers",
                marker=dict(size=7, color="#2ca02c", symbol="diamond"),
                text=[str(simplex) for simplex in birth_simplexes],
                hovertemplate="persistent birth %{text}<extra></extra>",
                visible=visible,
                name="persistent births",
                showlegend=False,
            ),
            row=1,
            col=1,
        )
        surface_highlight_indices.append(len(fig.data) - 1)
        if original_values is not None:
            _, _, original_zb = simplex_centers(
                birth_simplexes,
                field.coordinates,
                original_values,
                lift=4.0 * lift,
            )
            bindings.append(
                TraceHeightBinding(
                    trace_index=len(fig.data) - 1,
                    trace_type="scatter3d",
                    noisy_z=tuple(zb),
                    original_z=tuple(original_zb),
                )
            )

        death_simplexes = [item["death_simplex"] for item in selected]
        xd, yd, zd = simplex_centers(
            death_simplexes,
            field.coordinates,
            height_values,
            lift=4.5 * lift,
        )
        fig.add_trace(
            go.Scatter3d(
                x=xd,
                y=yd,
                z=zd,
                mode="markers",
                marker=dict(size=7, color="#d62728", symbol="x"),
                text=[str(simplex) for simplex in death_simplexes],
                hovertemplate="persistent death %{text}<extra></extra>",
                visible=visible,
                name="persistent deaths",
                showlegend=False,
            ),
            row=1,
            col=1,
        )
        surface_highlight_indices.append(len(fig.data) - 1)
        if original_values is not None:
            _, _, original_zd = simplex_centers(
                death_simplexes,
                field.coordinates,
                original_values,
                lift=4.5 * lift,
            )
            bindings.append(
                TraceHeightBinding(
                    trace_index=len(fig.data) - 1,
                    trace_type="scatter3d",
                    noisy_z=tuple(zd),
                    original_z=tuple(original_zd),
                )
            )

    fig.update_xaxes(range=[-axis_radius, axis_radius], title_text="birth value", row=1, col=2)
    fig.update_yaxes(range=[-axis_radius, axis_radius], title_text="death value", row=1, col=2)
    return persistence_trace_indices, surface_highlight_indices, tuple(bindings)


def _critical_count_text(simplexes: Sequence[Simplex]) -> str:
    counts = [sum(1 for simplex in simplexes if len(simplex) == dim + 1) for dim in range(3)]
    return f"{counts[0]} critical vertices, {counts[1]} critical edges, {counts[2]} critical triangles"


def _frame_trace(binding: TraceHeightBinding, *, use_original: bool) -> go.Mesh3d | go.Scatter3d:
    z = binding.original_z if use_original else binding.noisy_z
    if binding.trace_type == "mesh3d":
        intensity = binding.original_intensity if use_original else binding.noisy_intensity
        kwargs: dict[str, object] = {"z": list(z)}
        if intensity is not None:
            kwargs["intensity"] = list(intensity)
        return go.Mesh3d(**kwargs)
    return go.Scatter3d(z=list(z))


def _surface_mode_frames(bindings: Sequence[TraceHeightBinding]) -> tuple[go.Frame, go.Frame]:
    trace_indices = [binding.trace_index for binding in bindings]
    return (
        go.Frame(
            name="noisy",
            data=[_frame_trace(binding, use_original=False) for binding in bindings],
            traces=trace_indices,
        ),
        go.Frame(
            name="original",
            data=[_frame_trace(binding, use_original=True) for binding in bindings],
            traces=trace_indices,
        ),
    )


def _surface_mode_menu() -> dict[str, object]:
    animation_options = {
        "mode": "immediate",
        "frame": {"duration": 0, "redraw": True},
        "transition": {"duration": 0},
    }
    return {
        "type": "buttons",
        "direction": "right",
        "active": 0,
        "x": 0.01,
        "y": 1.165,
        "xanchor": "left",
        "yanchor": "top",
        "pad": {"r": 6, "t": 0},
        "bgcolor": "rgba(255,255,255,0.86)",
        "bordercolor": "rgba(80,80,80,0.25)",
        "borderwidth": 1,
        "font": {"size": 11},
        "buttons": [
            {
                "label": "Noisy",
                "method": "animate",
                "args": [["noisy"], animation_options],
            },
            {
                "label": "Original",
                "method": "animate",
                "args": [["original"], animation_options],
            },
        ],
    }


def plot_morse_surface_with_persistence(
    field: TriangulatedScalarField,
    *,
    algorithm: str = "f-max",
    surface_opacity: float = 0.42,
    persistence_steps: int = 12,
    show_coboundary: bool = True,
) -> go.Figure:
    """Create a two-panel interactive Plotly figure.

    Left panel: rotatable triangulated surface with critical simplexes and
    reference cocycle and coreference cycle supports. Right panel: persistence
    diagram with a threshold slider. A second slider controls the surface opacity.
    """

    frame = compute_morse_sequence_and_reference_map(field.complex, algorithm=algorithm)
    sequence = frame.sequence
    references = frame.references
    coreferences = compute_coreference_map(field.complex, sequence, algorithm=algorithm)
    diagram = compute_morse_persistence(field.complex, sequence, references)
    critical_simplexes = sequence.critical_simplices_as_simplices(field.complex)
    cocycle_supports = reference_cocycle_edges(field.complex, sequence, references)
    cycle_supports = coreference_cycle_edges(field.complex, sequence, coreferences)
    thresholds = thresholds_from_diagram(diagram, steps=persistence_steps)
    original_values = field.original_values

    z_values = list(field.values.values())
    z_span = max(z_values) - min(z_values) if z_values else 1.0
    lift = 0.025 * (z_span if z_span > 0.0 else 1.0)
    height_bindings: list[TraceHeightBinding] = []

    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "scene"}, {"type": "xy"}]],
        column_widths=[0.52, 0.32],
        subplot_titles=(
            "Surface",
            "Persistence",
        ),
        horizontal_spacing=0.16,
    )

    fig.add_trace(_surface_trace(field, opacity=surface_opacity, values=field.values), row=1, col=1)
    surface_trace_index = len(fig.data) - 1
    if original_values is not None:
        noisy_surface_values = _surface_vertex_values(field, field.values)
        original_surface_values = _surface_vertex_values(field, original_values)
        height_bindings.append(
            TraceHeightBinding(
                trace_index=surface_trace_index,
                trace_type="mesh3d",
                noisy_z=noisy_surface_values,
                original_z=original_surface_values,
                noisy_intensity=noisy_surface_values,
                original_intensity=original_surface_values,
            )
        )

    x_edges, y_edges, z_edges = edge_trace_points(
        triangulation_edges(field.complex),
        field.coordinates,
        field.values,
        lift=0.2 * lift,
    )
    fig.add_trace(
        go.Scatter3d(
            x=x_edges,
            y=y_edges,
            z=z_edges,
            mode="lines",
            line=dict(width=1, color="rgba(80,80,80,0.32)"),
            hoverinfo="skip",
            name="grid",
        ),
        row=1,
        col=1,
    )
    if original_values is not None:
        _, _, original_z_edges = edge_trace_points(
            triangulation_edges(field.complex),
            field.coordinates,
            original_values,
            lift=0.2 * lift,
        )
        height_bindings.append(
            TraceHeightBinding(
                trace_index=len(fig.data) - 1,
                trace_type="scatter3d",
                noisy_z=tuple(z_edges),
                original_z=tuple(original_z_edges),
            )
        )

    height_bindings.extend(
        _add_critical_traces(
            fig,
            field,
            critical_simplexes,
            lift=lift,
            values=field.values,
            original_values=original_values,
        )
    )
    cycle_cocycle_trace_groups, cycle_cocycle_bindings = _add_cycle_cocycle_threshold_traces(
        fig,
        field,
        sequence,
        diagram,
        thresholds,
        cycle_supports,
        cocycle_supports,
        lift=lift,
        show_coreference_cycles=show_coboundary,
        values=field.values,
        original_values=original_values,
    )
    height_bindings.extend(cycle_cocycle_bindings)
    persistence_trace_indices, surface_highlight_indices, persistence_bindings = _add_persistence_traces(
        fig,
        field,
        diagram,
        thresholds,
        lift=lift,
        values=field.values,
        original_values=original_values,
    )
    height_bindings.extend(persistence_bindings)

    cycle_cocycle_trace_indices = [
        trace_index
        for trace_group in cycle_cocycle_trace_groups
        for trace_index in trace_group
    ]
    dynamic_indices = (
        cycle_cocycle_trace_indices
        + persistence_trace_indices
        + surface_highlight_indices
    )
    persistence_slider_steps = []
    for threshold_index, threshold in enumerate(thresholds):
        visible = [True] * len(fig.data)
        for trace_index in dynamic_indices:
            visible[trace_index] = False
        for trace_index in cycle_cocycle_trace_groups[threshold_index]:
            visible[trace_index] = True
        visible[persistence_trace_indices[threshold_index]] = True
        visible[surface_highlight_indices[2 * threshold_index]] = True
        visible[surface_highlight_indices[2 * threshold_index + 1]] = True
        persistence_slider_steps.append(
            dict(
                method="update",
                args=[
                    {"visible": visible},
                    {"title": "Morse persistence demo"},
                ],
                label=_sparse_slider_label(threshold_index, len(thresholds), threshold),
            )
        )

    opacity_values = (0.15, 0.3, 0.45, 0.6, 0.8, 1.0)
    opacity_slider_steps = [
        dict(
            method="restyle",
            args=[{"opacity": [opacity]}, [surface_trace_index]],
            label=f"{opacity:.2g}",
        )
        for opacity in opacity_values
    ]

    fig.update_layout(
        autosize=True,
        title=dict(
            text="Morse persistence demo",
            x=0.5,
            xanchor="center",
            y=0.985,
            yanchor="top",
        ),
        height=850,
        margin=dict(l=12, r=18, t=90, b=175),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.93,
            xanchor="left",
            x=0.49,
            bgcolor="rgba(255,255,255,0.78)",
            bordercolor="rgba(80,80,80,0.24)",
            borderwidth=1,
            font=dict(size=10),
            itemsizing="constant",
            tracegroupgap=4,
        ),
        sliders=[
            dict(
                active=0,
                currentvalue=dict(prefix="Persistence: ", font=dict(size=11)),
                font=dict(size=9),
                x=0.05,
                y=-0.085,
                len=0.86,
                pad=dict(t=24, b=6),
                steps=persistence_slider_steps,
            ),
            dict(
                active=min(range(len(opacity_values)), key=lambda i: abs(opacity_values[i] - surface_opacity)),
                currentvalue=dict(prefix="Opacity: ", font=dict(size=11)),
                font=dict(size=9),
                x=0.06,
                y=-0.24,
                len=0.42,
                pad=dict(t=26, b=4),
                steps=opacity_slider_steps,
            ),
        ],
        scene=dict(
            xaxis_title="x",
            yaxis_title="y",
            zaxis_title="f",
            aspectmode="cube",
            camera=dict(eye=dict(x=1.35, y=-1.55, z=1.15)),
        ),
    )
    if original_values is not None:
        fig.frames = list(_surface_mode_frames(height_bindings))
        fig.update_layout(updatemenus=[_surface_mode_menu()])
    return fig
