#!/usr/bin/env python3
"""Summarize profile-vs-measured CSVs into candidate-gate proposals."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(frozen=True)
class GateProposal:
    group: str
    cases: int
    candidate_count: int
    coverage_percent: float
    avg_best_us: float
    winners: str
    proposed_candidates: str


def _optional_float(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def read_profile_rows(paths: Iterable[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(newline="") as handle:
            rows.extend(csv.DictReader(handle))
    return rows


def shape_bin_from_values(
    *,
    max_dimension: int,
    edge_density: float,
    triangle_density: float,
    largest_level_ratio: float,
    max_coboundary: int,
) -> str:
    if max_dimension < 0:
        return "unknown-shape"
    if max_dimension <= 1:
        return "graph-sparse" if edge_density < 0.1 else "graph-dense"
    if max_dimension == 2 and edge_density >= 0.9 and triangle_density >= 0.9:
        return "complete-2-skeleton"
    if max_dimension >= 3:
        return "higher-dimensional"
    if max_coboundary >= 100:
        return "high-coboundary"
    if largest_level_ratio >= 0.5:
        return "plateau-heavy"
    if edge_density >= 0.45 or triangle_density >= 0.08:
        return "dense-2d"
    return "moderate-2d"


def shape_bin(row: dict[str, str]) -> str:
    max_dimension = int(_optional_float(row, "shape_max_dimension", -1.0))
    edge_density = _optional_float(row, "shape_edge_density")
    triangle_density = _optional_float(row, "shape_triangle_density")
    largest_level_ratio = _optional_float(row, "shape_largest_level_ratio")
    max_coboundary = int(_optional_float(row, "shape_max_coboundary_size"))

    return shape_bin_from_values(
        max_dimension=max_dimension,
        edge_density=edge_density,
        triangle_density=triangle_density,
        largest_level_ratio=largest_level_ratio,
        max_coboundary=max_coboundary,
    )


def group_key(row: dict[str, str], group_by: str) -> str:
    normalized = group_by.lower().replace("_", "-")
    family = row.get("family", "unknown-family")
    size = row.get("size", "unknown-size")
    bin_name = shape_bin(row)
    if normalized == "all":
        return "all"
    if normalized == "family":
        return family
    if normalized == "family-size":
        return f"{family}:n={size}"
    if normalized == "shape-bin":
        return bin_name
    if normalized == "family-shape":
        return f"{family}:{bin_name}"
    raise ValueError(
        "Unknown grouping "
        f"{group_by!r}. Supported: all, family, family-size, shape-bin, family-shape."
    )


def propose_candidate_gates(
    rows: Iterable[dict[str, str]],
    *,
    group_by: str = "family-shape",
    max_candidates: int | None = None,
    include_candidates: dict[str, tuple[str, ...]] | None = None,
) -> list[GateProposal]:
    include_candidates = include_candidates or {}
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if not row.get("measured_best_algorithm"):
            continue
        grouped[group_key(row, group_by)].append(row)

    proposals: list[GateProposal] = []
    for key, group in sorted(grouped.items()):
        winners = Counter(row["measured_best_algorithm"] for row in group)
        ordered = [algorithm for algorithm, _ in winners.most_common()]
        selected = ordered if max_candidates is None else ordered[: max(1, max_candidates)]
        candidates = list(selected)
        for guard in include_candidates.get(key, ()) + include_candidates.get("all", ()):
            if guard not in candidates:
                candidates.append(guard)
        candidates = tuple(candidates)
        covered = sum(1 for row in group if row["measured_best_algorithm"] in candidates)
        avg_best_us = (
            1.0e6
            * sum(_optional_float(row, "measured_best_morse_seconds") for row in group)
            / len(group)
        )
        proposals.append(
            GateProposal(
                group=key,
                cases=len(group),
                candidate_count=len(candidates),
                coverage_percent=100.0 * covered / len(group),
                avg_best_us=avg_best_us,
                winners=";".join(f"{algorithm}:{count}" for algorithm, count in winners.most_common()),
                proposed_candidates=";".join(candidates),
            )
        )
    return proposals


def parse_include_candidates(specs: Iterable[str]) -> dict[str, tuple[str, ...]]:
    parsed: dict[str, list[str]] = defaultdict(list)
    for spec in specs:
        if "=" not in spec:
            raise ValueError(
                f"Invalid include-candidates value {spec!r}; expected GROUP=ALG[,ALG...]"
            )
        group, algorithms = spec.split("=", 1)
        group = group.strip()
        if not group:
            raise ValueError(f"Invalid include-candidates value {spec!r}: empty group")
        for algorithm in algorithms.split(","):
            algorithm = algorithm.strip()
            if algorithm and algorithm not in parsed[group]:
                parsed[group].append(algorithm)
    return {group: tuple(algorithms) for group, algorithms in parsed.items()}


def write_table(proposals: list[GateProposal], output: TextIO) -> None:
    headers = [
        "group",
        "cases",
        "cands",
        "coverage%",
        "best_us",
        "winners",
        "proposal",
    ]
    output.write(" ".join(f"{header:>18}" for header in headers) + "\n")
    for proposal in proposals:
        output.write(
            f"{proposal.group:>18} "
            f"{proposal.cases:18d} "
            f"{proposal.candidate_count:18d} "
            f"{proposal.coverage_percent:18.1f} "
            f"{proposal.avg_best_us:18.2f} "
            f"{proposal.winners:>18} "
            f"{proposal.proposed_candidates:>18}\n"
        )


def write_csv(proposals: list[GateProposal], output: TextIO) -> None:
    writer = csv.DictWriter(output, fieldnames=list(asdict(proposals[0]).keys()))
    writer.writeheader()
    for proposal in proposals:
        writer.writerow(asdict(proposal))


def write_json(proposals: list[GateProposal], output: TextIO) -> None:
    json.dump([asdict(proposal) for proposal in proposals], output, indent=2)
    output.write("\n")


def write_proposals(
    proposals: list[GateProposal],
    *,
    output_format: str,
    output_path: Path | None,
) -> None:
    if not proposals:
        raise ValueError("No gate proposals were generated.")

    if output_path is None:
        output = None
        handle = None
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        handle = output_path.open("w", newline="")
        output = handle

    try:
        target = output if output is not None else sys.stdout
        if output_format == "table":
            write_table(proposals, target)
        elif output_format == "csv":
            write_csv(proposals, target)
        elif output_format == "json":
            write_json(proposals, target)
        else:
            raise ValueError(f"Unknown output format: {output_format}")
    finally:
        if handle is not None:
            handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument(
        "--group-by",
        default="family-shape",
        choices=["all", "family", "family-size", "shape-bin", "family-shape"],
    )
    parser.add_argument(
        "--max-candidates",
        type=int,
        help="Keep only the most frequent winners in each group.",
    )
    parser.add_argument(
        "--include-candidates",
        action="append",
        default=[],
        metavar="GROUP=ALG[,ALG...]",
        help=(
            "Always include additional guard candidates for a proposal group. "
            "May be passed more than once; use group 'all' for a global guard."
        ),
    )
    parser.add_argument("--format", choices=["table", "csv", "json"], default="table")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = read_profile_rows(args.inputs)
    proposals = propose_candidate_gates(
        rows,
        group_by=args.group_by,
        max_candidates=args.max_candidates,
        include_candidates=parse_include_candidates(args.include_candidates),
    )
    write_proposals(proposals, output_format=args.format, output_path=args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
