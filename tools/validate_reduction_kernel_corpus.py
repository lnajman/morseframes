#!/usr/bin/env python3
"""Validate deterministic reduction-kernel behavior on a fixed corpus."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Callable, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import morseframes as mf  # noqa: E402


DEFAULT_CORPUS = ROOT / "tests" / "data" / "reduction_kernel_corpus.json"
SEQUENTIAL_ALGORITHM = mf.FLOODING_REDUCTION_KERNEL_SEQUENCE
PARALLEL_ALGORITHM = mf.FLOODING_REDUCTION_KERNEL_PARALLEL_SEQUENCE


def _require_list(value: object, *, field: str, case_name: str) -> list[Any]:
    if not isinstance(value, list):
        raise ValueError(f"{case_name}: {field!r} must be a list")
    return value


def load_corpus(path: Path) -> list[dict[str, Any]]:
    """Load and perform structural validation of a corpus manifest."""

    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise ValueError(f"corpus file does not exist: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in {path}: {exc}") from exc

    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise ValueError("corpus must be an object with schema_version 1")
    cases = _require_list(document.get("cases"), field="cases", case_name="corpus")
    if not cases:
        raise ValueError("corpus must contain at least one case")

    names: set[str] = set()
    validated: list[dict[str, Any]] = []
    for index, value in enumerate(cases):
        if not isinstance(value, dict):
            raise ValueError(f"case {index}: expected an object")
        name = value.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError(f"case {index}: name must be a non-empty string")
        if name in names:
            raise ValueError(f"duplicate corpus case name: {name}")
        names.add(name)
        if value.get("kind") not in {"constant_facets", "lower_star", "simplices"}:
            raise ValueError(f"{name}: unsupported kind {value.get('kind')!r}")
        validated.append(value)
    return validated


def build_complex(case: dict[str, Any]) -> mf.FilteredComplex:
    """Construct one finalized complex from its manifest entry."""

    name = str(case["name"])
    kind = case["kind"]
    if kind == "constant_facets":
        facets = _require_list(case.get("facets"), field="facets", case_name=name)
        filtration = float(case.get("filtration", 0.0))
        return mf.FilteredComplex.from_facets(facets, simplex_filtration=filtration)

    if kind == "lower_star":
        facets = _require_list(case.get("facets"), field="facets", case_name=name)
        values = _require_list(
            case.get("vertex_filtration"), field="vertex_filtration", case_name=name
        )
        vertex_filtration = {vertex: float(value) for vertex, value in enumerate(values)}
        return mf.FilteredComplex.from_lower_star(
            facets,
            vertex_filtration,
            dimension_offset=float(case.get("dimension_offset", 0.0)),
        )

    entries = _require_list(case.get("simplices"), field="simplices", case_name=name)
    simplices: list[tuple[Sequence[int], float]] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, list) or len(entry) != 2:
            raise ValueError(f"{name}: simplices[{index}] must be [vertices, filtration]")
        vertices = _require_list(
            entry[0], field=f"simplices[{index}][0]", case_name=name
        )
        simplices.append((vertices, float(entry[1])))
    return mf.FilteredComplex.from_simplices(simplices)


def _step_signature(
    complex_: mf.FilteredComplex, sequence: mf.MorseSequence
) -> tuple[tuple[object, ...], ...]:
    return tuple(
        (
            step.type,
            complex_.vertices(step.sigma),
            None if step.tau is None else complex_.vertices(step.tau),
            step.level,
        )
        for step in sequence.steps
    )


def _digest(signature: tuple[tuple[object, ...], ...]) -> str:
    payload = json.dumps(signature, separators=(",", ":"), sort_keys=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _barcode(diagram: mf.PersistenceDiagram) -> tuple[object, object]:
    # Morse cancellation may remove zero-persistence pairs, so the invariant is
    # the usual barcode with those diagonal intervals omitted.
    return diagram.finite_barcode(), diagram.essential_barcode()


def validate_case(
    case: dict[str, Any],
    workers: Sequence[int],
    *,
    progress: Callable[[str], None] | None = None,
) -> dict[str, object]:
    """Validate one case and return its compact report record."""

    name = str(case["name"])
    complex_ = build_complex(case)
    if progress is not None:
        progress(f"{name}: sequential sequence")
    sequential = mf.compute_morse_sequence(
        complex_, algorithm=SEQUENTIAL_ALGORITHM
    )
    expected_signature = _step_signature(complex_, sequential)
    if progress is not None:
        progress(f"{name}: persistence oracle")
    standard_barcode = _barcode(mf.compute_standard_persistence(complex_))
    sequential_barcode = _barcode(
        mf.compute_morse_persistence(complex_, sequence=sequential)
    )
    if sequential_barcode != standard_barcode:
        raise AssertionError(
            f"{name}: sequential Morse persistence differs from the standard reducer\n"
            f"  sequential={sequential_barcode}\n  standard={standard_barcode}"
        )

    for worker_count in workers:
        if progress is not None:
            progress(f"{name}: parallel sequence with workers={worker_count}")
        parallel = mf.compute_morse_sequence(
            complex_, algorithm=PARALLEL_ALGORITHM, max_workers=worker_count
        )
        actual_signature = _step_signature(complex_, parallel)
        if actual_signature != expected_signature:
            raise AssertionError(
                f"{name}: ordered Morse steps differ for max_workers={worker_count}"
            )
        parallel_barcode = _barcode(
            mf.compute_morse_persistence(complex_, sequence=parallel)
        )
        if parallel_barcode != standard_barcode:
            raise AssertionError(
                f"{name}: parallel Morse persistence differs from the standard reducer "
                f"for max_workers={worker_count}"
            )

    return {
        "name": name,
        "simplices": complex_.size,
        "levels": complex_.num_levels,
        "steps": len(sequential.steps),
        "critical": len(sequential.critical_simplices),
        "step_digest": _digest(expected_signature),
        "workers": list(workers),
    }


def _parse_workers(value: str) -> tuple[int, ...]:
    try:
        workers = tuple(int(item.strip()) for item in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("workers must be comma-separated integers") from exc
    if not workers or any(worker < 1 for worker in workers):
        raise argparse.ArgumentTypeError("workers must contain positive integers")
    if len(set(workers)) != len(workers):
        raise argparse.ArgumentTypeError("workers must not contain duplicates")
    return workers


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--workers", type=_parse_workers, default=(1, 2, 4))
    parser.add_argument("--json", action="store_true", help="emit a JSON report")
    args = parser.parse_args(argv)

    try:
        cases = load_corpus(args.corpus)
        progress = None if args.json else lambda message: print(message, flush=True)
        results = [
            validate_case(case, args.workers, progress=progress) for case in cases
        ]
    except (AssertionError, TypeError, ValueError) as exc:
        print(f"reduction-kernel corpus validation failed: {exc}", file=sys.stderr)
        return 1

    report = {
        "backend": "native" if mf.cpp_backend_available() else "pure-python",
        "corpus": str(args.corpus),
        "cases": results,
    }
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        for result in results:
            print(
                f"{result['name']}: ok "
                f"(simplices={result['simplices']}, levels={result['levels']}, "
                f"steps={result['steps']}, digest={result['step_digest']})"
            )
        print(
            f"validated {len(results)} cases with {report['backend']} backend "
            f"at workers={','.join(str(worker) for worker in args.workers)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
