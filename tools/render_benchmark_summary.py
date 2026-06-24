#!/usr/bin/env python3
"""Render compact Markdown benchmark tables into the public summary page."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = ROOT / "docs" / "benchmark_summary.md"
DEFAULT_SYNTHETIC_TABLE = ROOT / "docs" / "synthetic_scale_table.tex"
DEFAULT_NATIVE_TABLE = ROOT / "docs" / "native_gudhi_view_default_r30_table.tex"
DEFAULT_OVERHEAD_TABLE = ROOT / "docs" / "native_gudhi_overhead_summary_table.tex"
START_MARKER = "<!-- benchmark-summary-tables:start -->"
END_MARKER = "<!-- benchmark-summary-tables:end -->"


@dataclass(frozen=True)
class NativeRow:
    case: str
    strategy: str
    simplices: str
    direct_ms: str
    gudhi_ms: str
    gudhi_over_direct: str
    gudhi_over_reducer: str
    direct_ratio_median: float


def _table_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line.endswith(r"\\") or "&" not in line:
            continue
        if line.startswith(
            (
                r"\toprule",
                r"\midrule",
                r"\bottomrule",
                r"\caption",
                r"\label",
            )
        ) or " & " not in line:
            continue
        line = line[:-2].strip()
        cells = [cell.strip() for cell in line.split("&")]
        if not cells or cells[0] in {"Family", "Case", "Direct-view stage group"}:
            continue
        rows.append(cells)
    return rows


def _tex_to_markdown(text: str) -> str:
    text = text.strip()
    text = re.sub(r"\\texttt\{([^{}]+)\}", r"`\1`", text)
    text = text.replace(r"\slred", "same-level")
    text = text.replace(r"\%", "%")
    text = text.replace(r"\_", "_")
    text = text.replace("{", "").replace("}", "")
    text = text.replace("$", "")
    return text.strip()


def _format_ratio(text: str) -> str:
    match = re.fullmatch(r"\\ratioiqr\{([^{}]+)\}\{([^{}]+)\}\{([^{}]+)\}", text.strip())
    if match is None:
        return _tex_to_markdown(text)
    median, q1, q3 = match.groups()
    return f"{median} ({q1}-{q3})"


def _ratio_median(text: str) -> float:
    match = re.search(r"\{([0-9.]+)\}", text)
    if match is not None:
        return float(match.group(1))
    return float(_tex_to_markdown(text))


def render_synthetic_table(path: Path) -> str:
    family_names = {"rips": "Rips"}
    lines = [
        "### Synthetic Scale",
        "",
        "| Family | Strategy | Cases | Avg. simplices | Critical % | Morse time | Std/Morse |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for cells in _table_rows(path):
        if len(cells) != 7:
            continue
        family, strategy, cases, avg_simplices, critical, morse_us, ratio = cells
        family = family_names.get(_tex_to_markdown(family), _tex_to_markdown(family))
        lines.append(
            "| "
            + " | ".join(
                [
                    family,
                    _tex_to_markdown(strategy),
                    _tex_to_markdown(cases),
                    _tex_to_markdown(avg_simplices),
                    _tex_to_markdown(critical),
                    f"{_tex_to_markdown(morse_us)} us",
                    _format_ratio(ratio),
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def render_native_table(path: Path) -> str:
    best_by_case: dict[str, NativeRow] = {}
    for cells in _table_rows(path):
        if len(cells) != 8:
            continue
        case, strategy, simplices, direct_ms, _import_ms, gudhi_ms, direct_ratio, reducer_ratio = cells
        row = NativeRow(
            case=_tex_to_markdown(case),
            strategy=_tex_to_markdown(strategy),
            simplices=_tex_to_markdown(simplices),
            direct_ms=_tex_to_markdown(direct_ms),
            gudhi_ms=_tex_to_markdown(gudhi_ms),
            gudhi_over_direct=_format_ratio(direct_ratio),
            gudhi_over_reducer=_format_ratio(reducer_ratio),
            direct_ratio_median=_ratio_median(direct_ratio),
        )
        existing = best_by_case.get(row.case)
        if existing is None or row.direct_ratio_median > existing.direct_ratio_median:
            best_by_case[row.case] = row

    lines = [
        "### Native GUDHI View",
        "",
        "This compact table reports the best `GUDHI/Direct` strategy for each default",
        "30-repeat native case.",
        "",
        "| Case | Simplices | Best direct strategy | Direct | GUDHI | GUDHI/Direct | GUDHI/Reducer |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in best_by_case.values():
        lines.append(
            "| "
            + " | ".join(
                [
                    row.case,
                    row.simplices,
                    row.strategy,
                    f"{row.direct_ms} ms",
                    f"{row.gudhi_ms} ms",
                    row.gudhi_over_direct,
                    row.gudhi_over_reducer,
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def render_overhead_table(path: Path) -> str:
    lines = [
        "### Direct-Path Stage Split",
        "",
        "| Stage group | Share | Interpretation |",
        "| --- | ---: | --- |",
    ]
    for cells in _table_rows(path):
        if len(cells) != 3:
            continue
        stage, share, interpretation = cells
        lines.append(
            "| "
            + " | ".join(
                [
                    _tex_to_markdown(stage),
                    _tex_to_markdown(share),
                    _tex_to_markdown(interpretation),
                ]
            )
            + " |"
        )
    return "\n".join(lines)


def render_block(
    *,
    synthetic_table: Path,
    native_table: Path,
    overhead_table: Path,
) -> str:
    sections = [
        START_MARKER,
        "",
        "The tables below are short, rendered summaries of the current public benchmark",
        "fragments. This block is generated by `tools/render_benchmark_summary.py`.",
        "The full LaTeX fragments remain the source for reproducible benchmark tables in",
        "this software repository; they are not manuscript source.",
        "",
        render_synthetic_table(synthetic_table),
        "",
        render_native_table(native_table),
        "",
        render_overhead_table(overhead_table),
        "",
        END_MARKER,
    ]
    return "\n".join(sections)


def _replace_generated_block(text: str, block: str, path: Path) -> str:
    if START_MARKER not in text or END_MARKER not in text:
        raise ValueError(
            f"{path} must contain {START_MARKER!r} and {END_MARKER!r} markers."
        )
    before, rest = text.split(START_MARKER, 1)
    _old, after = rest.split(END_MARKER, 1)
    return before.rstrip() + "\n\n" + block.rstrip() + "\n\n" + after.lstrip()


def update_summary(path: Path, block: str) -> None:
    path.write_text(_replace_generated_block(path.read_text(), block, path))


def check_summary(path: Path, block: str) -> None:
    text = path.read_text()
    expected = _replace_generated_block(text, block, path)
    if text != expected:
        raise SystemExit(
            f"{path} is out of date. Run `python3 tools/render_benchmark_summary.py`."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render compact benchmark tables into docs/benchmark_summary.md."
    )
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--synthetic-table", type=Path, default=DEFAULT_SYNTHETIC_TABLE)
    parser.add_argument("--native-table", type=Path, default=DEFAULT_NATIVE_TABLE)
    parser.add_argument("--overhead-table", type=Path, default=DEFAULT_OVERHEAD_TABLE)
    parser.add_argument("--stdout", action="store_true")
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if the generated block in --summary is out of date.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    block = render_block(
        synthetic_table=args.synthetic_table,
        native_table=args.native_table,
        overhead_table=args.overhead_table,
    )
    if args.stdout:
        print(block)
    elif args.check:
        check_summary(args.summary, block)
    else:
        update_summary(args.summary, block)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
