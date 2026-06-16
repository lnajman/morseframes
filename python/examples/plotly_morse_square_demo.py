#!/usr/bin/env python3
"""Write a two-panel Plotly Morse persistence demo for a noisy sine square."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from morseframes.plotting import (  # noqa: E402
    build_noisy_sine_square,
    plot_morse_surface_with_persistence,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=18, help="number of grid squares per axis")
    parser.add_argument("--noise", type=float, default=0.035, help="Gaussian noise amplitude")
    parser.add_argument("--seed", type=int, default=7, help="random seed")
    parser.add_argument("--algorithm", default="f-max", help="Morse sequence strategy")
    parser.add_argument("--surface-opacity", type=float, default=0.42, help="initial surface opacity")
    parser.add_argument("--slider-steps", type=int, default=12, help="number of persistence thresholds")
    parser.add_argument(
        "--hide-coboundary",
        action="store_true",
        help="hide the coreference/coboundary links",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("morseframes/docs/plotly_morse_square_demo.html"),
        help="HTML output path",
    )
    parser.add_argument(
        "--cdn",
        action="store_true",
        help="load Plotly from the CDN instead of embedding it in the HTML",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    field = build_noisy_sine_square(size=args.size, noise=args.noise, seed=args.seed)
    fig = plot_morse_surface_with_persistence(
        field,
        algorithm=args.algorithm,
        surface_opacity=args.surface_opacity,
        persistence_steps=args.slider_steps,
        show_coboundary=not args.hide_coboundary,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(
        args.output,
        include_plotlyjs="cdn" if args.cdn else True,
        full_html=True,
        config={"responsive": True},
        default_width="100%",
        default_height="850px",
    )
    print(args.output.resolve())


if __name__ == "__main__":
    main()
