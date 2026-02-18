"""CLI entrypoint for orientation correlation analysis and plotting."""

import argparse
import glob
from pathlib import Path

import numpy as np

from pipeline import process_dataset
from plotting import save_dataset_plot, save_tau_summary_plot


def _build_parser():
    parser = argparse.ArgumentParser(description="Orientation correlation analysis")
    parser.add_argument("data_path", help="Path to folder containing .dat files")
    parser.add_argument("--tmax", type=float, default=None, help="Maximum time for shared grid")
    parser.add_argument("--n-grid", type=int, default=500, help="Number of log-spaced grid points")
    parser.add_argument("--col-time", type=int, default=1, help="Time column index")
    parser.add_argument("--col-corr", type=int, default=2, help="Correlation column index")
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory for plots and summary CSV (default: <data_path>/plots)",
    )
    parser.add_argument(
        "--simple-exp",
        action="store_true",
        help="Use simple exponential fit instead of stretched exponential",
    )
    parser.add_argument(
        "--plot-format",
        choices=["png", "pdf", "svg"],
        default="pdf",
        help="Plot output format (default: pdf for publication workflows)",
    )
    return parser


def _write_summary(results, output_csv):
    lines = ["file,tau_O,tau_err,A,xi,beta"]
    for row in results:
        lines.append(
            f"{row['name']},{row['tau_O']:.10g},{row['tau_err']:.10g},"
            f"{row['A']:.10g},{row['xi']:.10g},{row['beta']:.10g}"
        )
    output_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = _build_parser().parse_args()

    files = sorted(glob.glob(f"{args.data_path}/*.dat"))
    if not files:
        raise RuntimeError(f"No .dat files found in {args.data_path}")

    if args.tmax is None:
        raise ValueError("Please provide --tmax for non-interactive runs")

    t_grid = np.logspace(-1, np.log10(args.tmax), args.n_grid)

    output_dir = Path(args.output_dir or Path(args.data_path) / "plots")
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for file in files:
        name = Path(file).stem
        try:
            t, c, fit = process_dataset(
                file,
                t_grid,
                col_time=args.col_time,
                col_corr=args.col_corr,
                stretched=not args.simple_exp,
            )
            save_dataset_plot(t, c, fit, output_dir / f"{name}.{args.plot_format}", title=name)

            result = {
                "name": name,
                "tau_O": fit["tau_O"],
                "tau_err": fit["tau_err"],
                "A": fit["A"],
                "xi": fit["xi"],
                "beta": fit.get("beta", 1.0),
            }
            results.append(result)
            print(f"{name}: tau_O={fit['tau_O']:.6g} ± {fit['tau_err']:.3g}")
        except Exception as exc:
            print(f"{name} failed: {exc}")

    if results:
        _write_summary(results, output_dir / "orientation_corr_summary.csv")
        save_tau_summary_plot(
            results, output_dir / f"orientation_corr_tau_summary.{args.plot_format}"
        )


if __name__ == "__main__":
    main()
