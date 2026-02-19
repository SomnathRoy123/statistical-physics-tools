"""CLI entrypoint for orientation correlation analysis and plotting."""

import argparse
import glob
import re
from pathlib import Path

import numpy as np

from pipeline import process_dataset
from plotting import save_dataset_plot, save_metric_vs_x_plot, save_xi_vs_x_by_r_plot


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
        default="png",
        help="Plot output format (default: pdf for publication workflows)",
    )

    # --- PLOTTING OPTIONS ---
    parser.add_argument(
        "--plot-xi-vs-x",
        action="store_true",
        default=True,
        help="Plot xi vs R (enabled by default)",
    )
    parser.add_argument(
        "--plot-tau-O-vs-x",
        action="store_true",
        default=False,
        help="Plot tau_O vs R (disabled by default)",
    )
    parser.add_argument(
        "--no-dataset-plots",
        action="store_true",
        help="Disable per-dataset C_O(t) fit plots",
    )
    return parser


def _extract_value_from_filename(file_path, key):
    """Extract numeric token by key from filenames such as *_R3.3_* or *_X2.1_*."""
    stem = Path(file_path).stem
    tokens = re.split(r"[_\-]", stem)
    number_pattern = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    key = key.upper()

    for token in tokens:
        token_upper = token.upper()
        if token_upper.startswith(key):
            candidate = token[len(key) :]
            if re.fullmatch(number_pattern, candidate):
                return float(candidate)

    return None


def _extract_r_from_filename(file_path):
    """Extract R value from a filename token like R3.3."""
    return _extract_value_from_filename(file_path, "R")


def _write_summary(results, output_csv):
    lines = ["file,R,tau_O,tau_err,A,xi,beta"]
    for row in results:
        lines.append(
            f"{row['name']},{row['R']:.10g},{row['tau_O']:.10g},{row['tau_err']:.10g},"
            f"{row['A']:.10g},{row['xi']:.10g},{row['beta']:.10g}"
        )
    output_csv.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_vs_x_data(results, output_csv):
    lines = ["R,xi,xi_err,tau_O,tau_O_err"]
    for row in results:
        lines.append(
            f"{row['R']:.10g},{row['xi']:.10g},{row['xi_err']:.10g},"
            f"{row['tau_O']:.10g},{row['tau_err']:.10g}"
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
        r_value = _extract_r_from_filename(file)
        try:
            t, c, fit = process_dataset(
                file,
                t_grid,
                col_time=args.col_time,
                col_corr=args.col_corr,
                stretched=not args.simple_exp,
            )

            if not args.no_dataset_plots:
                save_dataset_plot(t, c, fit, output_dir / f"{name}.{args.plot_format}", title=name)

            if r_value is None:
                raise ValueError(f"Could not extract R from file name: {Path(file).name}")

            result = {
                "name": name,
                "R": r_value,
                "tau_O": fit["tau_O"],
                "tau_err": fit["tau_err"],
                "A": fit["A"],
                "xi": fit["xi"],
                "xi_err": fit["errors"][1],
                "beta": fit.get("beta", 1.0),
            }
            results.append(result)
            print(f"{name}: tau_O={fit['tau_O']:.6g} ± {fit['tau_err']:.3g}")
        except Exception as exc:
            print(f"{name} failed: {exc}")

    if results:
        results = sorted(results, key=lambda item: item["R"])
        _write_summary(results, output_dir / "orientation_corr_summary.csv")
        _write_vs_x_data(results, output_dir / "orientation_corr_vs_x_data.csv")

        if args.plot_xi_vs_x:
            save_xi_vs_x_by_r_plot(
                results,
                output_dir / f"orientation_corr_xi_vs_x.{args.plot_format}",
                x_label="R",
            )

        if args.plot_tau_O_vs_x:
            x_vals = [row["R"] for row in results]
            save_metric_vs_x_plot(
                x_vals,
                [row["tau_O"] for row in results],
                [row["tau_err"] for row in results],
                output_dir / f"orientation_corr_tau_O_vs_x.{args.plot_format}",
                y_label="$\\tau_O$",
                title="Orientation relaxation time $\\tau_O$ vs $R$",
            )


if __name__ == "__main__":
    main()
