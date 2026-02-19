"""Publication-quality plotting helpers for orientation-correlation analysis."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, LogLocator


def _apply_paper_style():
    """Apply a clean, soft-matter-friendly plotting style."""
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
            "axes.linewidth": 1.1,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.spines.top": True,
            "axes.spines.right": True,
            "legend.frameon": False,
            "mathtext.default": "regular",
        }
    )


def save_dataset_plot(t, c, fit, output_path, title=None):
    """Save one dataset plot with data and fitted curve in paper style."""
    _apply_paper_style()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5.2, 3.8))

    data_color = "#4C72B0"
    fit_color = "#D62728"

    ax.plot(
        t,
        c,
        "o",
        ms=3.6,
        mew=0.35,
        mec="white",
        alpha=0.9,
        color=data_color,
        label="data",
        zorder=3,
    )
    ax.plot(t, fit["fit"], "-", lw=2.0, color=fit_color, label="fit", zorder=4)

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=6))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$C_O(t)$")

    

    subtitle = (
    f"$\\tau_O={fit['tau_O']:.3g}\\pm{fit['tau_err']:.2g}$, "
    f"$\\xi={fit['xi']:.3g}$, $\\beta={fit.get('beta', 1.0):.3g}$"
    )
    
    ax.set_title(title or "Orientation autocorrelation")
    ax.text(
        0.03,
        0.04,
        subtitle,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.75, "linewidth": 0},
    )

    ax.grid(which="major", alpha=0.2, linewidth=0.8)
    ax.grid(which="minor", alpha=0.08, linewidth=0.5)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_combined_dataset_plot(datasets, output_path, title="Orientation autocorrelation"):
    """Save all orientation-correlation-vs-time datasets in one publication-style plot."""
    _apply_paper_style()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.6, 4.8))

    cmap = plt.cm.viridis
    n = max(len(datasets), 1)

    for idx, data in enumerate(datasets):
        color = cmap(0.12 + 0.78 * (idx / max(n - 1, 1)))
        t = data["t"]
        c = data["c"]
        fit_curve = data["fit"]["fit"]
        r_value = data["R"]

        ax.plot(t, c, "o", ms=3.4, mew=0.3, mec="white", alpha=0.35, color=color, zorder=2)
        ax.plot(t, fit_curve, "-", lw=2.0, color=color, alpha=0.95, label=f"R = {r_value:g}", zorder=3)

    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=7))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$C_O(t)$")
    ax.set_title(title)

    ax.grid(which="major", alpha=0.2, linewidth=0.8)
    ax.grid(which="minor", alpha=0.08, linewidth=0.5)
    ax.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="#cccccc", title=None)

    fig.tight_layout()
    fig.savefig(output_path, dpi=600)
    plt.close(fig)


def save_metric_vs_x_plot(x_values, y_values, y_errors, output_path, y_label, title):
    """Save generic paper-style metric-vs-x plot with error bars."""
    _apply_paper_style()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5.4, 3.9))
    ax.errorbar(
        x_values,
        y_values,
        yerr=y_errors,
        fmt="o",
        color="#222222",
        ecolor="#222222",
        elinewidth=1.1,
        capsize=3,
        capthick=1.1,
        markersize=4.8,
        zorder=4,
    )
    ax.plot(x_values, y_values, "-", lw=1.0, color="#999999", alpha=0.8, zorder=2)

    ax.set_xticks(x_values)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(axis="y", which="major", alpha=0.2, linewidth=0.8)
    ax.grid(axis="y", which="minor", alpha=0.08, linewidth=0.5)

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def save_xi_vs_x_by_r_plot(results, output_path, x_label="R"):
    """Save a single publication-style ξ-vs-X plot where X-axis values are R."""
    _apply_paper_style()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.2, 4.4))

    sorted_results = sorted(results, key=lambda row: row["R"])
    r_vals = [row["R"] for row in sorted_results]
    xi_vals = [row["xi"] for row in sorted_results]
    xi_errs = [row["xi_err"] for row in sorted_results]

    cmap = plt.cm.viridis
    colors = [cmap(0.15 + 0.7 * (idx / max(len(sorted_results) - 1, 1))) for idx in range(len(sorted_results))]

    ax.plot(r_vals, xi_vals, "-", color="#3a3a3a", linewidth=1.4, alpha=0.7, zorder=2)

    for r_val, xi_val, xi_err, color in zip(r_vals, xi_vals, xi_errs, colors):
        ax.errorbar(
            [r_val],
            [xi_val],
            yerr=[xi_err],
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=1.1,
            capsize=3,
            capthick=1.0,
            markersize=5.8,
            label=f"R = {r_val:g}",
            zorder=4,
        )

    ax.set_xlabel(x_label)
    ax.set_ylabel(r"$\xi$")
    ax.set_title(r"Correlation length $\xi$ vs $R$")

    ax.grid(which="major", alpha=0.22, linewidth=0.85)
    ax.grid(which="minor", alpha=0.08, linewidth=0.5)
    ax.minorticks_on()

    ax.legend(loc="upper right", frameon=True, framealpha=0.92, edgecolor="#cccccc", title=None)

    fig.tight_layout()
    fig.savefig(output_path, dpi=600)
    plt.close(fig)
