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
        rf"$\tau_O={fit['tau_O']:.3g}\pm{fit['tau_err']:.2g}$, "
        rf"$\xi={fit['xi']:.3g}$, $\beta={fit.get('beta', 1.0):.3g}$"
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

    ax.set_xlabel(r"$R$")
    ax.set_ylabel(y_label)
    ax.set_title(title)

    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(axis="both", which="major", alpha=0.2, linewidth=0.8)
    ax.grid(axis="both", which="minor", alpha=0.08, linewidth=0.5)

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)
