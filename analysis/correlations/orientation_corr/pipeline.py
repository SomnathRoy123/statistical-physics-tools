"""Analysis pipeline for orientation autocorrelation datasets."""

from pathlib import Path

from dataio import load_correlation_file
from fitting import compute_tau_O, fit_correlation
from preprocessing import map_to_shared_grid


def process_dataset(file_path, t_grid, col_time=1, col_corr=2, stretched=True):
    """Load, preprocess, fit and summarize one dataset.

    Returns interpolated time/correlation arrays and fit metadata.
    """
    file_path = Path(file_path)

    # ---- Load data ----
    t_raw, c_raw = load_correlation_file(file_path, col_time, col_corr)

    # ---- Interpolate ----
    t, c = map_to_shared_grid(t_raw, c_raw, t_grid)

    if len(t) < 5:
        raise ValueError(f"{file_path} has insufficient valid points after interpolation")

    # ---- Fit correlation ----
    fit = fit_correlation(t, c, stretched=stretched)

    # ---- Compute tau_O ----
    xi = fit["xi"]
    xi_err = fit["errors"][1]

    if stretched:
        beta = fit["beta"]
        beta_err = fit["errors"][2]
    else:
        beta = 1.0
        beta_err = 0.0

    tau, tau_err = compute_tau_O(xi, beta, xi_err, beta_err)
    fit["tau_O"] = tau
    fit["tau_err"] = tau_err

    return t, c, fit
