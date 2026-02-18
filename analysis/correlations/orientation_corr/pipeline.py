# pipeline.py
import numpy as np
from preprocessing import map_to_shared_grid
from fitting import fit_correlation, compute_tau_O

def process_dataset(file_path, t_grid, col_time=1, col_corr=2):

    # ---- Load data ----
    data = np.loadtxt(file_path)

    if data.shape[1] <= max(col_time, col_corr):
        raise ValueError(f"{file_path} has insufficient columns")

    t_raw = data[:, col_time]
    C_raw = data[:, col_corr]

    # ---- Interpolate ----
    t, C = map_to_shared_grid(t_raw, C_raw, t_grid)

    if len(t) < 5:
        raise ValueError("Not enough data after interpolation")

    # ---- Fit correlation ----
    fit = fit_correlation(t, C, stretched=True)

    # ---- Compute tau_O ----
    xi = fit["xi"]
    beta = fit["beta"]
    xi_err = fit["errors"][1]
    beta_err = fit["errors"][2]

    tau, tau_err = compute_tau_O(xi, beta, xi_err, beta_err)

    fit["tau_O"] = tau
    fit["tau_err"] = tau_err

    return t, C, fit
