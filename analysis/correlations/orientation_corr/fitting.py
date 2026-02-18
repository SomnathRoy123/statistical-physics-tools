# fitting.py
import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit
import warnings
from scipy.special import gamma, psi
from models import stretched_exponential, exponential

def compute_tau_O(xi, beta, xi_err, beta_err):
    val = (xi / beta) * gamma(1 / beta)

    d_dxi = (1 / beta) * gamma(1 / beta)
    term1 = - (xi / beta**2) * gamma(1 / beta)
    term2 = - (xi / beta**3) * gamma(1 / beta) * psi(1 / beta)
    d_dbeta = term1 + term2

    err = np.sqrt((d_dxi**2 * xi_err**2) + (d_dbeta**2 * beta_err**2))
    return val, err


def fit_correlation(t, C, stretched=True):
    """Fit orientation correlation with optional stretched exponential model."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)

        if stretched:
            popt, pcov = curve_fit(
                stretched_exponential,
                t,
                C,
                p0=[1.0, t.max() / 5, 0.8],
                bounds=([0, 0, 0.1], [2, np.inf, 2]),
                maxfev=10000,
            )
            A, xi, beta = popt
            errors = np.sqrt(np.maximum(np.diag(pcov), 0))
            return {
                "A": A,
                "xi": xi,
                "beta": beta,
                "errors": errors,
                "fit": stretched_exponential(t, *popt),
            }

        popt, pcov = curve_fit(exponential, t, C, p0=[1.0, t.max() / 5], maxfev=10000)
        return {
            "A": popt[0],
            "xi": popt[1],
            "errors": np.sqrt(np.maximum(np.diag(pcov), 0)),
            "fit": exponential(t, *popt),
        }
