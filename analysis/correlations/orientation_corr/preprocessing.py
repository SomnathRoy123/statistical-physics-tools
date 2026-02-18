# preprocessing.py
import numpy as np
from scipy.interpolate import interp1d

def map_to_shared_grid(t_raw, C_raw, t_grid):
    t_raw, idx = np.unique(t_raw, return_index=True)
    C_raw = C_raw[idx]

    interp = interp1d(t_raw, C_raw, bounds_error=False, fill_value=np.nan)
    C_grid = interp(t_grid)

    mask = ~np.isnan(C_grid)
    return t_grid[mask], C_grid[mask]
