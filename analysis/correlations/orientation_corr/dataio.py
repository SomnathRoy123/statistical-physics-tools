"""Data loading helpers for correlation files."""

import numpy as np


def load_correlation_file(path, col_time, col_corr):
    data = np.loadtxt(path)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] <= max(col_time, col_corr):
        raise ValueError(f"{path} does not have required columns {col_time}, {col_corr}")

    return data[:, col_time], data[:, col_corr]
