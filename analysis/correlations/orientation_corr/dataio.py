# io.py
import numpy as np

def load_correlation_file(path, col_time, col_corr):
    data = np.loadtxt(path)

    if data.shape[1] <= max(col_time, col_corr):
        raise ValueError("Not enough columns")

    return data[:,col_time], data[:,col_corr]
