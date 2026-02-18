# models.py
import numpy as np

def exponential(t, A, xi):
    return A * np.exp(-t / xi)

def stretched_exponential(t, A, xi, beta):
    t = np.abs(t)
    return A * np.exp(-(t / xi)**beta)
