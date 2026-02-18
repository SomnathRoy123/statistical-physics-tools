# config.py
from dataclasses import dataclass

@dataclass
class GridConfig:
    t_min: float = 0.1
    n_points: int = 500

@dataclass
class FitConfig:
    stretched: bool = True
    calculate_tau: bool = True
    initial_beta: float = 0.8

@dataclass
class ColumnConfig:
    time: int = 1
    corr: int = 2
