# main.py
import numpy as np
import glob
import argparse
from pipeline import process_dataset

parser = argparse.ArgumentParser(description="Orientation correlation analysis")
parser.add_argument("data_path", help="Path to folder containing .dat files")
parser.add_argument("--tmax", type=float, default=None, help="Maximum time for grid")

args = parser.parse_args()

# --- find files ---
files = sorted(glob.glob(f"{args.data_path}/*.dat"))

if len(files) == 0:
    raise RuntimeError(f"No .dat files found in {args.data_path}")

# --- time grid ---
if args.tmax is None:
    Tmax = float(input("Tmax: "))
else:
    Tmax = args.tmax

t_grid = np.logspace(-1, np.log10(Tmax), 500)

# --- run pipeline ---
for file in files:
    try:
        t, C, res = process_dataset(file, t_grid)
        print(file, "tau_O =", res["tau_O"])
    except Exception as e:
        print(file, "failed:", e)
