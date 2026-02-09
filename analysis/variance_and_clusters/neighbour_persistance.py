import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.spatial import cKDTree

# ==========================================
# 1. CONFIGURATION (EDIT THIS SECTION)
# ==========================================

# --- Path Settings ---
ROOT_DIR      = Path("/home/somnath2/Codes/Trial_102")      # Parent directory
OUTPUT_ROOT   = Path("/home/somnath2/Codes/mean_sqaure_displacement/persistence_results") # Where to save output
EVO_SUBFOLDER = "evo"                                       # Subfolder containing .dat files

# --- Run Settings ---
# List of specific subdirectories to analyze. 
# Set to None if you want to auto-detect all folders.
SUBDIRS = ["TA400_R2.6","TA400_R2.8"]

# --- Physics Constants (From your C++ Code) ---
OUTPUT_EVERY_N = 10000        # static const long DEFAULT_OUTPUT_EVERY_N
INTEGRATION_DT = 1e-5         # static const double DEFAULT_INTEGRATION_DT
DR             = 1.0          # static const double DEFAULT_DR

# --- Derived Units for Normalization ---
# Time elapsed between two consecutive data files (frames)
TIME_PER_FRAME = OUTPUT_EVERY_N * INTEGRATION_DT   # = 0.1 simulation units

# Characteristic Persistence Time (tau_p = 1 / Dr)
# We will use this to normalize the x-axis.
TAU_P          = 1.0 / DR                          # = 1.0 simulation units

# --- Analysis Parameters ---
BOX_SIZE   = 180.00001    # Simulation box length (for Periodic Boundary Conditions)
R_CUT      = 1.5      # Bond cutoff distance (Set this based on your g(r) minimum)
START_TIME = 500      # Start frame (numeric part of filename)
END_TIME   = 2000     # End frame
MAX_LAG    = 500      # Max lag (in frames) to calculate. 500 frames = 50.0 time units.

# --- File Pattern ---
# Matches "time_500.dat", "time_1000.dat", etc.
FILE_REGEX = r"time_(\d+)\.dat"
pat = re.compile(FILE_REGEX)

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================

def list_time_files_in_range(folder: Path):
    """
    Finds and sorts data files based on the numeric timestamp in the filename.
    Returns: List of tuples (time_step, full_path)
    """
    keep = []
    if not folder.exists():
        return []
        
    for fname in os.listdir(folder):
        m = pat.match(fname)
        if m:
            t = int(m.group(1))
            if START_TIME <= t <= END_TIME:
                keep.append((t, folder / fname))
    
    # Sort files naturally by time step
    keep.sort(key=lambda x: x[0])
    return keep

def get_bonds_for_frame(filepath, box_size, r_cut):
    """
    Reads X,Y data, finds neighbors using Periodic Boundary Conditions (PBC),
    and returns a set of unique bond tuples.
    """
    try:
        # Load data. Assume columns are [X, Y] or [X, Y, Theta]
        data = np.loadtxt(filepath)
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return set()

    if data.ndim == 1:
        data = data.reshape(1, -1)
    
    # Extract only X and Y (first 2 columns)
    positions = data[:, :2]
    n_particles = len(positions)
    
    # cKDTree handles the heavy lifting of finding neighbors + PBC
    # boxsize argument enables the periodic wrap-around logic
    tree = cKDTree(positions, boxsize=box_size)
    
    # Query for all pairs within R_CUT
    # k=12 limits the search to nearest 12 neighbors (optimization)
    distances, indices = tree.query(positions, k=12, distance_upper_bound=r_cut)
    
    bonds = set()
    
    for i, neighbors in enumerate(indices):
        # i is the particle ID (row index)
        for neighbor_idx in neighbors:
            # cKDTree returns N_particles as a placeholder for "no neighbor found"
            if neighbor_idx == n_particles:
                continue
            
            # Store bond as (min, max) tuple to ensure (A, B) == (B, A)
            # Only store if i < neighbor_idx to avoid double counting
            if i < neighbor_idx:
                bond = (i, neighbor_idx)
                bonds.add(bond)
                
    return bonds

def calculate_run_persistence(evo_dir):
    """
    Core algorithm: Calculates bond survival probability P(t) over time.
    """
    # 1. Gather all files
    files = list_time_files_in_range(evo_dir)
    if not files:
        print(f"   [SKIP] No matching files in {evo_dir}")
        return None, None

    print(f"   -> Processing {len(files)} frames (t={files[0][0]} to {files[-1][0]})...")

    # 2. Pre-load all bonds into memory (Frame -> Set of Bonds)
    # This is faster than re-reading files for every lag calculation
    history = {}
    time_steps = []
    
    for t, filepath in files:
        bonds = get_bonds_for_frame(filepath, BOX_SIZE, R_CUT)
        history[t] = bonds
        time_steps.append(t)

    # 3. Compute Persistence for specific lags
    # We check lags: 0, 1, 2, ... up to MAX_LAG (in terms of file index steps)
    
    # Assuming files are roughly contiguous (500, 501, 502...) or stepped (500, 510...)
    # We calculate the lag in "number of frames" relative to the sorted list
    persistence_results = []
    lag_frames_list = []
    
    # Create indices for the files list [0, 1, 2, ... N]
    file_indices = np.arange(len(time_steps))
    
    # We loop over "index lag" (1 file apart, 2 files apart, etc.)
    # Be careful: This assumes constant time spacing between files.
    max_index_lag = min(MAX_LAG, len(time_steps) - 1)
    
    for index_lag in range(max_index_lag + 1):
        surviving_count = 0
        total_start_count = 0
        
        # Slide the window across the trajectory
        # Compare file[i] vs file[i + lag]
        for i in range(len(time_steps) - index_lag):
            t_start = time_steps[i]
            t_end   = time_steps[i + index_lag]
            
            bonds_start = history[t_start]
            bonds_end   = history[t_end]
            
            if len(bonds_start) > 0:
                # Intersection: Bonds that exist at START and END
                preserved = bonds_start.intersection(bonds_end)
                
                surviving_count += len(preserved)
                total_start_count += len(bonds_start)
        
        if total_start_count > 0:
            P_bond = surviving_count / total_start_count
            
            # Calculate the actual numerical time difference for this lag
            # (Just taking the difference of the first pair is usually sufficient if constant dt)
            actual_time_diff_frames = time_steps[index_lag] - time_steps[0]
            
            persistence_results.append(P_bond)
            lag_frames_list.append(actual_time_diff_frames)

    return np.array(lag_frames_list), np.array(persistence_results)

# ==========================================
# 3. MAIN EXECUTION
# ==========================================

def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    
    # Initialize Plot
    plt.figure(figsize=(8, 6))
    
    # Determine which directories to process
    if SUBDIRS is None:
        run_list = sorted([p.name for p in ROOT_DIR.iterdir() if (p / EVO_SUBFOLDER).is_dir()])
    else:
        run_list = SUBDIRS

    print(f"Starting analysis for {len(run_list)} runs...")
    print(f"Normalization parameters: Dr={DR}, dt={INTEGRATION_DT}, OutputEvery={OUTPUT_EVERY_N}")
    print(f"Time per frame = {TIME_PER_FRAME:.4f}, Tau_p = {TAU_P:.4f}")

    for run_name in run_list:
        run_dir = ROOT_DIR / run_name
        evo_dir = run_dir / EVO_SUBFOLDER
        
        print(f"\n[ANALYZING] {run_name}")
        
        if not evo_dir.is_dir():
            print(f"   [WARN] Evo directory not found: {evo_dir}")
            continue

        lags_raw_units, p_bond = calculate_run_persistence(evo_dir)
        
        if lags_raw_units is not None and len(lags_raw_units) > 0:
            # --- TIME NORMALIZATION ---
            # 1. Convert file-number lag to Simulation Time
            # (Note: lags_raw_units is already t_end - t_start from the filename numbers)
            # In your filename "time_500.dat", 500 is likely the step count / output_N?
            # Or is it just an index?
            # CASE A: If filename "500" means "Step 500":
            # physical_time = lags_raw_units * INTEGRATION_DT
            
            # CASE B: If filename "500" is an index (0, 1, 2...) OR if it is "Frame Number":
            # Based on your constants, let's assume the filename number corresponds to
            # the loop iteration or steps / 10000? 
            # Safest approach: Multiply the numeric difference by your known time-per-unit.
            
            # Assuming "time_500.dat" means "File index 500" or similar sequential counter
            # If your files are: time_1.dat, time_2.dat... then lag=1 means 1 * TIME_PER_FRAME
            
            # Let's assume the calculate_run_persistence returns the difference in the NUMBER
            # found in the filename.
            
            # If filename is just an index (1, 2, 3...)
            physical_time = lags_raw_units * TIME_PER_FRAME
            
            # If filename is actual Steps (10000, 20000...) change this line to:
            # physical_time = lags_raw_units * INTEGRATION_DT
             
            # 2. Normalize by Active Timescale (1/Dr)
            normalized_time = physical_time / TAU_P
            
            # --- SAVE DATA ---
            out_file = OUTPUT_ROOT / f"persistence_{run_name}.txt"
            header_txt = (f"Run: {run_name}\n"
                          f"Dr: {DR}, dt: {INTEGRATION_DT}, SaveRate: {OUTPUT_EVERY_N}\n"
                          f"Columns: Lag_Raw_Diff  Physical_Time  Normalized_Time(t*Dr)  P_bond")
            
            np.savetxt(out_file, 
                       np.column_stack((lags_raw_units, physical_time, normalized_time, p_bond)),
                       header=header_txt)
            print(f"   [SAVED] {out_file}")
            
            # --- PLOT ---
            # Plot against Normalized Time
            plt.plot(normalized_time, p_bond, label=run_name, linewidth=2, marker='o', markersize=3, markevery=5)

    # --- FINALIZE PLOT ---
    plt.xlabel(r"Normalized Time $t \times D_r$", fontsize=12)
    plt.ylabel(r"Neighbor Persistence $P_{bond}(t)$", fontsize=12)
    plt.title(f"Bond Lifetime Analysis (Dr={DR})", fontsize=14)
    
    # Formatting for publication quality
    plt.ylim(0, 1.05)
    plt.xlim(left=0)
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.legend(title="Run ID")
    
    # Save Plot
    plot_path = OUTPUT_ROOT / "persistence_comparison_Dr_normalized.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"\n[DONE] Plot saved to: {plot_path}")
    
    # Show Plot
    plt.show()

if __name__ == "__main__":
    main()