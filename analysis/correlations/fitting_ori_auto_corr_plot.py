import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from scipy.optimize import curve_fit
from scipy.special import gamma, psi

# ------------------- CONFIGURATION -------------------

directory_path = '/home/somnath2/Codes/polar_order/orient_time_corr/extra/'
output_filename = '/home/somnath2/Codes/polar_order/orientation_corr_plot_160_2.1_4_555.png'
output_plot_directory = '/home/somnath2/Codes/plots/' # Directory for xi and tau_O plots

start_row = 1
end_row = None  # None means use full data

block_size = 1 # Block size for averaging

# --- CHOOSE FITTING AND CALCULATION OPTIONS ---

# Set to True for stretched exponential fit: C(r) = A * exp(-(r / xi)**beta)
USE_STRETCHED_EXPONENTIAL = True

# Set to True to calculate the average orientation correlation time (τO).
# This requires USE_STRETCHED_EXPONENTIAL to be True.
CALCULATE_TAU_O = True

# --- NEW: ACTIVATE PLOTTING OPTIONS ---
# Set to True to plot the fitted correlation length (ξ) vs. a user-defined x-axis.
plot_xi_vs_x = False
# Set to True to plot the calculated average time (τO) vs. a user-defined x-axis.
plot_tau_O_vs_x = True
# ---------------------------------------------------

# --- SCRIPT VALIDATION ---
if (CALCULATE_TAU_O or plot_tau_O_vs_x) and not USE_STRETCHED_EXPONENTIAL:
    print("Error: To calculate or plot τO, you must set USE_STRETCHED_EXPONENTIAL = True.")
    print("The parameter 'beta' is required.")
    exit()

# --- FUNCTION DEFINITIONS ---
def log_binned_average(x, y, num_bins=100):
    """
    Calculates the average of y-values within logarithmically spaced bins of x.
    """
    # Ensure x starts from a positive value for log spacing
    min_x = x[x > 0].min()
    max_x = x.max()
    
    # Create logarithmically spaced bin edges
    bins = np.logspace(np.log10(min_x), np.log10(max_x), num=num_bins)
    
    # Use numpy.digitize to find which bin each x-value belongs to
    bin_indices = np.digitize(x, bins)
    
    # Calculate the mean for each bin
    binned_x = [x[bin_indices == i].mean() for i in range(1, len(bins))]
    binned_y = [y[bin_indices == i].mean() for i in range(1, len(bins))]
    
    # Filter out empty bins which result in NaN
    binned_x = [val for val in binned_x if not np.isnan(val)]
    binned_y = [val for val in binned_y if not np.isnan(val)]
    
    return np.array(binned_x), np.array(binned_y)



def exponential_fit(r, A, xi):
    """Normal exponential decay function C(r) = A * exp(-r / xi)"""
    return A * np.exp(-r / xi)

def stretched_exponential_fit(r, A, xi, beta):
    """Stretched exponential decay function C(r) = A * exp(-(r / xi)**beta)"""
    return A * np.exp(-(r / xi)**beta)

def block_average(r, C_r, block_size=1):
    """Perform block averaging on r and C_r."""
    n_blocks = len(C_r) // block_size
    r_trim = r[:n_blocks * block_size]
    C_trim = C_r[:n_blocks * block_size]
    r_blocks = r_trim.reshape((n_blocks, block_size))
    C_blocks = C_trim.reshape((n_blocks, block_size))
    r_avg = np.mean(r_blocks, axis=1)
    C_avg = np.mean(C_blocks, axis=1)
    return r_avg, C_avg

def calculate_tau_O(xi, beta, xi_err, beta_err):
    """
    Calculates the average orientation correlation time (τO) and its error.
    """
    try:
        tau_O = (xi / beta) * gamma(1 / beta)
        d_dxi = (1 / beta) * gamma(1 / beta)
        term1 = - (xi / beta**2) * gamma(1 / beta)
        term2 = - (xi / beta**3) * gamma(1 / beta) * psi(1 / beta)
        d_dbeta = term1 + term2
        tau_O_err = np.sqrt((d_dxi**2 * xi_err**2) + (d_dbeta**2 * beta_err**2))
        return tau_O, tau_O_err
    except Exception as e:
        print(f"  Could not calculate Tau_O. Error: {e}")
        return np.nan, np.nan

# --- MAIN SCRIPT ---

file_list = sorted(glob.glob(os.path.join(directory_path, '*.dat')))

if not file_list:
    print(f"Error: No '.dat' files found in '{directory_path}'")
    exit()

print(f"Found {len(file_list)} files. Proceeding to process them...\n")

plot_details = []
for dat_file in file_list:
    filename = os.path.basename(dat_file)
    print(f"\nFile: {filename}")
    try:
        T0_val = float(input(f"  Enter value for T_0 for '{filename}': ").strip())
        R_val = float(input(f"  Enter value for R for '{filename}': ").strip())
    
    # Create the LaTeX formatted legend label
        legend_label = fr"$\mathcal{{T}}_0 = {T0_val}, R = {R_val}$"
        legend_color = input(f"  Enter color for '{filename}': ").strip()

        legend_order = int(input(f"  Enter plotting order for '{filename}': ").strip())
        plot_details.append({
            'path': dat_file, 'label': legend_label, 'color': legend_color, 'order': legend_order
        })
    except Exception as e:
        print(f"  Error reading input for {filename}: {e}")
        continue

plot_details.sort(key=lambda x: x['order'])

# Initialize lists for secondary plots
x_axis_values = []
xi_values = []
tau_O_values = []
tau_O_err_values = []

if plot_xi_vs_x or plot_tau_O_vs_x:
    print("\nYou have activated secondary plotting (ξ and/or τO vs x-axis).")

plt.figure(figsize=(12, 8))

for details in plot_details:
    file_path = details['path']
    label = details['label']
    color = details['color']

    # Initialize loop-specific variables
    tau_O_this_loop, tau_O_err_this_loop = np.nan, np.nan
    xi_fit_this_loop = np.nan
    
    try:
        column1, column2 = np.loadtxt(file_path, unpack=True)
        # r, C_r = block_average(cou, C_r, block_size=block_size)

        r = column1[start_row - 1:end_row]
        C_r = column2[start_row - 1:end_row]


        r_for_fit = r - r[0]
        if len(r_for_fit) < 2:
            print(f"Skipping {label} due to insufficient data points.")
            continue
            
        A_guess = C_r[0]
        # xi_guess = (r[-1] - r[0]) / 3 if len(r) > 1 else 1.0
        xi_guess = 280


        if not USE_STRETCHED_EXPONENTIAL:
            print(f"\n{label}: Using NORMAL exponential fit...")
            popt, pcov = curve_fit(exponential_fit, r_for_fit, C_r, p0=[A_guess, xi_guess], maxfev=5000)
            A_fit, xi_fit = popt
            A_err, xi_err = np.sqrt(np.diag(pcov))
            fitted_values_for_gof = exponential_fit(r_for_fit, A_fit, xi_fit)
            print(f"   Fitted \u03BE = {xi_fit:.3f} \u00B1 {xi_err:.3f}")
            xi_fit_this_loop = xi_fit

        else: # Use stretched exponential
            print(f"\n{label}: Using STRETCHED exponential fit...")
            popt, pcov = curve_fit(stretched_exponential_fit, r_for_fit, C_r, p0=[A_guess, xi_guess, 0.49], maxfev=5000)
            A_fit, xi_fit, beta_fit = popt
            A_err, xi_err, beta_err = np.sqrt(np.diag(pcov))
            fitted_values_for_gof = stretched_exponential_fit(r_for_fit, A_fit, xi_fit, beta_fit)
            print(f"   Fitted \u03BE = {xi_fit:.3f} \u00B1 {xi_err:.3f}")
            print(f"   Fitted \u03B2 = {beta_fit:.3f} \u00B1 {beta_err:.3f}")
            xi_fit_this_loop = xi_fit
            
            if CALCULATE_TAU_O:
                tau_O_this_loop, tau_O_err_this_loop = calculate_tau_O(xi_fit, beta_fit, xi_err, beta_err)
                print(f"   \033[94mCalculated \u03C4O = {tau_O_this_loop:.3f} \u00B1 {tau_O_err_this_loop:.3f}\033[0m")

        residuals = C_r - fitted_values_for_gof
        r_squared = 1 - (np.sum(residuals**2) / np.sum((C_r - np.mean(C_r))**2))
        print(f"   R^2 = {r_squared:.4f}")

        # MODIFIED: Ask for x-value if any secondary plot is enabled
        if plot_xi_vs_x or plot_tau_O_vs_x:
            filename = os.path.basename(file_path)
            while True:
                try:
                    x_val = float(input(f"  Enter x-axis value for '{filename}': ").strip())
                    x_axis_values.append(x_val)
                    break
                except ValueError:
                    print("  Invalid input. Please enter a numeric value.")
            
            # Store values for the plots
            if plot_xi_vs_x:
                xi_values.append(xi_fit_this_loop)
            if plot_tau_O_vs_x:
                tau_O_values.append(tau_O_this_loop)
                tau_O_err_values.append(tau_O_err_this_loop)
        

   

        # 1. Calculate the log-binned average of the data
        # This remains the crucial step for a clean log-log plot.
        r_avg, C_r_avg = log_binned_average(r_for_fit, C_r, num_bins=100) 

        # 2. Plot ONLY the averaged data. 
        # You can make this transparent if you like, but opaque (alpha=1) is usually clearer.
        # The main legend label is assigned here.
        plt.loglog(r_avg, C_r_avg, marker='o', markersize=5, linestyle='None', color=color, label=f"{label} data", alpha=0.5)
        # 3. Plot the fit line (this remains the same)
        r_plot_fit = np.linspace(r_for_fit[0], r_for_fit[-1], 100000)

        if not USE_STRETCHED_EXPONENTIAL:
            C_fit = exponential_fit(r_plot_fit, *popt)
            fit_label = f"{label} fit ($\\xi$={xi_fit:.2f}, $R^2$={r_squared:.3f})"
        else:
            C_fit = stretched_exponential_fit(r_plot_fit, *popt)
            fit_label = f"{label} fit ($\\xi$={xi_fit:.2f}, $\\beta$={beta_fit:.2f}, $R^2$={r_squared:.3f})"

        plt.loglog(r_plot_fit, C_fit, linestyle='-', color=color, alpha=0.9, label=fit_label)
# --- END MODIFIED PLOTTING LOGIC ---
   
    #     # Plot data and fit on the main graph
    #     plt.loglog(r_for_fit, C_r, marker='.', linestyle='None', color=color, label=f"{label} data")
    #     r_plot_fit = np.linspace(r[0], r[-1], 100000)
    #     r_calc_fit = r_plot_fit - r[0]

    #     if not USE_STRETCHED_EXPONENTIAL:
    #         C_fit = exponential_fit(r_calc_fit, *popt)
    #         fit_label = f"{label} fit ($\\xi$={xi_fit:.2f}, $R^2$={r_squared:.3f})"
    #     else:
    #         C_fit = stretched_exponential_fit(r_calc_fit, *popt)
    #         fit_label = f"{label} fit ($\\xi$={xi_fit:.2f}, $\\beta$={beta_fit:.2f}, $R^2$={r_squared:.3f})"
        
    #     plt.loglog(r_calc_fit, C_fit, linestyle='-', color=color, alpha=0.7, label=fit_label)

    except Exception as e:
        print(f"Error processing {label}: {e}")

# Final formatting for correlation plot
title_suffix = "Stretched Exponential Fit" if USE_STRETCHED_EXPONENTIAL else "Exponential Fit"
# plt.title(f"Orientation Correlation C(r) with {title_suffix}")
plt.xlabel("Correlation Time")
plt.ylabel("Orientation Correlation C(r)")
plt.xscale('log')
plt.yscale('log')
plt.ylim(bottom = 0.025)
# plt.xlim(left = 0.001)
# plt.xlim(right=  1000)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend(title="Legend", loc="best")
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"\nCorrelation plot saved successfully as '{output_filename}'")
plt.show()

# --- SECONDARY PLOTTING SECTION ---

os.makedirs(output_plot_directory, exist_ok=True) # Ensure output directory exists

# Plot xi vs x-axis if activated
if plot_xi_vs_x and xi_values:
    x_plot = np.array(x_axis_values)
    xi_plot = np.array(xi_values)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, xi_plot, marker='o', linestyle='-')
    # plt.title("Fitted Correlation Length (\u03BE) vs User-Defined X-Axis")
    plt.xlabel("Radius of interaction")
    plt.ylabel("\u03BE (relaxation time)")
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.tight_layout()
    output_path = os.path.join(output_plot_directory, 'xi_vs_x.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"ξ vs x-axis plot saved successfully as '{output_path}'")
    plt.show()

# NEW: Plot tau_O vs x-axis if activated
if plot_tau_O_vs_x and tau_O_values:
    # Convert to numpy arrays for filtering
    x_plot = np.array(x_axis_values)
    tau_plot = np.array(tau_O_values)
    err_plot = np.array(tau_O_err_values)

    # Create a mask to filter out NaN values (from non-stretched fits)
    mask = ~np.isnan(tau_plot)
    if not np.any(mask):
        print("\nCould not generate τO plot: No valid τO values were calculated.")
    else:
        plt.figure(figsize=(10, 6))
        plt.errorbar(x_plot[mask], tau_plot[mask], yerr=err_plot[mask], 
                     marker='s', linestyle='--', capsize=5, ecolor='red', label='Calculated τO')
        # plt.title("Average Orientation Time (\beta) vs User-Defined X-Axis")
        plt.yscale('log')
        plt.xlabel("Interaction Radius")
        plt.ylabel("Tau_O(Average Orientation Time)")
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.legend()
        plt.tight_layout()
        output_path = os.path.join(output_plot_directory, 'tau_O_vs_x.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"τO vs x-axis plot saved successfully as '{output_path}'")
        plt.show()