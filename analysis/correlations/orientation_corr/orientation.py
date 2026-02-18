# plotting.py
import matplotlib.pyplot as plt

def plot_correlation(ax, t, C, fit, color, label):
    step = max(1, len(t)//30)
    ax.plot(t[::step], C[::step], 'o', mfc='none', color=color, label=label)
    ax.plot(t, fit, '-', color=color)

def finalize_corr_plot(ax, tmin, tmax):
    ax.set_xscale('log')
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(1e-3,1.1)
    ax.set_xlabel(r"Time $t$")
    ax.set_ylabel(r"$C(t)$")
    ax.legend()
