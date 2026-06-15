import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from data import TurbulentFlatData

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(os.path.dirname(FILE_PATH), "data")

# Parse command line arguments
parser = argparse.ArgumentParser(description="Plot data from single phase turbulent flat validation case.")
parser.add_argument("--Re", type=int, default=180, help="Reynolds number (180, 395, or 934)")
parser.add_argument("--plot-type", type=str, default="both", 
                    choices=["mean-velocity", "rms-velocities", "both"],
                    help="Which plots to generate")
args = parser.parse_args()

# Load data dynamically based on Re number
data_handler = TurbulentFlatData(args.Re)

# Mapping for DNS author labels
DNS_LABELS = {
    180: "DNS - Kim et al.",
    395: "DNS - Moser et al.",
    934: "DNS - Hoyas et al."
}

def plot_mean_velocity():
    """Plot mean velocity profile."""
    data = data_handler.get_mean_velocity()
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Only plot every other data point
    y_plus = data.iloc[::2, 0]
    u_plus = data.iloc[::2, 1]
    
    dns_label = DNS_LABELS.get(args.Re, f"DNS - Re{args.Re}")
    ax.plot(y_plus, u_plus, 'ko', label=dns_label, markersize=8, fillstyle='none', markeredgewidth=1)
    
    # Add reference lines with appropriate ranges
    y_all = data.iloc[:, 0]
    y_viscous = y_all[y_all < 1.34e1]
    ax.plot(y_viscous, y_viscous, 'k--', linewidth=1.5, label='$u^+ = y^+$')
    
    y_loglaw = y_all[y_all > 9]
    ax.plot(y_loglaw, 2.5 * np.log(y_loglaw) + 5.5, 'k-.', linewidth=1.5, label='$u^+ = 2.5 \\ln(y^+) + 5.5$')    
    
    ax.set_xlabel('$y^+$', fontsize=16)
    ax.set_ylabel('$u^+$', fontsize=16)
    ax.set_xscale('log')
    ax.set_title(f'Mean Velocity Profile - Re$_\\tau$ = {args.Re}', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # Set y limits based on Re number
    if args.Re == 180:
        ax.set_ylim(0, 20)
    else:
        ax.set_ylim(0, 25)
    
    ax.legend(fontsize=14)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig

def plot_rms_velocities():
    """Plot RMS velocity profiles."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for component in ['U', 'V', 'W']:
        
        data = data_handler.get_rms_velocity(component)
        ax.plot(data.iloc[:, 0], data.iloc[:, 1], 'o', markersize=8, markeredgewidth=1, fillstyle='none', label=f'DNS - ${component}$')
    
    ax.set_xlabel('$y^+$', fontsize=16)
    ax.set_ylabel(r'$u_{rms}$, $v_{rms}$, $w_{rms}$', fontsize=16)
    ax.set_title(f'RMS Velocity Profiles - Re$_\\tau$ = {args.Re}', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(fontsize=14)
    ax.set_ylim(0, 3)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig

# Generate plots
if args.plot_type in ["mean-velocity", "both"]:
    fig_mean = plot_mean_velocity()
    plt.show()

if args.plot_type in ["rms-velocities", "both"]:
    fig_rms = plot_rms_velocities()
    plt.show()    

