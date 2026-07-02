"""
Plotting functions for turbulent channel validation data.

Can be used standalone to plot experimental data:
    python plot_data.py --Re 180

Or imported from post_process.py to overlay simulation results:
    from plot_data import plot_mean_velocity_profile, plot_rms_velocity_profiles
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from data import TurbulentFlatData

FILE_PATH = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(os.path.dirname(FILE_PATH), "data")

# Mapping for DNS author labels
DNS_LABELS = {
    180: "DNS - Kim et al.",
    395: "DNS - Moser et al.",
    934: "DNS - Hoyas et al."
}


def plot_mean_velocity_profile(re_number, yplus_sim=None, Uplus_sim=None, 
                                label_sim="kynema-sgf", outpath=None):
    """
    Plot mean velocity profile with optional simulation data overlay.
    
    Parameters
    ----------
    re_number : int
        Reynolds number for experimental data (180, 395, 934)
    yplus_sim : np.ndarray, optional
        Simulation y+ values to overlay
    Uplus_sim : np.ndarray, optional
        Simulation U+ values to overlay
    label_sim : str
        Label for simulation data (default "kynema-sgf")
    outpath : Path or str, optional
        Path to save figure. If None, returns figure without saving.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    # Load experimental data
    data_handler = TurbulentFlatData(re_number)
    data = data_handler.get_mean_velocity()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot experimental data (every other point for clarity)
    y_plus_data = data.iloc[::2, 0]
    u_plus_data = data.iloc[::2, 1]
    
    dns_label = DNS_LABELS.get(re_number, f"DNS - Re{re_number}")
    ax.semilogx(y_plus_data, u_plus_data, 'ko', label=dns_label, markersize=8, 
                fillstyle='none', markeredgewidth=1)
    
    # Plot simulation data if provided
    if yplus_sim is not None and Uplus_sim is not None:
        ax.semilogx(yplus_sim, Uplus_sim, "-", label=label_sim, linewidth=2)
    
    # Add reference lines with appropriate ranges
    y_all = data.iloc[:, 0]
    y_viscous = y_all[y_all < 1.34e1]
    ax.semilogx(y_viscous, y_viscous, 'k--', alpha=0.6, linewidth=1.5, 
                label='$u^+ = y^+$ (viscous sublayer)')
    
    y_loglaw = y_all[y_all > 9]
    ax.semilogx(y_loglaw, 2.5 * np.log(y_loglaw) + 5.5, 'k:', alpha=0.6, linewidth=1.5,
                label='$u^+ = 2.5 \\ln(y^+) + 5.5$ (log law)')    
    
    ax.set_xlabel('$y^+$', fontsize=16)
    ax.set_ylabel('$u^+$', fontsize=16)
    ax.set_title(f'Mean Velocity Profile - Re$_\\tau$ = {re_number}', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    
    # Set y limits and x limits
    if re_number == 180:
        ax.set_ylim(0, 20)
    else:
        ax.set_ylim(0, 25)
    ax.set_xlim(1, max(200.0, y_all.max()))
    
    ax.legend(fontsize=12, loc='lower right')
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    
    if outpath is not None:
        fig.savefig(outpath, dpi=150)
    
    return fig


def plot_rms_velocity_profiles(re_number, yplus_sim=None, urms_sim=None, 
                                vrms_sim=None, wrms_sim=None, 
                                label_sim="kynema-sgf", outpath=None):
    """
    Plot RMS velocity profiles with optional simulation data overlay.
    
    Parameters
    ----------
    re_number : int
        Reynolds number for experimental data (180, 395, 934)
    yplus_sim : np.ndarray, optional
        Simulation y+ values to overlay
    urms_sim : np.ndarray, optional
        Simulation u_rms+ values to overlay
    vrms_sim : np.ndarray, optional
        Simulation v_rms+ values to overlay
    wrms_sim : np.ndarray, optional
        Simulation w_rms+ values to overlay
    label_sim : str
        Label for simulation data (default "kynema-sgf")
    outpath : Path or str, optional
        Path to save figure. If None, returns figure without saving.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    # Load experimental data
    data_handler = TurbulentFlatData(re_number)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Collect max y+ from experimental data for xlim calculation
    max_yplus = 0
    
    # Plot experimental data for all components and capture colors
    colors = {}
    for component in ['U', 'V', 'W']:
        data = data_handler.get_rms_velocity(component)
        max_yplus = max(max_yplus, data.iloc[:, 0].max())
        dns_label_comp = DNS_LABELS.get(re_number, f"DNS - Re{re_number}")
        line = ax.plot(data.iloc[:, 0], data.iloc[:, 1], 'o', markersize=8, 
                       markeredgewidth=1, fillstyle='none', label=f'{dns_label_comp} - ${component}$')
        colors[component] = line[0].get_color()
    
    # Plot simulation data if provided, using same colors
    if yplus_sim is not None:
        if urms_sim is not None:
            ax.plot(yplus_sim, urms_sim, "-", color=colors['U'], label=f'{label_sim} - $u$')
        if vrms_sim is not None:
            ax.plot(yplus_sim, vrms_sim, "-", color=colors['V'], label=f'{label_sim} - $v$')
        if wrms_sim is not None:
            ax.plot(yplus_sim, wrms_sim, "-", color=colors['W'], label=f'{label_sim} - $w$')
    
    ax.set_xlabel('$y^+$', fontsize=16)
    ax.set_ylabel(r'$u_{rms}$, $v_{rms}$, $w_{rms}$', fontsize=16)
    ax.set_title(f'RMS Velocity Profiles - Re$_\\tau$ = {re_number}', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_ylim(0, 3)
    ax.set_xlim(0, max_yplus * 1.02)
    ax.legend(fontsize=12, loc='upper right')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    
    if outpath is not None:
        fig.savefig(outpath, dpi=150)
    
    return fig


def main() -> int:
    """Standalone CLI interface for plotting experimental data."""
    parser = argparse.ArgumentParser(
        description="Plot data from single phase turbulent flat validation case."
    )
    parser.add_argument("--Re", type=int, default=180, 
                       help="Reynolds number (180, 395, or 934)")
    parser.add_argument("--plot-type", type=str, default="both", 
                       choices=["mean-velocity", "rms-velocities", "both"],
                       help="Which plots to generate")
    args = parser.parse_args()
    
    # Generate plots
    if args.plot_type in ["mean-velocity", "both"]:
        fig_mean = plot_mean_velocity_profile(args.Re)
        plt.show()
    
    if args.plot_type in ["rms-velocities", "both"]:
        fig_rms = plot_rms_velocity_profiles(args.Re)
        plt.show()
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

