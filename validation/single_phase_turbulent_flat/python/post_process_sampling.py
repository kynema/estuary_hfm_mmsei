"""
Post-process particle sampling statistics and plot against experimental data.

This script reads particle sampling statistics (from process_stats.py output)
and computes wall-normal velocity profiles for comparison with DNS data.
Does not require yt; all data is read directly from particle_stats.txt.

Usage:
    python post_process_sampling.py --Re 180 --DNS
    python post_process_sampling.py --Re 180 --LES
    python post_process_sampling.py --Re 180 --LES --IB --drag og
    python post_process_sampling.py --Re 180 --LES --IB --drag tf1
"""

import numpy as np
import pandas as pd
import argparse
import sys
import os
import subprocess
from pathlib import Path

from case_setup import domain_and_flow, get_pressure_gradient
from plot_data import plot_mean_velocity_profile, plot_rms_velocity_profiles
from data import build_case_dir_name, get_output_dir


def find_particle_stats_file(Re: int, DNS: bool = False, LES: bool = False,
                             IB: bool = False, drag: str = None,
                             nprocs: int = 1) -> Path:
    """
    Find the particle_stats.txt file for the given Reynolds number.
    If it doesn't exist, run process_stats.py to generate it.
    
    Parameters
    ----------
    Re : int
        Stress Reynolds number
    DNS : bool
        If True, look for DNS variant
    LES : bool
        If True, look for LES variant
    IB : bool
        If True, look for immersed boundary variant
    drag : str
        Drag model: 'og' or 'tf1' (only used with IB)
    nprocs : int
        Number of parallel workers passed to process_stats.py
    
    Returns
    -------
    Path
        Path to particle_stats.txt
    """
    case_dir = Path(__file__).parent.parent / "cases"
    case_name = build_case_dir_name(Re, DNS=DNS, LES=LES, IB=IB, drag=drag)
    stats_file = case_dir / case_name / "particle_stats.txt"
    
    # If particle_stats.txt doesn't exist, generate it
    if not stats_file.exists():
        print(f"Particle statistics file not found. Generating from sampling data...")
        
        # Build command to run process_stats.py
        python_dir = Path(__file__).parent
        cmd = ["python", str(python_dir / "process_stats.py"), "--Re", str(Re)]
        if DNS:
            cmd.append("--DNS")
        if LES:
            cmd.append("--LES")
        if IB:
            cmd.append("--IB")
        if drag:
            cmd.extend(["--drag", drag])
        if nprocs > 1:
            cmd.extend(["--nprocs", str(nprocs)])
        
        try:
            result = subprocess.run(cmd, cwd=str(python_dir), 
                                   capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"process_stats.py failed:\n{result.stderr}")
            print(result.stdout)
        except Exception as e:
            raise RuntimeError(f"Failed to generate particle statistics: {e}")
        
        # Verify the file was created
        if not stats_file.exists():
            raise FileNotFoundError(f"Particle statistics file still not found after running process_stats.py: {stats_file}")
    
    return stats_file


def load_particle_stats(stats_file: Path) -> tuple:
    """
    Load particle statistics from file.
    
    Parameters
    ----------
    stats_file : Path
        Path to particle_stats.txt
    
    Returns
    -------
    tuple
        (z_coord, u_mean, u_var, v_sgf_mean, v_sgf_var, w_sgf_mean, w_sgf_var)
    """
    # Read the file, skipping the header comment
    data = np.loadtxt(stats_file, skiprows=1)
    
    # Extract columns
    z_coord = data[:, 0]
    u_mean = data[:, 1]
    u_var = data[:, 2]
    v_sgf_mean = data[:, 3]
    v_sgf_var = data[:, 4]
    w_sgf_mean = data[:, 5]
    w_sgf_var = data[:, 6]
    
    # Compute RMS (standard deviation)
    u_rms = np.sqrt(u_var)
    v_sgf_rms = np.sqrt(v_sgf_var)
    w_sgf_rms = np.sqrt(w_sgf_var)
    
    return z_coord, u_mean, u_rms, v_sgf_mean, v_sgf_rms, w_sgf_mean, w_sgf_rms


def rescale_to_wall_units(z_coord, u_mean, u_rms, v_sgf_rms, w_sgf_rms, 
                         config, Re, DNS=False, LES=False, IB=False, drag=None):
    """
    Map z-coordinates to wall units and rescale velocities.
    
    Parameters
    ----------
    z_coord : array
        z-coordinates in domain coordinates [-delta, +delta] where -delta=bottom wall, 0=centerline
    u_mean : array
        Mean streamwise velocity
    u_rms : array
        Streamwise velocity RMS
    v_sgf_rms : array
        Spanwise velocity RMS (in kynema-sgf coordinates)
    w_sgf_rms : array
        Wall-normal velocity RMS (in kynema-sgf coordinates)
    config : dict
        Configuration dictionary from domain_and_flow()
    Re : int
        Stress Reynolds number for metadata
    DNS : bool
        DNS variant
    LES : bool
        LES variant
    IB : bool
        Immersed boundary variant
    drag : str
        Drag model: 'og' or 'tf1'
    
    Returns
    -------
    dict
        Dictionary with wall-unit quantities
    """
    delta = config['delta']
    u_tau = config['u_tau']
    l_nu = config['mu'] / (config['density'] * u_tau)
    
    # z_coord is in domain coordinates [-delta, +delta]
    # Convert to wall-normal distance from bottom wall: y = z + delta
    y_phys = z_coord + delta
    print(f"y_phys (wall-normal distance from bottom wall): {y_phys}")
    
    # Convert to wall units: y+ = y / l_nu
    y_plus = y_phys / l_nu
    
    # Rescale velocities to wall units
    U_plus = u_mean / u_tau
    urms_plus = u_rms / u_tau
    
    # Coordinate system swap for DNS interpretation:
    # Kynema-sgf uses z as wall-normal, but DNS data uses y as wall-normal
    # In kynema-sgf: v_sgf = spanwise, w_sgf = wall-normal
    # In DNS: V = spanwise, W = wall-normal
    # So: vrms_sim = v_sgf_rms (spanwise), wrms_sim = w_sgf_rms (wall-normal)
    vrms_plus = v_sgf_rms / u_tau  # spanwise RMS
    wrms_plus = w_sgf_rms / u_tau  # wall-normal RMS
    
    return {
        'y_phys': y_phys,
        'y_plus': y_plus,
        'U_plus': U_plus,
        'urms_plus': urms_plus,
        'vrms_plus': vrms_plus,
        'wrms_plus': wrms_plus,
    }


def save_profiles(output_dir: Path, rescaled_data: dict):
    """
    Save velocity profiles to CSV file.
    
    Parameters
    ----------
    output_dir : Path
        Directory to save to
    rescaled_data : dict
        Dictionary with y_plus, U_plus, urms_plus, vrms_plus, wrms_plus
    """
    df = pd.DataFrame({
        'z_phys': rescaled_data['y_phys'],
        'y_plus': rescaled_data['y_plus'],
        'U_plus': rescaled_data['U_plus'],
        'urms_plus': rescaled_data['urms_plus'],
        'vrms_plus': rescaled_data['vrms_plus'],
        'wrms_plus': rescaled_data['wrms_plus'],
    })
    
    output_file = output_dir / "profiles_sampling.csv"
    df.to_csv(output_file, index=False)
    print(f"Saved profiles to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Post-process particle sampling statistics and plot against DNS data"
    )
    parser.add_argument("--Re", type=int, default=180,
                        help="Stress Reynolds number (default: 180)")
    parser.add_argument("--DNS", action="store_true",
                        help="Use DNS variant")
    parser.add_argument("--LES", action="store_true",
                        help="Use LES variant")
    parser.add_argument("--IB", action="store_true",
                        help="Use immersed boundary variant")
    parser.add_argument("--drag", type=str, choices=['og', 'tf1'],
                        help="Drag model for IB cases: 'og' or 'tf1'")
    parser.add_argument("--nprocs", type=int,
                        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)),
                        help="Number of parallel workers for reading binary data "
                             "(default: $SLURM_CPUS_PER_TASK or 1)")
    parser.add_argument("--utau-source", choices=['gradP', 'gradU', 'loglaw'], default=None,
                        help="Source for friction velocity: 'gradU' (quadratic wall gradient, default for non-IB), "
                             "'loglaw' (log-law from wall model using cell k+1, default for IB), "
                             "or 'gradP' (pressure gradient)")
    
    args = parser.parse_args()
    
    # Validate options
    if not args.DNS and not args.LES:
        parser.error("Either --DNS or --LES must be specified")
    
    if args.DNS and args.LES:
        parser.error("Cannot specify both --DNS and --LES")
    
    if args.IB and not args.LES:
        parser.error("--IB can only be used with --LES")
    
    if args.IB and not args.drag:
        parser.error("--drag option required when using --IB")
    
    if args.drag and not args.IB:
        parser.error("--drag can only be used with --IB")
    
    # Load configuration
    config = domain_and_flow(Re=args.Re, IB=args.IB)
    
    # Find and load particle statistics
    try:
        stats_file = find_particle_stats_file(args.Re, DNS=args.DNS, LES=args.LES, 
                                             IB=args.IB, drag=args.drag,
                                             nprocs=args.nprocs)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    print(f"Loading particle statistics from: {stats_file}")
    
    try:
        z_coord, u_mean, u_rms, v_sgf_mean, v_sgf_rms, w_sgf_mean, w_sgf_rms = load_particle_stats(stats_file)
    except Exception as e:
        print(f"Error loading particle statistics: {e}", file=sys.stderr)
        return 1
    
    # Default utau_source based on case type:
    #   DNS:  gradU (quadratic fit, valid for resolved no-slip wall)
    #   LES:  loglaw (Monin-Obukhov log-law, consistent with wall model)
    if args.utau_source is not None:
        utau_source = args.utau_source
    elif args.DNS:
        utau_source = 'gradU'
        args.utau_source = utau_source
    elif args.Re < 190.0:  # LES (with or without IB)
        utau_source = 'gradU'
        args.utau_source = utau_source
    else:
        utau_source = 'loglaw'
        args.utau_source = utau_source

    # Calculate u_tau
    if utau_source == 'gradP':
        pass  # config['u_tau'] already set from pressure gradient
    else:
        mu = config['mu']
        density = config['density']
        yw = z_coord + config['delta']
        order = np.argsort(yw)
        j0, j1, j2 = order[0], order[1], order[2]
        if utau_source == 'loglaw':
            # Log-law consistent with IB og wall model (PDF section 3, eq. 22):
            # u* = M_{k+1} * kappa / ln(y_{k+1}/z0)
            # j0,j1 are the two fine-grid sub-cells of the base-mesh drag cell;
            # j2 is the first fine-grid cell of k+1 (first cell above drag cell).
            kappa = 0.41
            z0 = 1e-5  # ABL.surface_roughness_z0
            config['u_tau'] = abs(u_mean[j2]) * kappa / np.log(yw[j2] / z0)
        else:  # gradU
            print("\nCalculating u_tau from quadratic fit to wall gradient...")
            y0, y1 = yw[j0], yw[j1]
            U0, U1 = u_mean[j0], u_mean[j1]
            print(f"   u_tau from gradP: {config['u_tau']:.6e}")
            print(f"   Using points: (y0={y0:.6e}, U0={U0:.6e}), (y1={y1:.6e}, U1={U1:.6e})")
            # Quadratic fit U(y) = a*y + b*y^2 through (0,0), (y0,U0), (y1,U1)
            # dU/dy|_wall = a = (U0*y1^2 - U1*y0^2) / (y0*y1*(y1 - y0))
            dU_dy_wall = (U0 * y1**2 - U1 * y0**2) / (y0 * y1 * (y1 - y0))
            tau_w = mu * abs(dU_dy_wall)
            config['u_tau'] = np.sqrt(tau_w / density)
            print(f"   u_tau from gradU: {config['u_tau']:.6e}")

    # Rescale to wall units
    rescaled = rescale_to_wall_units(
        z_coord, u_mean, u_rms, v_sgf_rms, w_sgf_rms,
        config, args.Re, DNS=args.DNS, LES=args.LES, IB=args.IB, drag=args.drag
    )
    
    # Get output directory
    output_dir = get_output_dir(args.Re, DNS=args.DNS, LES=args.LES, 
                               IB=args.IB, drag=args.drag)
    
    # Save profiles
    save_profiles(output_dir, rescaled)
    
    # Plot mean velocity profile
    print("Plotting mean velocity profile...")
    plot_file = output_dir / f"Uplus_sampling_{args.utau_source}.png"
    plot_mean_velocity_profile(
        args.Re,
        yplus_sim=rescaled['y_plus'],
        Uplus_sim=rescaled['U_plus'],
        label_sim="kynema-sgf",
        outpath=str(plot_file)
    )
    print(f"Saved plot to {plot_file}")
    
    # Plot RMS velocity profiles
    print("Plotting RMS velocity profiles...")

    plot_file = output_dir / f"VelRMSplus_sampling_{args.utau_source}.png"
    plot_rms_velocity_profiles(
        args.Re,
        yplus_sim=rescaled['y_plus'],
        urms_sim=rescaled['urms_plus'],
        vrms_sim=rescaled['vrms_plus'],
        wrms_sim=rescaled['wrms_plus'],
        label_sim="kynema-sgf",
        outpath=str(plot_file)
    )
    print(f"Saved plot to {plot_file}")
    
    print("\nSuccessfully processed particle sampling statistics!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
