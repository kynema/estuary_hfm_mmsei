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
    python post_process_sampling.py --Re 180 --LES --all
    python post_process_sampling.py --Re 180 --LES --all --drag none,og,tf1 --utau-source gradU,gradP,gradP
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
        Drag model: 'og' or 'tf1' (only used with IB, default tf1)
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


def process_case(Re: int, DNS: bool, LES: bool, IB: bool, drag: str, 
                 utau_source: str, nprocs: int = 1) -> tuple:
    """
    Load particle statistics for a single case, compute u_tau via the given
    method, and rescale to wall units.
    
    Parameters
    ----------
    Re : int
        Stress Reynolds number
    DNS : bool
        DNS variant
    LES : bool
        LES variant
    IB : bool
        Immersed boundary variant
    drag : str or None
        Drag model: 'og' or 'tf1' (only used with IB)
    utau_source : str
        'gradP', 'gradU', or 'loglaw'
    nprocs : int
        Number of parallel workers passed to process_stats.py
    
    Returns
    -------
    tuple
        (config, rescaled) where config is the domain_and_flow() dict (with
        u_tau updated) and rescaled is the dict returned by rescale_to_wall_units()
    """
    config = domain_and_flow(Re=Re, IB=IB)
    
    stats_file = find_particle_stats_file(Re, DNS=DNS, LES=LES, IB=IB, drag=drag,
                                          nprocs=nprocs)
    print(f"Loading particle statistics from: {stats_file}")
    
    z_coord, u_mean, u_rms, v_sgf_mean, v_sgf_rms, w_sgf_mean, w_sgf_rms = load_particle_stats(stats_file)
    
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
    
    rescaled = rescale_to_wall_units(
        z_coord, u_mean, u_rms, v_sgf_rms, w_sgf_rms,
        config, Re, DNS=DNS, LES=LES, IB=IB, drag=drag
    )
    
    return config, rescaled


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


def run_all(Re: int, drag_list: list, utau_list: list, nprocs: int = 1) -> int:
    """
    Generate comparison plots for LES across drag models.
    
    Parameters
    ----------
    Re : int
        Stress Reynolds number
    drag_list : list of str
        Drag models to compare, in order. 'none' means no IB; otherwise 'og' or 'tf1'.
    utau_list : list of str
        Friction velocity source ('gradP', 'gradU', or 'loglaw') to use for each
        entry in drag_list, matched by position.
    nprocs : int
        Number of parallel workers passed to process_stats.py
    
    Returns
    -------
    int
        0 on success (with at least one case processed), 1 if no cases could be processed
    """
    cases = []
    for drag, utau_source in zip(drag_list, utau_list):
        if drag == 'none':
            cases.append({'IB': False, 'drag': None, 'utau_source': utau_source, 'label': 'kynema-sgf (no IB)'})
        else:
            cases.append({'IB': True, 'drag': drag, 'utau_source': utau_source, 'label': f'kynema-sgf (IB: {drag.upper()})'})
    
    mean_sim_list = []
    rms_sim_list = []
    
    for case in cases:
        try:
            config, rescaled = process_case(Re, DNS=False, LES=True, IB=case['IB'],
                                            drag=case['drag'], utau_source=case['utau_source'],
                                            nprocs=nprocs)
        except FileNotFoundError as e:
            print(f"Warning: skipping '{case['label']}' - {e}", file=sys.stderr)
            continue
        
        output_dir = get_output_dir(Re, LES=True, IB=case['IB'], drag=case['drag'])
        save_profiles(output_dir, rescaled)
        
        mean_sim_list.append({
            'yplus': rescaled['y_plus'],
            'Uplus': rescaled['U_plus'],
            'label': case['label'],
        })
        rms_sim_list.append({
            'yplus': rescaled['y_plus'],
            'urms': rescaled['urms_plus'],
            'vrms': rescaled['vrms_plus'],
            'wrms': rescaled['wrms_plus'],
            'label': case['label'],
        })
    
    if not mean_sim_list:
        print("Error: no cases could be processed for --all", file=sys.stderr)
        return 1
    
    output_dir = Path(__file__).parent.parent / "figures" / f"ReTau{Re}_LES_ALL"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Plotting mean velocity profile comparison...")
    plot_file = output_dir / "Uplus_sampling_all.png"
    plot_mean_velocity_profile(Re, sim_list=mean_sim_list, outpath=str(plot_file))
    print(f"Saved plot to {plot_file}")
    
    print("Plotting RMS velocity profile comparison...")
    plot_file = output_dir / "VelRMSplus_sampling_all.png"
    plot_rms_velocity_profiles(Re, sim_list=rms_sim_list, outpath=str(plot_file))
    print(f"Saved plot to {plot_file}")
    
    print("\nSuccessfully processed all LES drag model comparisons!")
    return 0


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
    parser.add_argument("--drag", type=str, default=None,
                        help="Drag model. Single-case mode: 'og' or 'tf1', only with --IB "
                             "(default: 'tf1' when --IB is set). With --all: comma-separated list "
                             "combining 'none' (no IB), 'og', 'tf1', e.g. 'none,og,tf1' "
                             "(default: 'none,og,tf1').")
    parser.add_argument("--all", action="store_true",
                        help="Generate a comparison plot across LES drag models for the given --Re. "
                             "Requires --LES; incompatible with --DNS and --IB. Use --drag and "
                             "--utau-source (comma-separated, matched by position) to customize which "
                             "models are compared and which u_tau source each uses.")
    parser.add_argument("--nprocs", type=int,
                        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)),
                        help="Number of parallel workers for reading binary data "
                             "(default: $SLURM_CPUS_PER_TASK or 1)")
    parser.add_argument("--utau-source", type=str, default=None,
                        help="Source for friction velocity. Single-case mode: 'gradU' (quadratic wall "
                             "gradient, default for non-IB), 'loglaw' (log-law from wall model, default "
                             "for IB), or 'gradP' (pressure gradient). With --all: comma-separated list "
                             "matching --drag, e.g. 'gradU,gradP,gradP' (default: 'gradU,gradP,gradP').")
    
    args = parser.parse_args()
    
    VALID_UTAU = {'gradP', 'gradU', 'loglaw'}
    VALID_DRAG_ALL = {'none', 'og', 'tf1'}
    VALID_DRAG_SINGLE = {'og', 'tf1'}
    
    # Validate options
    if args.all:
        if not args.LES:
            parser.error("--all requires --LES")
        if args.DNS or args.IB:
            parser.error("--all cannot be combined with --DNS or --IB")
        
        drag_list = args.drag.split(',') if args.drag else ['none', 'og', 'tf1']
        utau_list = args.utau_source.split(',') if args.utau_source else ['gradU', 'gradP', 'gradP']
        
        if len(drag_list) != len(utau_list):
            parser.error(f"--drag and --utau-source must have the same number of comma-separated "
                        f"values ({len(drag_list)} vs {len(utau_list)})")
        for d in drag_list:
            if d not in VALID_DRAG_ALL:
                parser.error(f"Invalid --drag value '{d}' for --all; choices are: none, og, tf1")
        for u in utau_list:
            if u not in VALID_UTAU:
                parser.error(f"Invalid --utau-source value '{u}'; choices are: gradP, gradU, loglaw")
        
        return run_all(args.Re, drag_list, utau_list, args.nprocs)
    
    if not args.DNS and not args.LES:
        parser.error("Either --DNS or --LES must be specified")
    
    if args.DNS and args.LES:
        parser.error("Cannot specify both --DNS and --LES")
    
    if args.IB and not args.LES:
        parser.error("--IB can only be used with --LES")
    
    if args.drag is not None and args.drag not in VALID_DRAG_SINGLE:
        parser.error(f"--drag must be 'og' or 'tf1' (got '{args.drag}')")
    if args.IB and not args.drag:
        args.drag = 'tf1'
    if args.drag and not args.IB:
        parser.error("--drag can only be used with --IB")
    
    if args.utau_source is not None and args.utau_source not in VALID_UTAU:
        parser.error(f"--utau-source must be one of gradP, gradU, loglaw (got '{args.utau_source}')")
    
    # Default utau_source based on case type:
    #   DNS:  gradU (quadratic fit, valid for resolved no-slip wall)
    #   LES:  loglaw (Monin-Obukhov log-law, consistent with wall model)
    if args.utau_source is not None:
        utau_source = args.utau_source
    elif args.DNS:
        utau_source = 'gradU'
    elif args.Re < 190.0:  # LES (with or without IB)
        utau_source = 'gradU'
    else:
        utau_source = 'loglaw'
    args.utau_source = utau_source
    
    try:
        config, rescaled = process_case(args.Re, DNS=args.DNS, LES=args.LES, IB=args.IB,
                                        drag=args.drag, utau_source=utau_source,
                                        nprocs=args.nprocs)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
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
