"""
Post-process pre-averaged pltAvg files (velocity_mean_averaging*) and plot against DNS data.

This script reads the LAST pltAvg file from a case, which contains
time-averaged velocities (velocity_mean_averagingx, etc.) computed
directly by kynema-sgf during simulation.
"""

from __future__ import annotations

import argparse
import glob
import math
import yt 
import sys
from pathlib import Path
import numpy as np

import case_setup
from plot_data import plot_mean_velocity_profile, plot_rms_velocity_profiles
from data import build_case_dir_name, get_output_dir

# --------------------------------------------------------------------
# Plotfile loading
# --------------------------------------------------------------------

def find_last_pltavg_file(Re: int, DNS: bool = False, LES: bool = False,
                          IB: bool = False, drag: str = None) -> Path:
    """
    Find the LAST pltAvg file for the given Reynolds number and variant.
    
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
    
    Returns
    -------
    Path
        Path to the LAST pltAvg Header file
    """
    # Construct case directory path
    case_dir = Path(__file__).parent.parent / "cases"
    case_name = build_case_dir_name(Re, DNS=DNS, LES=LES, IB=IB, drag=drag)
    pltavg_pattern = case_dir / case_name / "pltAvg*" / "Header"
    
    files = sorted(glob.glob(str(pltavg_pattern)))
    
    if not files:
        raise FileNotFoundError(
            f"No pltAvg files found in {case_dir / case_name}/"
        )
    
    # Return the LAST file (latest time step)
    return Path(files[-1])


def load_pltavg_averaged_velocity(plt_path: Path, delta: float, IB: bool = False):
    """
    Load pre-averaged velocities and Reynolds stresses from a kynema-sgf pltAvg file via yt.
    
    The file should contain velocity_mean_averagingx, velocity_mean_averagingy,
    velocity_mean_averagingz (time-averaged means) and velocity_reynolds_stress_averaging
    components (time-averaged Reynolds stresses: <u'u'>, <v'v'>, <w'w'>).
    
    For IB cases, masks solid cells outside channel bounds [-delta/2, +delta/2].
    
    Parameters
    ----------
    plt_path : Path
        Path to pltAvg directory or its Header file
    delta : float
        Channel half-height in meters
    IB : bool
        If True, use z-location to mask solid cells
    
    Returns
    -------
    u, v, w : np.ndarray
        Time-averaged velocity components (shape Nx, Ny, Nz)
    uu, vv, ww : np.ndarray
        Time-averaged normal Reynolds stress components <u'u'>, <v'v'>, <w'w'> (shape Nx, Ny, Nz)
    prob_lo, prob_hi : tuple of length 3
        Domain bounds
    time : float
        Simulation time
    """

    yt.set_log_level("error")
    
    # Allow user to pass the pltAvg dir or the Header inside it
    p = Path(plt_path)
    if p.is_file() and p.name == "Header":
        p = p.parent
    
    ds = yt.load(str(p))
    
    # Get dimensions at the finest AMR level
    max_level = ds.max_level
    ref_ratio = 2  # AMRex default refinement ratio
    fine_dims = ds.domain_dimensions * (ref_ratio ** max_level)
    
    # Create covering grid at finest level resolution
    cg = ds.covering_grid(level=max_level,
                          left_edge=ds.domain_left_edge,
                          dims=fine_dims)

    # Load time-averaged velocity components
    u = np.asarray(cg["boxlib", "velocity_mean_averagingx"])
    v = np.asarray(cg["boxlib", "velocity_mean_averagingy"])
    w = np.asarray(cg["boxlib", "velocity_mean_averagingz"])
    
    # Load time-averaged Reynolds stress components (normal stresses)
    # velocity_reynolds_stress_averaging0 = <u'u'>
    # velocity_reynolds_stress_averaging3 = <v'v'>
    # velocity_reynolds_stress_averaging5 = <w'w'>
    uu = np.asarray(cg["boxlib", "velocity_reynolds_stress_averaging0"])
    vv = np.asarray(cg["boxlib", "velocity_reynolds_stress_averaging3"])
    ww = np.asarray(cg["boxlib", "velocity_reynolds_stress_averaging5"])
    
    # For IB cases, mask out solid cells based on z-location
    # The channel extends from -delta/2 to +delta/2 (centered at z=0)
    # IB solid regions are outside these bounds
    if IB:
        z_coords = np.asarray(cg["boxlib", "z"])
        z_min = -delta / 2.0
        z_max = delta / 2.0
        mask_solid = (z_coords < z_min) | (z_coords > z_max)
        u[mask_solid] = np.nan
        v[mask_solid] = np.nan
        w[mask_solid] = np.nan
        uu[mask_solid] = np.nan
        vv[mask_solid] = np.nan
        ww[mask_solid] = np.nan
        n_masked = np.sum(mask_solid)
        n_total = mask_solid.size
        print(f"  Masked {n_masked}/{n_total} cells outside channel bounds [{z_min:.6f}, {z_max:.6f}]", file=sys.stderr)
    
    prob_lo = tuple(float(x) for x in ds.domain_left_edge)
    prob_hi = tuple(float(x) for x in ds.domain_right_edge)
    time = float(ds.current_time)
    
    return u, v, w, uu, vv, ww, prob_lo, prob_hi, time

# --------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------

def wall_normal_stats(u, v, w, uu, vv, ww, periodic_axes=(0, 1), IB: bool = False):
    """
    Spatially average over periodic axes (x, y by default) to get
    per-z profiles of mean and RMS.
    
    For IB cases, uses nanmean to ignore NaN cells.
    
    Parameters
    ----------
    u, v, w : np.ndarray
        Time-averaged mean velocity components
    uu, vv, ww : np.ndarray
        Time-averaged Reynolds stress normal components (<u'u'>, <v'v'>, <w'w'>)
        pre-computed by kynema-sgf
    periodic_axes : tuple
        Axes over which to average spatially
    IB : bool
        If True, use nanmean to skip solid cells
    
    Returns
    -------
    dict
        Dicts of arrays indexed by wall-normal index k:
        Keys: U_bar, V_bar, W_bar, U_rms, V_rms, W_rms
        Values: 1D arrays (length Nz for non-IB, slightly less for IB due to NaN masking)
    """
    if IB:
        # Use nanmean to skip solid cells
        U_bar = np.nanmean(u, axis=periodic_axes)
        V_bar = np.nanmean(v, axis=periodic_axes)
        W_bar = np.nanmean(w, axis=periodic_axes)
        # RMS = sqrt(<u'u'>)
        U_rms_prof = np.sqrt(np.nanmean(uu, axis=periodic_axes))
        V_rms_prof = np.sqrt(np.nanmean(vv, axis=periodic_axes))
        W_rms_prof = np.sqrt(np.nanmean(ww, axis=periodic_axes))
    else:
        U_bar = u.mean(axis=periodic_axes)
        V_bar = v.mean(axis=periodic_axes)
        W_bar = w.mean(axis=periodic_axes)
        # RMS = sqrt(<u'u'>)
        U_rms_prof = np.sqrt(uu.mean(axis=periodic_axes))
        V_rms_prof = np.sqrt(vv.mean(axis=periodic_axes))
        W_rms_prof = np.sqrt(ww.mean(axis=periodic_axes))
    
    return {
        "U_bar": U_bar, "V_bar": V_bar, "W_bar": W_bar,
        "U_rms": U_rms_prof, "V_rms": V_rms_prof, "W_rms": W_rms_prof
    }


# --------------------------------------------------------------------
# Main and plotting
# --------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--Re", type=int, default=180,
                    help="Stress Reynolds number: 180, 395, or 934 "
                         "(default: 180)")
    ap.add_argument("--DNS", action="store_true",
                    help="Use DNS variant")
    ap.add_argument("--LES", action="store_true",
                    help="Use LES variant")
    ap.add_argument("--IB", action="store_true",
                    help="Use immersed boundary variant")
    ap.add_argument("--drag", type=str, choices=['og', 'tf1'],
                    help="Drag model for IB cases: 'og' or 'tf1'")
    args = ap.parse_args()

    # Validate options
    if not args.DNS and not args.LES:
        ap.error("Either --DNS or --LES must be specified")
    
    if args.DNS and args.LES:
        ap.error("Cannot specify both --DNS and --LES")
    
    if args.IB and not args.LES:
        ap.error("--IB can only be used with --LES")
    
    if args.IB and not args.drag:
        ap.error("--drag option required when using --IB")
    
    if args.drag and not args.IB:
        ap.error("--drag can only be used with --IB")

    # Find the last pltAvg file and set output dir
    try:
        pltavg_file = find_last_pltavg_file(args.Re, DNS=args.DNS, LES=args.LES,
                                            IB=args.IB, drag=args.drag)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    outdir = get_output_dir(args.Re, DNS=args.DNS, LES=args.LES,
                           IB=args.IB, drag=args.drag)

    variant = ""
    if args.DNS:
        variant = "DNS"
    elif args.LES:
        variant = "LES"
        if args.IB:
            variant += f" + IB ({args.drag.upper()})"
    
    print(f"Post-processing Re_tau = {args.Re}, variant = {variant}")
    print(f"Using pre-averaged pltAvg file: {pltavg_file.parent.name}")
    print(f"Output directory: {outdir}\n")

    # Get flow parameters from case_setup
    config = case_setup.domain_and_flow(Re=args.Re, IB=args.IB)
    delta = config['delta']
    density = config['density']
    mu = config['mu']
    dpdx_magnitude = abs(config['dpdx'])  # case_setup stores as signed; we need magnitude
    
    print(f"Configuration:")
    print(f"  delta = {delta:.6f} m")
    print(f"  density = {density:.6f} kg/m^3")
    print(f"  mu = {mu:.6e} Pa·s")
    print(f"  |dP/dx| = {dpdx_magnitude:.2f} Pa/m\n")

    # Load the pre-averaged velocity file
    u, v, w, uu, vv, ww, prob_lo, prob_hi, t = load_pltavg_averaged_velocity(pltavg_file, 
                                                                                delta=delta, 
                                                                                IB=args.IB)
    
    print(f"  loaded {pltavg_file.parent.name}  t={t:.4e}  shape={u.shape}")
    
    stats = wall_normal_stats(u, v, w, uu, vv, ww, periodic_axes=(0, 1), IB=args.IB)
    print(f"  computed wall-normal statistics\n")

    Nz = stats["U_bar"].size
    shape = u.shape

    # Wall-normal physical coordinate (z)
    # In kynema-sgf: walls at z = ±delta (prob_lo[2] and prob_hi[2])
    # Measure distance from lower wall only (to centerline)
    k = np.arange(Nz)
    dz = (prob_hi[2] - prob_lo[2]) / shape[2]
    z_cc = prob_lo[2] + (k + 0.5) * dz  # cell-center coordinates
    
    # Wall distance: distance from lower wall
    zw = z_cc - prob_lo[2]
    
    # Keep only lower wall to centerline (zw <= delta)
    mask = zw <= delta
    zw = zw[mask]
    U_bar_filt = stats["U_bar"][mask]
    V_bar_filt = stats["V_bar"][mask]
    W_bar_filt = stats["W_bar"][mask]
    U_rms_filt = stats["U_rms"][mask]
    V_rms_filt = stats["V_rms"][mask]
    W_rms_filt = stats["W_rms"][mask]

    # Friction velocity (from pressure gradient)
    u_tau = math.sqrt(dpdx_magnitude * delta / density)
    
    nu = mu / density
    print(f"Results:")
    print(f"  u_tau = {u_tau:.5g}")
    print(f"  nu    = {nu:.5g}")
    print(f"  delta = {delta:.5g}")
    print(f"  Re_tau = {u_tau * delta / nu:.3f}\n")

    # Wall-units rescale
    # RMS is computed from time-averaged Reynolds stress: RMS = sqrt(<u'u'>)
    yplus = zw * u_tau / nu  # y+ in data coordinate frame
    Uplus = U_bar_filt / u_tau
    urms_p = U_rms_filt / u_tau  # sqrt(<u'u'>) / u_tau
    vrms_p = V_rms_filt / u_tau
    wrms_p = W_rms_filt / u_tau

    # Sort by y+ for clean plotting
    order = np.argsort(yplus)
    yplus_o = yplus[order]
    Uplus_o = Uplus[order]
    urms_o = urms_p[order]
    vrms_o = vrms_p[order]
    wrms_o = wrms_p[order]

    # Output
    outdir.mkdir(parents=True, exist_ok=True)

    # CSV file with profiles
    csv_path = outdir / "profiles_avg.csv"
    with csv_path.open("w") as fh:
        fh.write("# y_plus, U_plus, urms_plus, vrms_plus, wrms_plus\n")
        fh.write(f"# u_tau = {u_tau:.6g}  nu = {nu:.6g}  "
                 f"Re_tau = {u_tau * delta / nu:.3f}\n")
        for k in range(len(yplus_o)):
            fh.write(f"{yplus_o[k]:.6e}, {Uplus_o[k]:.6e}, "
                     f"{urms_o[k]:.6e}, {vrms_o[k]:.6e}, {wrms_o[k]:.6e}\n")
    print(f"Output files:")
    print(f"  wrote {csv_path}")

    # Plot with experimental data overlay
    plot_mean_velocity_profile(
        args.Re,
        yplus_sim=yplus_o,
        Uplus_sim=Uplus_o,
        label_sim="kynema-sgf",
        outpath=outdir / "Uplus_avg.png"
    )
    print(f"  wrote {outdir / 'Uplus_avg.png'}")

    # Plot RMS velocity profiles with simulation data overlay
    plot_rms_velocity_profiles(
        args.Re,
        yplus_sim=yplus_o,
        urms_sim=urms_o,
        vrms_sim=vrms_o,  # V_rms = sqrt(<v'v'>) from velocity_reynolds_stress_averaging3
        wrms_sim=wrms_o,  # W_rms = sqrt(<w'w'>) from velocity_reynolds_stress_averaging5
        label_sim="kynema-sgf",
        outpath=outdir / "VelRMSplus_avg.png"
    )
    print(f"  wrote {outdir / 'VelRMSplus_avg.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
