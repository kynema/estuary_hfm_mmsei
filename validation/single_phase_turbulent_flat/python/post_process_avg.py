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
    
    For IB cases, masks solid cells outside channel bounds [-delta, +delta].
    
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
    z_coords_1d : np.ndarray or None
        1D array of z cell center coordinates (for IB cases only, else None)
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

    # Load z-coordinates
    z_coords = np.asarray(cg["boxlib", "z"])

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
    
    # For IB cases, mask out solid cells based on terrain_blank
    # terrain_blank < 1 means fluid, >= 1 means solid (IB)
    z_coords_1d = None  # Will store 1D array of z cell centers for IB

    print(f"DEBUG: z before masking: shape={z_coords.shape}, z range=[{np.min(z_coords):.6f}, {np.max(z_coords):.6f}]", file=sys.stderr)
    if IB:
        terrain_blank = np.asarray(cg["boxlib", "terrain_blank"])
        mask = terrain_blank < 1  # True for fluid, False for solid
        z_coords[~mask] = np.nan

        # Reshape to remove nans in z-direction
        valid_z_mask = np.any(np.isnan(z_coords) == False, axis=(0, 1))

        # Get the index of the first and last valid z-slices
        valid_z_indices = np.where(valid_z_mask)[0]
        first_valid_z = valid_z_indices[0]
        last_valid_z = valid_z_indices[-1]

        # Create a dynamic slice for the z-axis (axis 2)
        # We add +1 to last_valid_z because Python slicing is exclusive at the stop index
        z_slice = slice(first_valid_z, last_valid_z + 1)

        # Slice all of your 3D arrays to get the trimmed shape
        u        = u[:, :, z_slice]
        v        = v[:, :, z_slice]
        w        = w[:, :, z_slice]
        uu       = uu[:, :, z_slice]
        vv       = vv[:, :, z_slice]
        ww       = ww[:, :, z_slice]
        z_coords = z_coords[:, :, z_slice]

        print(f"Trimmed Shape: {z_coords.shape}") 

        # Diagnostics
        n_masked = np.sum(~mask)
        n_fluid = np.sum(mask)
        print(f"  Masked {n_masked} solid cells, {n_fluid} fluid cells remain", file=sys.stderr)

    print(f"DEBUG: z after masking: shape={z_coords.shape}, z range=[{np.min(z_coords):.6f}, {np.max(z_coords):.6f}]", file=sys.stderr)
    # Now remove z not in 0 <= z <=  (lower wall to centerline)
    # z varies only along axis 2, so build a 1D index mask and slice
    z_1d = z_coords[0, 0, :]
    z_idx = np.where((z_1d >= -delta) & (z_1d <= 0))[0]
    z_coords = z_coords[:, :, z_idx]
    u  = u[:, :, z_idx]
    v  = v[:, :, z_idx]
    w  = w[:, :, z_idx]
    uu = uu[:, :, z_idx]
    vv = vv[:, :, z_idx]
    ww = ww[:, :, z_idx]
    print(f"DEBUG: z after limiting z: shape={z_coords.shape}, z range=[{np.nanmin(z_coords):.6f}, {np.nanmax(z_coords):.6f}]", file=sys.stderr)

    # Extract 1D z profile (same for all x,y)
    z_coords_1d = z_coords[0, 0, :]
    print(f"DEBUG: z-coordinates (1D) range: {z_coords_1d}")
    time = float(ds.current_time)
    
    return u, v, w, uu, vv, ww, time, z_coords_1d

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
        U_bar = np.mean(u, axis=periodic_axes)
        V_bar = np.mean(v, axis=periodic_axes)
        W_bar = np.mean(w, axis=periodic_axes)
        # RMS = sqrt(<u'u'>)
        U_rms_prof = np.sqrt(np.mean(uu, axis=periodic_axes))
        V_rms_prof = np.sqrt(np.mean(vv, axis=periodic_axes))
        W_rms_prof = np.sqrt(np.mean(ww, axis=periodic_axes))
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
    ap.add_argument("--utau-source", choices=['gradP', 'gradU'], default='gradP',
                    help="Source for friction velocity: 'gradP' (pressure gradient, default) "
                         "or 'gradU' (velocity gradient at wall)")
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

    # Load the pre-averaged simulation data
    u, v, w, uu, vv, ww, t, zw = load_pltavg_averaged_velocity(pltavg_file, delta=delta, IB=args.IB)
    
    print(f"  loaded {pltavg_file.parent.name}  t={t:.4e}  shape={u.shape}")

    
    stats = wall_normal_stats(u, v, w, uu, vv, ww, periodic_axes=(0, 1), IB=args.IB)
    print(f"  computed wall-normal statistics\n")
    print(f"  U_bar has {np.sum(~np.isnan(stats['U_bar']))} valid values out of {stats['U_bar'].size}", file=sys.stderr)
    

    nu = mu / density

    # Friction velocity
    if args.utau_source == 'gradP':
        u_tau = math.sqrt(dpdx_magnitude * delta / density)
    else:
        # Velocity gradient at wall: tau_w = mu * dU/dy, no-slip => U(0) = 0
        yw = zw + delta
        order = np.argsort(yw)
        j0, j1 = order[0], order[1]
        y0, y1 = yw[j0], yw[j1]
        U0, U1 = stats["U_bar"][j0], stats["U_bar"][j1]
        # Quadratic fit U(y) = a*y + b*y^2 through (0,0), (y0,U0), (y1,U1)
        # dU/dy|_wall = a = (U0*y1^2 - U1*y0^2) / (y0*y1*(y1 - y0))
        dU_dy_wall = (U0 * y1**2 - U1 * y0**2) / (y0 * y1 * (y1 - y0))
        tau_w = mu * abs(dU_dy_wall)
        u_tau = math.sqrt(tau_w / density)
    print(f"Results:")
    print(f"  u_tau = {u_tau:.5g}")
    print(f"  nu    = {nu:.5g}")
    print(f"  delta = {delta:.5g}")
    print(f"  Re_tau = {u_tau * delta / nu:.3f}\n")

    # Wall-units rescale
    # RMS is computed from time-averaged Reynolds stress: RMS = sqrt(<u'u'>)
    # Wall distance from lower wall: y_phys = z + delta (z runs from -delta to 0)
    yplus = (zw + delta) * u_tau / nu
    Uplus = stats["U_bar"] / u_tau
    urms_p = stats["U_rms"] / u_tau  # sqrt(<u'u'>) / u_tau
    vrms_p = stats["V_rms"] / u_tau
    wrms_p = stats["W_rms"] / u_tau

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
