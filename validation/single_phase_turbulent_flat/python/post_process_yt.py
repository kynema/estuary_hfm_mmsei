#!/usr/bin/env python3
"""
Post-processor for kynema-sgf channel-flow pltAvg files.

Reads pltAvg output, computes spatial and temporal statistics,
rescales to wall units, and overlays simulation results on experimental data.

Workflow
========

1. Automatically loads pltAvg files from cases/ReTau{Re}/ (or ReTau{Re}_IB/).
2. Extract velocity components (u, v, w) and wall-normal coordinate (z).
3. For IB cases, identify first fluid cell using terrain_blank field.
4. Spatially average in periodic directions (x, y by default) to get
   U_mean(z), V_mean(z), W_mean(z) plus variances.
5. Time-average across all input pltAvg files.
6. Compute friction velocity from pressure gradient via case_setup.
7. Wall-units rescale and overlay on experimental data via plot_data.

Examples
========

Re=180 case without IB::

    python post_process.py --Re 180

Re=180 case with IB::

    python post_process.py --Re 180 --IB

Output directories:
    - figures/ReTau180/          (for Re=180 without IB)
    - figures/ReTau180_IB/       (for Re=180 with IB)

Outputs (profiles.csv, Uplus.png, VelRMSplus.png):
  - profiles.csv       : tabular z_phys, y+, U+, urms+, vrms+, wrms+
  - Uplus.png          : U+ vs y+ overlaid on experimental DNS data
  - VelRMSplus.png     : urms+, vrms+, wrms+ vs y+ overlaid on experimental data
"""

from __future__ import annotations

import argparse
import glob
import math
import sys
from pathlib import Path

import numpy as np

import case_setup
from plot_data import plot_mean_velocity_profile, plot_rms_velocity_profiles

# --------------------------------------------------------------------
# Plotfile loading
# --------------------------------------------------------------------

def find_pltavg_files(Re: int, IB: bool = False) -> list[Path]:
    """
    Find all pltAvg files for the given Reynolds number and IB variant.
    
    Parameters
    ----------
    Re : int
        Stress Reynolds number
    IB : bool
        If True, look for ReTau{Re}_IB/ directory
    
    Returns
    -------
    list[Path]
        Sorted list of pltAvg Header file paths
    """
    # Construct case directory path
    case_dir = Path(__file__).parent.parent / "cases"
    
    if IB:
        pltavg_pattern = case_dir / f"ReTau{Re}_IB" / "pltAvg*" / "Header"
    else:
        pltavg_pattern = case_dir / f"ReTau{Re}" / "pltAvg*" / "Header"
    
    files = sorted(glob.glob(str(pltavg_pattern)))
    
    if not files:
        case_variant = f"ReTau{Re}_IB" if IB else f"ReTau{Re}"
        raise FileNotFoundError(
            f"No pltAvg files found in {case_dir / case_variant}/"
        )
    
    return [Path(f) for f in files]


def get_output_dir(Re: int, IB: bool = False) -> Path:
    """
    Get the output directory for figures.
    
    Parameters
    ----------
    Re : int
        Stress Reynolds number
    IB : bool
        If True, use ReTau{Re}_IB/ directory
    
    Returns
    -------
    Path
        Output directory path
    """
    figures_dir = Path(__file__).parent.parent / "figures"
    
    if IB:
        return figures_dir / f"ReTau{Re}_IB"
    else:
        return figures_dir / f"ReTau{Re}"


def load_pltavg_velocity(plt_path: Path, IB: bool = False):
    """
    Load (u, v, w) from a kynema-sgf pltAvg file via yt.
    
    For IB cases, extracts velocities from first fluid cell (not cut by IB).
    
    Parameters
    ----------
    plt_path : Path
        Path to pltAvg directory or its Header file
    IB : bool
        If True, use terrain_blank to identify first fluid cell
    
    Returns
    -------
    u, v, w : np.ndarray
        Velocity components (shape Nx, Ny, Nz)
    prob_lo, prob_hi : tuple of length 3
        Domain bounds
    time : float
        Simulation time
    """
    import yt 

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

    # Load velocity components
    u = np.asarray(cg["boxlib", "velocityx"])  # flow direction
    v = np.asarray(cg["boxlib", "velocityy"])  # width (maps to z in data coords)
    w = np.asarray(cg["boxlib", "velocityz"])  # wall-normal (maps to y in data coords)
    
    # For IB cases, mask out solid cells using terrain_blank
    if IB:
        try:
            terrain_blank = np.asarray(cg["boxlib", "terrain_blank"])
            # terrain_blank = 1 for fluid, 0 for solid
            # Set solid cell velocities to NaN for proper averaging
            u[terrain_blank == 0] = np.nan
            v[terrain_blank == 0] = np.nan
            w[terrain_blank == 0] = np.nan
        except Exception as e:
            print(f"Warning: Could not load terrain_blank field: {e}", file=sys.stderr)
    
    prob_lo = tuple(float(x) for x in ds.domain_left_edge)
    prob_hi = tuple(float(x) for x in ds.domain_right_edge)
    time = float(ds.current_time)
    
    return u, v, w, prob_lo, prob_hi, time

# --------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------

def wall_normal_stats(u, v, w, periodic_axes=(0, 1), IB: bool = False):
    """
    Spatially average over periodic axes (x, y by default) to get
    per-z profiles of mean and variance.
    
    For IB cases, uses nanmean/nanvar to ignore NaN cells.
    
    Returns
    -------
    dict
        Dicts of arrays indexed by wall-normal index k:
        Keys: U_bar, V_bar, W_bar, U_var, V_var, W_var
        Values: 1D arrays (length Nz for non-IB, slightly less for IB due to NaN masking)
    """
    if IB:
        # Use nanmean/nanvar to skip solid cells
        U_bar = np.nanmean(u, axis=periodic_axes)
        V_bar = np.nanmean(v, axis=periodic_axes)
        W_bar = np.nanmean(w, axis=periodic_axes)
        U_var = np.nanvar(u, axis=periodic_axes)
        V_var = np.nanvar(v, axis=periodic_axes)
        W_var = np.nanvar(w, axis=periodic_axes)
    else:
        U_bar = u.mean(axis=periodic_axes)
        V_bar = v.mean(axis=periodic_axes)
        W_bar = w.mean(axis=periodic_axes)
        U_var = u.var(axis=periodic_axes)
        V_var = v.var(axis=periodic_axes)
        W_var = w.var(axis=periodic_axes)
    
    return {
        "U_bar": U_bar, "V_bar": V_bar, "W_bar": W_bar,
        "U_var": U_var, "V_var": V_var, "W_var": W_var
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
    ap.add_argument("--IB", action="store_true",
                    help="Use immersed boundary variant (reads from ReTau{Re}_IB/)")
    ap.add_argument("--utau-source", choices=("gradP", "wall_stress"),
                    default="gradP",
                    help="how to compute u_tau; gradP uses pressure gradient "
                         "from case_setup (default); wall_stress uses dU/dz at wall")
    args = ap.parse_args()

    # Find pltAvg files and set output dir
    try:
        pltavg_files = find_pltavg_files(args.Re, IB=args.IB)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    outdir = get_output_dir(args.Re, IB=args.IB)

    print(f"Post-processing Re_tau = {args.Re}, IB = {args.IB}")
    print(f"Found {len(pltavg_files)} pltAvg files")
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

    # Load and accumulate statistics
    accum = None
    accum_count = 0
    geom = None
    
    for pf in pltavg_files:
        u, v, w, prob_lo, prob_hi, t = load_pltavg_velocity(pf, IB=args.IB)
        
        if geom is None:
            geom = (u.shape, prob_lo, prob_hi)
        else:
            if u.shape != geom[0]:
                print(f"shape mismatch for {pf}: {u.shape} vs {geom[0]}",
                      file=sys.stderr)
                return 2
        
        stats = wall_normal_stats(u, v, w, periodic_axes=(0, 1), IB=args.IB)
        
        if accum is None:
            accum = {k: v.copy() for k, v in stats.items()}
        else:
            for k in accum:
                accum[k] += stats[k]
        
        accum_count += 1
        print(f"  loaded {pf.parent.name}  t={t:.4e}  shape={u.shape}")

    if accum is None or accum_count == 0:
        print("no pltAvg files loaded", file=sys.stderr)
        return 2

    # Time-average
    for k in accum:
        accum[k] /= accum_count

    print(f"  time-averaged {accum_count} files\n")

    Nz = accum["U_bar"].size
    shape, prob_lo, prob_hi = geom

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
    U_bar_filt = accum["U_bar"][mask]
    V_bar_filt = accum["V_bar"][mask]
    W_bar_filt = accum["W_bar"][mask]
    U_var_filt = accum["U_var"][mask]
    V_var_filt = accum["V_var"][mask]
    W_var_filt = accum["W_var"][mask]

    # Friction velocity
    if args.utau_source == "gradP":
        u_tau = math.sqrt(dpdx_magnitude * delta / density)
    else:
        # Use velocity gradient at wall: tau_w = mu * dU/dz
        order = np.argsort(zw)
        k0 = order[0]
        dU_dz_wall = U_bar_filt[k0] / zw[k0]
        tau_w = mu * abs(dU_dz_wall)
        u_tau = math.sqrt(tau_w / density)
    
    nu = mu / density
    print(f"Results:")
    print(f"  u_tau = {u_tau:.5g}  (source: {args.utau_source})")
    print(f"  nu    = {nu:.5g}")
    print(f"  delta = {delta:.5g}")
    print(f"  Re_tau = {u_tau * delta / nu:.3f}\n")

    # Wall-units rescale
    yplus = zw * u_tau / nu  # y+ in data coordinate frame
    Uplus = U_bar_filt / u_tau
    urms_p = np.sqrt(np.maximum(U_var_filt, 0.0)) / u_tau
    vrms_p = np.sqrt(np.maximum(V_var_filt, 0.0)) / u_tau
    wrms_p = np.sqrt(np.maximum(W_var_filt, 0.0)) / u_tau

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
    csv_path = outdir / "profiles.csv"
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
        outpath=outdir / "Uplus.png"
    )
    print(f"  wrote {outdir / 'Uplus.png'}")

    # Plot RMS velocity profiles with simulation data overlay
    plot_rms_velocity_profiles(
        args.Re,
        yplus_sim=yplus_o,
        urms_sim=urms_o,
        vrms_sim=vrms_o,  # V_var = spanwise (velocityy in kynema-sgf)
        wrms_sim=wrms_o,  # W_var = wall-normal (velocityz in kynema-sgf)
        label_sim="kynema-sgf",
        outpath=outdir / "VelRMSplus.png"
    )
    print(f"  wrote {outdir / 'VelRMSplus.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())