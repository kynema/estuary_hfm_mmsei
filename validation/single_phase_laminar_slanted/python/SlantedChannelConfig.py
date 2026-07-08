"""
Generate input parameters for slanted channel configuration based on wall 
angle theta.

This script calculates:
1. Domain bounds (geometry.prob_lo/hi)
2. Segment parameters (centered vertically in channel)
3. amr.n_cell for isotropic mesh
"""

import numpy as np
import argparse


def calculate_slanted_channel_config(
    theta_deg=None,
    theta_rad=None,
    x_hi=4.0,
    y_hi=0.0625,
    z_hi=2.0,
    H=1.0,
    nx=256,
    blocking_factor=4,
    channel_shift=None,
    rho=1.0,
    mu=1.0e-2,
    Umax=1.0,
    align='default',
):
    """
    Calculate slanted channel configuration.
    
    Args:
        theta_deg: Channel angle in degrees (relative to x-axis)
        theta_rad: Channel angle in radians
        x_hi: Domain length in x-direction
        y_hi: Domain length in y-direction (typically very small for 2D)
        z_hi: Domain length in z-direction
        H: Channel height (perpendicular distance between walls)
        nx: Number of cells in x-direction (default: 256)
        blocking_factor: AMR blocking factor (all n_cell must be divisible by this)
        channel_shift: Vertical shift to move channel up from domain bottom (default: 0.1)
        rho: Fluid density (default: 1.0)
        mu: Dynamic viscosity (default: 1.0e-2)
        Umax: Maximum velocity (default: 1.0)
        align: Grid alignment for channel bottom when theta=0 (default or cf for cell face)
    
    Returns:
        dict with configuration parameters
    """
    
    # Convert angle to radians if needed
    if theta_deg is not None:
        theta = np.deg2rad(theta_deg)
    elif theta_rad is not None:
        theta = theta_rad
    else:
        raise ValueError("Must specify either theta_deg or theta_rad")
    
    # Physical parameters
    x_lo = 0.0
    y_lo = 0.0
    z_lo = 0.0
    y_mid = y_hi / 2.0
    
    # Mesh resolution: compute cell size from nx
    cell_size = (x_hi - x_lo) / nx
    
    # Calculate initial n_cell
    n_x = nx
    n_y = int(np.round(y_hi / cell_size))
    n_z = int(np.round((z_hi - z_lo) / cell_size))
    
    # Ensure n_y is at least blocking_factor
    n_y = max(n_y, blocking_factor)
    
    # Round all n_cell to be divisible by blocking_factor
    def round_to_blocking_factor(n):
        return int(np.ceil(n / blocking_factor) * blocking_factor)
    n_x = round_to_blocking_factor(n_x)
    n_y = round_to_blocking_factor(n_y)
    n_z = round_to_blocking_factor(n_z)
    
    # Recalculate actual cell sizes based on rounded n_cell
    dx = (x_hi - x_lo) / n_x
    
    # Adjust y_hi and z_hi to ensure uniform grid spacing using dx
    y_hi = dx * n_y
    z_hi = z_lo + dx * n_z
    
    dy = y_hi / n_y
    dz = (z_hi - z_lo) / n_z
    
    # Segment points: run from bottom-left to upper-right corner of domain
    # For well-defined channel at corners, segment should be within domain bounds
    # Shift centerline up by H/10 to avoid boundary issues
    # Extend x coordinates beyond domain by H*tan(theta) on each side
    x_start = x_lo - H * np.tan(theta)
    x_end = x_hi + H * np.tan(theta)
    
    # Shift channel up by a configurable amount to avoid boundary issues
    # and ensure terrain cells. Default is 0.1.
    if channel_shift is None:
        channel_shift = 0.1

    # z rises along the slant: tan(theta) per unit x
    # Center the channel at z = H/(2*cos(theta)) + channel_shift at x = x_lo
    z_at_xlo = H / (2.0 * np.cos(theta)) + channel_shift
    
    # At x_start, z should be:
    z_start = z_at_xlo - (x_lo - x_start) * np.tan(theta)
    
    # At x_end, z should be:
    z_end = z_at_xlo + (x_end - x_lo) * np.tan(theta)
    
    # Apply grid alignment for horizontal channel (theta ≈ 0)
    if np.abs(theta) < 1e-6 and align == 'cf':
        # For horizontal channel, z_start and z_end should be equal
        z_center = (z_start + z_end) / 2.0
        # Cell face alignment: snap to nearest cell face (multiple of dz)
        z_aligned = np.round(z_center / dz) * dz
        z_start = z_aligned
        z_end = z_aligned
    
    seg_start = [x_start, y_mid, z_start]
    seg_end = [x_end, y_mid, z_end]

    # Slope of the channel
    seg_slope = (z_end - z_start) / (x_end - x_start)
    seg_intercept = z_start - seg_slope * x_start

    # For trapezoidal segments in ChannelBuilder, the velocity profile uses
    # vertical distance (z-coordinate) to compute the parabolic profile
    height_start = H
    height_end = H
    
    # Calculate BodyForce magnitude
    # BodyForce.magnitude = 1/rho * 8 * mu * Umax / H^2 * (cos(theta), 0, sin(theta))
    body_force_coeff = (8.0 * mu * Umax) / (rho * H * H)
    body_force_x = body_force_coeff * np.cos(theta)
    body_force_y = 0.0
    body_force_z = body_force_coeff * np.sin(theta)
    
    config = {
        'theta_deg': np.rad2deg(theta) if theta_rad is None else theta_deg,
        'theta_rad': theta,
        'domain_lo': [x_lo, y_lo, z_lo],
        'domain_hi': [x_hi, y_hi, z_hi],
        'segment_start': seg_start,
        'segment_end': seg_end,
        'segment_slope': seg_slope,
        'segment_intercept': seg_intercept,
        'n_cell': [n_x, n_y, n_z],
        'nx': nx,
        'blocking_factor': blocking_factor,
        'cell_size': cell_size,
        'dx': dx,
        'dy': dy,
        'dz': dz,
        'H': H,
        'channel_shift': channel_shift,
        'height_start': height_start,
        'height_end': height_end,
        'rho': rho,
        'mu': mu,
        'Umax': Umax,
        'body_force_x': body_force_x,
        'body_force_y': body_force_y,
        'body_force_z': body_force_z,
        'align': align,
    }
    
    return config


def print_config(config):
    """Print configuration in a readable format."""
    print("\n" + "="*70)
    print("SLANTED CHANNEL CONFIGURATION")
    print("="*70)
    
    print(f"\nAngle:")
    print(f"  θ = {config['theta_deg']:.2f}° = {config['theta_rad']:.6f} rad")
    if np.abs(config['theta_rad']) < 1e-6:
        align_str = 'None' if config['align'] == 'default' else 'Cell Face'
        print(f"  Grid alignment (horizontal channel): {align_str}")
    
    print(f"\nDomain:")
    lo = config['domain_lo']
    hi = config['domain_hi']
    print(f"  X: [{lo[0]:.4f}, {hi[0]:.4f}]")
    print(f"  Y: [{lo[1]:.4f}, {hi[1]:.6f}]")
    print(f"  Z: [{lo[2]:.4f}, {hi[2]:.4f}]")
    
    print(f"\nChannel Height (perpendicular): H = {config['H']:.4f}")
    print(f"Vertical extent (z): H/cos(θ) = {config['H']/np.cos(config['theta_rad']):.4f}")
    
    print(f"\nSegment Points:")
    s = config['segment_start']
    e = config['segment_end']
    print(f"  Start: ({s[0]:.4f}, {s[1]:.6f}, {s[2]:.4f})")
    print(f"  End:   ({e[0]:.4f}, {e[1]:.6f}, {e[2]:.4f})")

    print(f"\nSegment Line:")
    print(f"  z = {config['segment_slope']:.6f} * x + {config['segment_intercept']:.6f}")
    print(f"  tan(θ) = {np.tan(config['theta_rad']):.6f}")
    
    print(f"\nMesh Resolution (isotropic cell size):")
    print(f"  amr.n_cell = {config['n_cell'][0]} {config['n_cell'][1]} {config['n_cell'][2]}")
    print(f"  Cell sizes:")
    print(f"    Δx = {config['dx']:.6f}")
    print(f"    Δy = {config['dy']:.6f}")
    print(f"    Δz = {config['dz']:.6f}")
    
    print("\n" +"="*70)
    print(f"Configurable Parameters for .inp file:")
    print("="*70)
    print(f"\ngeometry.prob_lo = {config['domain_lo'][0]:.1f} {config['domain_lo'][1]:.6f} {config['domain_lo'][2]:.1f}")
    print(f"geometry.prob_hi = {config['domain_hi'][0]:.1f} {config['domain_hi'][1]:.6f} {config['domain_hi'][2]:.4f}")
    print(f"\namr.n_cell = {config['n_cell'][0]} {config['n_cell'][1]} {config['n_cell'][2]}")
    print(f"\nChannelBuilder.s1.segment_start_point = {config['segment_start'][0]:.16f} {config['segment_start'][1]:.6f} {config['segment_start'][2]:.16f}")
    print(f"ChannelBuilder.s1.segment_end_point = {config['segment_end'][0]:.16f} {config['segment_end'][1]:.6f} {config['segment_end'][2]:.16f}")
    print(f"\nChannelBuilder.s1.height_start = {config['height_start']:.4f}")
    print(f"ChannelBuilder.s1.height_end = {config['height_end']:.4f}")
    
    print(f"\nFluid Properties:")
    print(f"  ρ (density) = {config['rho']:.4f}")
    print(f"  μ (viscosity) = {config['mu']:.6f}")
    print(f"  Umax = {config['Umax']:.4f}")
    
    print(f"\nBody Force (from Poiseuille equation):")
    print(f"  Formula: 1/ρ * 8*μ*Umax/H² * (cos(θ), 0, sin(θ))")
    print(f"  BodyForce.magnitude = {config['body_force_x']:.16f} {config['body_force_y']:.1f} {config['body_force_z']:.16f}")
    
    print("\n" + "="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate slanted channel configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python SlantedChannelConfig.py --theta_deg 7.5
  python SlantedChannelConfig.py --theta_deg 15 --nx 256
  python SlantedChannelConfig.py --theta_rad 0.2618 --nx 512
  python SlantedChannelConfig.py --theta_deg 30 --x_hi 5.0 --H 2.0 --nx 1024
        """
    )
    
    parser.add_argument('--theta_deg', type=float, help='Wall angle in degrees')
    parser.add_argument('--theta_rad', type=float, help='Wall angle in radians')
    parser.add_argument('--x_hi', type=float, default=4.0, help='Domain length in x (default: 4.0)')
    parser.add_argument('--y_hi', type=float, default=0.0625, help='Domain length in y (default: 0.0625)')
    parser.add_argument('--z_hi', type=float, default=2.0, help='Domain length in z (default: 2.0)')
    parser.add_argument('--H', type=float, default=1.0, help='Channel height perpendicular (default: 1.0)')
    parser.add_argument('--nx', type=int, default=512, help='Number of cells in x-direction (default: 512)')
    parser.add_argument('--blocking_factor', type=int, default=4, help='AMR blocking factor (default: 4)')
    parser.add_argument('--channel_shift', type=float, default=0.1, help='Vertical shift of channel up from domain bottom (default: 0.1)')
    parser.add_argument('--rho', type=float, default=1.0, help='Fluid density (default: 1.0)')
    parser.add_argument('--mu', type=float, default=1.0e-2, help='Dynamic viscosity (default: 1.0e-2)')
    parser.add_argument('--Umax', type=float, default=1.0, help='Maximum velocity (default: 1.0)')
    parser.add_argument('--align', type=str, default='default', choices=['default', 'cf'], help='Grid alignment when theta=0: default (none) or cf (cell face)')
    
    args = parser.parse_args()
    
    if args.theta_deg is None and args.theta_rad is None:
        parser.print_help()
        return
    
    config = calculate_slanted_channel_config(
        theta_deg=args.theta_deg,
        theta_rad=args.theta_rad,
        x_hi=args.x_hi,
        y_hi=args.y_hi,
        H=args.H,
        nx=args.nx,
        blocking_factor=args.blocking_factor,
        channel_shift=args.channel_shift,
        rho=args.rho,
        mu=args.mu,
        Umax=args.Umax,
        align=args.align,
    )
    
    print_config(config)


if __name__ == "__main__":
    main()
