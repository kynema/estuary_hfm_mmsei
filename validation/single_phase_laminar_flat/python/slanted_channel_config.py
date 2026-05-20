"""
Generate slanted channel configuration based on wall angle theta.

This script calculates:
1. Domain bounds (geometry.prob_lo/hi)
2. Segment start/end points (centered vertically in channel)
3. amr.n_cell for isotropic mesh
"""

import numpy as np
import argparse


def calculate_slanted_channel_config(
    theta_deg=None,
    theta_rad=None,
    x_hi=4.0,
    y_hi=0.0625,
    H=1.0,
    nx=512,
    blocking_factor=4,
):
    """
    Calculate slanted channel configuration.
    
    Args:
        theta_deg: Channel angle in degrees (relative to x-axis)
        theta_rad: Channel angle in radians
        x_hi: Domain length in x-direction
        y_hi: Domain extent in y-direction (typically very small for 2D)
        H: Channel height (perpendicular distance between walls)
        nx: Number of cells in x-direction (default: 512)
        blocking_factor: AMR blocking factor (all n_cell must be divisible by this)
    
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
    
    # Shift channel up by H/10 to avoid boundary issues and ensure terrain cells
    channel_shift = H / 10.0
    
    # Domain extent in z (vertical)
    # z_hi = (x_hi - x_lo)*tan(theta) + H/cos(theta) + shift
    # Extended to give extra margin
    z_hi = 1.2 * ((x_hi - x_lo) * np.tan(theta) + H / np.cos(theta) + channel_shift)
    
    # Mesh resolution: compute cell size from nx
    cell_size = (x_hi - x_lo) / nx
    
    # Calculate initial n_cell
    n_x = nx
    n_y = int(np.round(y_hi / cell_size))
    n_z = int(np.round((z_hi - z_lo) / cell_size))
    
    # Ensure n_y is at least blocking_factor (typical for thin domains)
    n_y = max(n_y, blocking_factor)
    
    # Round all n_cell to be divisible by blocking_factor
    n_x = int(np.ceil(n_x / blocking_factor) * blocking_factor)
    n_y = int(np.ceil(n_y / blocking_factor) * blocking_factor)
    n_z = int(np.ceil(n_z / blocking_factor) * blocking_factor)
    
    # Recalculate actual cell sizes based on rounded n_cell
    actual_dx = (x_hi - x_lo) / n_x
    actual_dy = y_hi / n_y
    actual_dz = (z_hi - z_lo) / n_z
    
    # Segment points: run from bottom-left to upper-right corner of domain
    # For well-defined channel at corners, segment should be within domain bounds
    # Shift centerline up by H/10 to avoid boundary issues
    # Extend x coordinates beyond domain by H*tan(theta) on each side
    x_start = x_lo - H * np.tan(theta)
    x_end = x_hi + H * np.tan(theta)
    
    # z rises along the slant: tan(theta) per unit x
    # Center the channel at z = H/(2*cos(theta)) + channel_shift at x = x_lo
    z_at_xlo = H / (2.0 * np.cos(theta)) + channel_shift
    
    # At x_start, z should be:
    z_start = z_at_xlo - (x_lo - x_start) * np.tan(theta)
    
    # At x_end, z should be:
    z_end = z_at_xlo + (x_end - x_lo) * np.tan(theta)
    
    seg_start = [x_start, y_mid, z_start]
    seg_end = [x_end, y_mid, z_end]
    
    # For sloped segment, the perpendicular cross-section height must be
    # H / cos(theta) to achieve z-extent of H in world coordinates.
    # Keep this constant along the entire segment.
    height_start = H / np.cos(theta)
    height_end = H / np.cos(theta)
    
    # Initial velocity aligned with segment direction
    # For a segment sloped at angle theta in the xz-plane,
    # a unit vector along the segment is (cos(theta), 0, sin(theta))
    # With flow_speed = 1.0, velocity components are:
    velocity_x = 1.0 * np.cos(theta)
    velocity_y = 0.0
    velocity_z = 1.0 * np.sin(theta)
    
    config = {
        'theta_deg': np.rad2deg(theta) if theta_rad is None else theta_deg,
        'theta_rad': theta,
        'domain_lo': [x_lo, y_lo, z_lo],
        'domain_hi': [x_hi, y_hi, z_hi],
        'segment_start': seg_start,
        'segment_end': seg_end,
        'n_cell': [n_x, n_y, n_z],
        'nx': nx,
        'blocking_factor': blocking_factor,
        'cell_size': cell_size,
        'actual_dx': actual_dx,
        'actual_dy': actual_dy,
        'actual_dz': actual_dz,
        'H': H,
        'channel_shift': channel_shift,
        'height_start': height_start,
        'height_end': height_end,
        'velocity_x': velocity_x,
        'velocity_y': velocity_y,
        'velocity_z': velocity_z,
    }
    
    return config


def print_config(config):
    """Print configuration in a readable format."""
    print("\n" + "="*70)
    print("SLANTED CHANNEL CONFIGURATION")
    print("="*70)
    
    print(f"\nAngle:")
    print(f"  θ = {config['theta_deg']:.2f}° = {config['theta_rad']:.6f} rad")
    
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
    
    
    print(f"\nMesh Resolution (isotropic cell size):")
    print(f"  nx = {config['nx']}")
    print(f"  amr.n_cell = {config['n_cell'][0]} {config['n_cell'][1]} {config['n_cell'][2]}")
    print(f"  Cell sizes:")
    print(f"    Δx = {config['actual_dx']:.6f}")
    print(f"    Δy = {config['actual_dy']:.6f}")
    print(f"    Δz = {config['actual_dz']:.6f}")
    
    print(f"\nConfigurable Parameters for .inp file:")
    print(f"\n  geometry.prob_lo = {config['domain_lo'][0]:.1f} {config['domain_lo'][1]:.6f} {config['domain_lo'][2]:.1f}")
    print(f"  geometry.prob_hi = {config['domain_hi'][0]:.1f} {config['domain_hi'][1]:.6f} {config['domain_hi'][2]:.4f}")
    print(f"\n  incflo.gravity = 0.0 0.0 0.0")
    print(f"  incflo.velocity = {config['velocity_x']:.4f} {config['velocity_y']:.1f} {config['velocity_z']:.4f}  # Aligned with segment")
    print(f"\n  amr.n_cell = {config['n_cell'][0]} {config['n_cell'][1]} {config['n_cell'][2]}")
    print(f"\n  ChannelBuilder.s1.segment_start_point = {config['segment_start'][0]:.16f} {config['segment_start'][1]:.16f} {config['segment_start'][2]:.16f}")
    print(f"  ChannelBuilder.s1.segment_end_point = {config['segment_end'][0]:.16f} {config['segment_end'][1]:.16f} {config['segment_end'][2]:.16f}")
    print(f"\n  ChannelBuilder.s1.height_start = {config['height_start']:.4f}")
    print(f"  ChannelBuilder.s1.height_end = {config['height_end']:.4f}")
    
    print("\n" + "="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate slanted channel configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python slanted_channel_config.py --theta_deg 7.5
  python slanted_channel_config.py --theta_deg 15 --nx 256
  python slanted_channel_config.py --theta_rad 0.2618 --nx 512
  python slanted_channel_config.py --theta_deg 30 --x_hi 5.0 --H 2.0 --nx 1024
        """
    )
    
    parser.add_argument('--theta_deg', type=float, help='Wall angle in degrees')
    parser.add_argument('--theta_rad', type=float, help='Wall angle in radians')
    parser.add_argument('--x_hi', type=float, default=4.0, help='Domain length in x (default: 4.0)')
    parser.add_argument('--y_hi', type=float, default=0.0625, help='Domain extent in y (default: 0.0625)')
    parser.add_argument('--H', type=float, default=1.0, help='Channel height perpendicular (default: 1.0)')
    parser.add_argument('--nx', type=int, default=512, help='Number of cells in x-direction (default: 512)')
    parser.add_argument('--blocking_factor', type=int, default=4, help='AMR blocking factor (default: 4)')
    
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
    )
    
    print_config(config)


if __name__ == "__main__":
    main()
