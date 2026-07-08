#!/usr/bin/env python3
"""
Generate input parameters for flat (theta=0) or slanted (theta=45) channel configurations.

For --flat: Single horizontal channel at z_s=0.625
For --slanted: Three diagonal channels arranged in triangular pattern at 45 degrees
"""

import numpy as np
import argparse


def calculate_channel_config(
    flat_mode=False,
    slanted_mode=False,
    x_hi=4.0,
    y_hi=0.0625,
    z_hi=None,
    H=1.0,
    nx=256,
    blocking_factor=4,
    channel_shift=None,
    rho=1.0,
    mu=1.0e-2,
    Umax=1.0,
    align='default',
    nx_align=64,
):
    """
    Calculate channel configuration for flat or slanted modes.
    
    For flat mode: Single channel at theta=0, z_hi=2.0
    For slanted mode: Three channels at theta=45°, z_hi=x_hi (square domain)
    
    Returns:
        dict with configuration parameters and channel segments
    """
    
    # Validate modes
    if not (flat_mode or slanted_mode):
        raise ValueError("Must specify either --flat or --slanted")
    if flat_mode and slanted_mode:
        raise ValueError("Cannot specify both --flat and --slanted")
    
    # Set domain parameters based on mode
    if flat_mode:
        theta_deg = 0.0
        theta_rad = 0.0
        if z_hi is None:
            z_hi = 2.0
    else:  # slanted_mode
        theta_deg = 45.0
        theta_rad = np.deg2rad(45.0)
        if z_hi is None:
            z_hi = x_hi  # Square domain for 45-degree slant
    
    # Physical parameters
    x_lo = 0.0
    y_lo = 0.0
    z_lo = 0.0
    y_mid = y_hi / 2.0
    
    # Mesh resolution: compute cell size from nx
    cell_size = (x_hi - x_lo) / nx
    
    # For convergence studies with theta=0 and cell face alignment,
    # ensure H also aligns with grid cells at all refinement levels
    if np.abs(theta_rad) < 1e-6 and align == 'cf':
        # Calculate reference cell size based on alignment grid
        cell_size_ref = (x_hi - x_lo) / nx_align
        # Round H to be an integer multiple of reference cell size
        n_cells_H = H / cell_size_ref
        n_cells_H_rounded = int(np.round(n_cells_H))
        H_aligned = n_cells_H_rounded * cell_size_ref
        if np.abs(H_aligned - H) > 1e-10:
            print(f"Warning: Aligning H from {H:.6f} to {H_aligned:.6f} for grid convergence")
        H = H_aligned
    
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
    
    # Calculate BodyForce magnitude
    # BodyForce.magnitude = 1/rho * 8 * mu * Umax / H^2 * (cos(theta), 0, sin(theta))
    body_force_coeff = (8.0 * mu * Umax) / (rho * H * H)
    body_force_x = body_force_coeff * np.cos(theta_rad)
    body_force_y = 0.0
    body_force_z = body_force_coeff * np.sin(theta_rad)
    
    # Generate channel segments
    segments = []
    
    if flat_mode:
        # Single flat channel
        if channel_shift is None:
            channel_shift = 0.1
        
        # Calculate centerline z-position for flat channel
        z_center = H / (2.0 * np.cos(theta_rad)) + channel_shift
        
        # Apply grid alignment for cell face alignment
        if align == 'cf':
            cell_size_align = (x_hi - x_lo) / nx_align
            z_center = np.round(z_center / cell_size_align) * cell_size_align
        
        segments.append({
            'label': 's1',
            'start': [x_lo, y_mid, z_center],
            'end': [x_hi, y_mid, z_center],
            'height': H,
        })
    
    else:  # slanted_mode
        # Three slanted channels at 45 degrees
        # Lx = x_hi, Lz = z_hi = x_hi (square domain)
        Lx = x_hi
        Lz = z_hi
        
        # Channel 1: (0, 0, 0) -> (Lx, 0, Lz)
        segments.append({
            'label': 's1',
            'start': [0.0, y_mid, 0.0],
            'end': [Lx, y_mid, Lz],
            'height': H,
        })
        
        # Channel 2: (-Lx/2, 0, Lz/2) -> (Lx/2, 0, 3Lz/2)
        segments.append({
            'label': 's2',
            'start': [-Lx/2.0, y_mid, Lz/2.0],
            'end': [Lx/2.0, y_mid, 3.0*Lz/2.0],
            'height': H,
        })
        
        # Channel 3: (Lx/2, 0, -Lz/2) -> (3Lx/2, 0, Lz/2)
        segments.append({
            'label': 's3',
            'start': [Lx/2.0, y_mid, -Lz/2.0],
            'end': [3.0*Lx/2.0, y_mid, Lz/2.0],
            'height': H,
        })
    
    config = {
        'theta_deg': theta_deg,
        'theta_rad': theta_rad,
        'flat_mode': flat_mode,
        'slanted_mode': slanted_mode,
        'domain_lo': [x_lo, y_lo, z_lo],
        'domain_hi': [x_hi, y_hi, z_hi],
        'n_cell': [n_x, n_y, n_z],
        'nx': nx,
        'blocking_factor': blocking_factor,
        'dx': dx,
        'dy': dy,
        'dz': dz,
        'H': H,
        'channel_shift': channel_shift,
        'rho': rho,
        'mu': mu,
        'Umax': Umax,
        'body_force_x': body_force_x,
        'body_force_y': body_force_y,
        'body_force_z': body_force_z,
        'align': align,
        'nx_align': nx_align,
        'segments': segments,
    }
    
    return config


def print_config(config):
    """Print configuration in a readable format."""
    mode_str = "FLAT" if config['flat_mode'] else "SLANTED (45°)"
    
    print("\n" + "="*70)
    print(f"CHANNEL CONFIGURATION ({mode_str})")
    print("="*70)
    
    print(f"\nAngle:")
    print(f"  θ = {config['theta_deg']:.2f}° = {config['theta_rad']:.6f} rad")
    if np.abs(config['theta_rad']) < 1e-6 and config['align'] == 'cf':
        print(f"  Grid alignment: Cell Face (nx_align={config['nx_align']})")
    
    print(f"\nDomain:")
    lo = config['domain_lo']
    hi = config['domain_hi']
    print(f"  X: [{lo[0]:.4f}, {hi[0]:.4f}]")
    print(f"  Y: [{lo[1]:.4f}, {hi[1]:.6f}]")
    print(f"  Z: [{lo[2]:.4f}, {hi[2]:.4f}]")
    
    print(f"\nChannel Height (perpendicular): H = {config['H']:.6f}")
    
    print(f"\nMesh Resolution (isotropic cell size):")
    print(f"  amr.n_cell = {config['n_cell'][0]} {config['n_cell'][1]} {config['n_cell'][2]}")
    print(f"  Cell sizes:")
    print(f"    Δx = {config['dx']:.6f}")
    print(f"    Δy = {config['dy']:.6f}")
    print(f"    Δz = {config['dz']:.6f}")
    
    print("\n" + "="*70)
    print("Configurable Parameters for .inp file:")
    print("="*70)
    print(f"\ngeometry.prob_lo = {config['domain_lo'][0]:.1f} {config['domain_lo'][1]:.6f} {config['domain_lo'][2]:.1f}")
    print(f"geometry.prob_hi = {config['domain_hi'][0]:.1f} {config['domain_hi'][1]:.6f} {config['domain_hi'][2]:.4f}")
    
    print(f"\namr.n_cell = {config['n_cell'][0]} {config['n_cell'][1]} {config['n_cell'][2]}")
    
    print(f"\nChannelBuilder.segment_labels =", " ".join([seg['label'] for seg in config['segments']]))
    
    for seg in config['segments']:
        label = seg['label']
        start = seg['start']
        end = seg['end']
        h = seg['height']
        
        print(f"\nChannelBuilder.{label}.flow_speed = 1.0")
        print(f"ChannelBuilder.{label}.velocity_profile = Parabolic")
        print(f"ChannelBuilder.{label}.type = Trapezoid")
        print(f"ChannelBuilder.{label}.top_width_start = 1.0")
        print(f"ChannelBuilder.{label}.top_width_end = 1.0")
        print(f"ChannelBuilder.{label}.bottom_width_start = 1.0")
        print(f"ChannelBuilder.{label}.bottom_width_end = 1.0")
        print(f"ChannelBuilder.{label}.height_start = {h:.4f}")
        print(f"ChannelBuilder.{label}.height_end = {h:.4f}")
        print(f"ChannelBuilder.{label}.segment_start_point = {start[0]:.16f} {start[1]:.6f} {start[2]:.16f}")
        print(f"ChannelBuilder.{label}.segment_end_point = {end[0]:.16f} {end[1]:.6f} {end[2]:.16f}")
    
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
        description="Generate flat or slanted channel configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Flat channel (theta=0)
  python SlantedChannelConfig.py --flat --nx 128
  python SlantedChannelConfig.py --flat --align cf --nx_align 64
  
  # Slanted channel (theta=45, three segments)
  python SlantedChannelConfig.py --slanted --nx 32
  python SlantedChannelConfig.py --slanted --nx 64
        """
    )
    
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--flat', action='store_true', help='Flat channel configuration (theta=0)')
    mode_group.add_argument('--slanted', action='store_true', help='Slanted channel configuration (theta=45, 3 segments)')
    
    parser.add_argument('--H', type=float, default=1.0, help='Channel height perpendicular (default: 1.0)')
    parser.add_argument('--nx', type=int, default=None, help='Number of cells in x-direction (default: 64 for flat, 32 for slanted)')
    parser.add_argument('--blocking_factor', type=int, default=4, help='AMR blocking factor (default: 4)')
    parser.add_argument('--channel_shift', type=float, default=0.1, help='Vertical shift of flat channel (default: 0.1)')
    parser.add_argument('--rho', type=float, default=1.0, help='Fluid density (default: 1.0)')
    parser.add_argument('--mu', type=float, default=1.0e-2, help='Dynamic viscosity (default: 1.0e-2)')
    parser.add_argument('--Umax', type=float, default=1.0, help='Maximum velocity (default: 1.0)')
    parser.add_argument('--align', type=str, default='default', choices=['default', 'cf'], help='Grid alignment for flat channel: default (none) or cf (cell face)')
    parser.add_argument('--nx_align', type=int, default=64, help='Reference nx for grid alignment in convergence studies (default: 64)')
    
    args = parser.parse_args()
    
    # Set mode-specific defaults
    if args.flat:
        x_hi = 4.0
        y_hi = 0.0625
        z_hi = 2.0
        nx = args.nx if args.nx is not None else 64
    else:  # slanted
        x_hi = 2.0
        y_hi = 0.0625
        z_hi = 2.0
        nx = args.nx if args.nx is not None else 32
    
    config = calculate_channel_config(
        flat_mode=args.flat,
        slanted_mode=args.slanted,
        x_hi=x_hi,
        y_hi=y_hi,
        z_hi=z_hi,
        H=args.H,
        nx=nx,
        blocking_factor=args.blocking_factor,
        channel_shift=args.channel_shift,
        rho=args.rho,
        mu=args.mu,
        Umax=args.Umax,
        align=args.align,
        nx_align=args.nx_align,
    )
    
    print_config(config)


if __name__ == "__main__":
    main()
