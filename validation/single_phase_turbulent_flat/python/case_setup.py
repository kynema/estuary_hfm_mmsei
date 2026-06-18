#!/usr/bin/env python3
"""
Setup domain parameters for turbulent channel flow validation case.
Based on DNS configuration at Re_tau = 180

This script calculates:
1. Domain bounds (geometry.prob_lo/hi)
2. Grid resolution (amr.n_cell) - divisible by blocking factor
3. Cell spacing maintained isotropic in all directions
"""

import argparse


def get_pressure_gradient(Re):
    """
    Get pressure gradient (dpdx) for a given Reynolds number.
    """

    Re_to_dpdx = {
        180: -709.79,
        395: -426.121,
        934: -2382.50,
    }
    
    if Re not in Re_to_dpdx:
        raise ValueError(
            f"Reynolds number {Re} not in database. "
            f"Available: {list(Re_to_dpdx.keys())}"
        )
    
    return Re_to_dpdx[Re]


def get_refinement_boxes(delta, Lx, Ly, density, mu, u_tau, IB=False):
    """
    Define AMR refinement boxes for boundary layer capture.
    Refinement extents scale with Reynolds number via viscous length scale.
    
    Two levels of refinement:
    - Level 1: Outer wall boxes capturing buffer layer (y+ ~ 50)
    - Level 2: Inner wall boxes capturing viscous sublayer (y+ ~ 5)
    
    The refinement extent is based on the viscous length scale:
    l_nu = nu / u_tau
    
    Parameters
    ----------
    delta : float
        Channel half width in meters
    Lx, Ly : float
        Domain extent in x and y
    density, mu, u_tau : float
        Fluid properties (density, dynamic viscosity, friction velocity)
    IB : bool
        Use immersed boundary
    
    Returns
    -------
    dict
        Refinement box configuration for input file
    """
    
    # Calculate viscous length scale from provided parameters
    nu = mu / density  # Kinematic viscosity
    l_nu = nu / u_tau  # Viscous length scale
    
    # Refinement box extents in wall units (y+ value)
    # Level 1: Capture buffer layer up to y+ ~ 40-60
    # Level 2: Capture viscous sublayer up to y+ ~ 5
    yplus_level1 = 50.0
    yplus_level2 = 5.0
    
    # Convert from wall units to physical distance
    z_extent_level1 = yplus_level1 * l_nu
    z_extent_level2 = yplus_level2 * l_nu
    
    # Ensure not larger than delta
    z_extent_level1 = min(z_extent_level1, 0.5 * delta)
    z_extent_level2 = min(z_extent_level2, 0.25 * delta)
    
    if IB:
        z_max = 0.0052
    else:
        z_max = delta
    
    boxes = {
        'level1_bottom': {
            'origin': [0.0, 0.0, -z_max],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_level1],
        },
        'level1_top': {
            'origin': [0.0, 0.0, z_max - z_extent_level1],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_level1],
        },
        'level2_bottom': {
            'origin': [0.0, 0.0, -z_max],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_level2],
        },
        'level2_top': {
            'origin': [0.0, 0.0, z_max - z_extent_level2],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_level2],
        },
        'viscous_length_scale': l_nu,
        'z_extent_level1_yplus': yplus_level1,
        'z_extent_level2_yplus': yplus_level2,
    }
    
    return boxes


def domain_and_flow(delta=0.005, Nx=384, IB=False, blocking_factor=8, Re=180):
    """
    Compute domain boundaries and grid parameters for turbulent channel case.
    
    Parameters
    ----------
    delta : float
        Channel half width in meters (default 0.005 m)
    Nx : int
        Number of cells in x-direction
    IB : bool
        Use immersed boundary (extends domain in z)
    blocking_factor : int
        AMR blocking factor (default 8; must divide Nx, Ny, Nz)
    Re : int
        Stress Reynolds number Re_tau (default 180)
    
    Returns
    -------
    dict
        Dictionary with domain parameters
    """
    
    # Domain dimensions from Re_tau = 180 DNS case
    # x: 6.24δ, y: 3.12δ, z: 2.0δ (no IB) or 2.08δ (with IB)
    Lx = 6.24 * delta
    Ly = 3.12 * delta
    Lz_no_ib = 2.0 * delta
    Lz_ib = 2.08 * delta
    
    Lz = Lz_ib if IB else Lz_no_ib
    
    # Compute grid spacing and cell counts
    dx = Lx / Nx
    
    # Compute Ny and Nz based on maintaining uniform spacing
    Ny = int(round(Ly / dx))
    Nz = int(round(Lz / dx))
    
    # Adjust to be divisible by blocking factor
    Nx = (Nx // blocking_factor) * blocking_factor
    Ny = (Ny // blocking_factor) * blocking_factor
    Nz = (Nz // blocking_factor) * blocking_factor
    
    # Domain boundaries
    # Center domain in y (periodic), asymmetric in z (walls)
    if IB:
        z_extent = 0.0052  # Slightly larger than delta for IB
        prob_lo = [0.0, 0.0, -z_extent]
        prob_hi = [Lx, Ly, z_extent]
    else:
        prob_lo = [0.0, 0.0, -delta]
        prob_hi = [Lx, Ly, delta]
    
    # Verify blocking factor divisibility
    assert Nx % blocking_factor == 0, f"Nx={Nx} not divisible by {blocking_factor}"
    assert Ny % blocking_factor == 0, f"Ny={Ny} not divisible by {blocking_factor}"
    assert Nz % blocking_factor == 0, f"Nz={Nz} not divisible by {blocking_factor}"
    
    # Get flow characteristics for specified stress Reynolds number
    density = 0.468793
    mu = 3.57816e-5  # Dynamic viscosity (air at 750K)
    u_tau = Re * mu / (density * delta)  # Friction velocity
    tau_w = density * u_tau**2  # Wall shear stress
    t_star = delta / u_tau  # Time scale based on friction velocity and channel half-width
    dpdx = get_pressure_gradient(Re)

    # Body force acceleration: F = |dp/dx| / rho (units: m/s^2)
    body_force = -dpdx / density
    
    # Get refinement boxes (using computed physical parameters)
    refinement_boxes = get_refinement_boxes(delta, Lx, Ly, density, mu, u_tau, IB=IB)
    
    return {
        'delta': delta,
        'Nx': Nx,
        'Ny': Ny,
        'Nz': Nz,
        'blocking_factor': blocking_factor,
        'IB': IB,
        'Lx': Lx,
        'Ly': Ly,
        'Lz': Lz,
        'dx': dx,
        'prob_lo': prob_lo,
        'prob_hi': prob_hi,
        'Re': Re,
        'density': density,
        'mu': mu,
        'u_tau': u_tau,
        'tau_w': tau_w,
        'dpdx': dpdx,
        't_star': t_star,
        'body_force': body_force,
        'refinement_boxes': refinement_boxes,
    }


def print_config(config):
    """Print configuration in a readable format."""
    print("\n" + "="*60)
    print("TURBULENT CHANNEL CONFIGURATION")
    print("="*60)
    
    print(f"\nFlow characteristics:")
    print(f"  Stress Reynolds number (Re_τ): {config['Re']}")
    print(f"  Friction velocity (u_τ):       {config['u_tau']:.4f} m/s")
    print(f"  Wall shear stress (τ_w):       {config['tau_w']:.4f} Pa")
    print(f"  Pressure gradient (dp/dx):     {config['dpdx']:.2f} Pa/m")
    print(f"  Time scale (t* = δ/u_τ):       {config['t_star']:.4e} s")
    print(f"  Density:                  {config['density']} kg/m^3")
    print(f"  Dynamic viscosity:        {config['mu']:.5e} Pa.s")
    

    
    
    print(f"\nDomain Dimensions:")
    lo = config['prob_lo']
    hi = config['prob_hi']
    print(f"  Immersed Boundary (IB):     {config['IB']}")
    print(f"  Channel half-width (δ):     {config['delta']:.6f} m")
    print(f"  X: [{lo[0]:.6f}, {hi[0]:.6f}] m, Lx = {config['Lx']/config['delta']:.2f}δ")
    print(f"  Y: [{lo[1]:.6f}, {hi[1]:.6f}] m, Ly = {config['Ly']/config['delta']:.2f}δ")
    print(f"  Z: [{lo[2]:.6f}, {hi[2]:.6f}] m, Lz = {config['Lz']/config['delta']:.2f}δ")
    
    print(f"\nMesh Resolution:")
    print(f"  Blocking Factor: {config['blocking_factor']}")
    print(f"  L0 = {config['Nx']} {config['Ny']} {config['Nz']}")
    print(f"  L1 = {config['Nx']*2} {config['Ny']*2} {config['Nz']*2}")
    print(f"  L2 = {config['Nx']*4} {config['Ny']*4} {config['Nz']*4}")
    print(f"  Cell spacing (L0): dx = dy = dz = {config['dx']:.8f} m")
    
    print("\n" + "="*60)
    print("Parameters for base-turbulent-flat.inp file:")
    print("="*60)
    print(f"\ngeometry.prob_lo = {lo[0]:.1f} {lo[1]:.6f} {lo[2]:.6f}")
    print(f"geometry.prob_hi = {hi[0]:.6f} {hi[1]:.6f} {hi[2]:.6f}")
    print(f"\namr.n_cell = {config['Nx']} {config['Ny']} {config['Nz']}")
    
    print(f"\nICNS.source_terms = BodyForce")

    print("\n" + "="*60)
    fname = f"turbulent-flat-re-{config['Re']}"
    if config['IB']:
        fname += "-IB"
    print(f"Parameters for {fname}.inp file:")
    print("="*60)
    print(f"time.stop_time = {20*config['t_star']:.6f}  # ~20 flow-through time")
    print(f"#time.stop_time = {30*config['t_star']:.6f}  # ~30 flow-through time for statistics")
    print(f"\nBodyForce.magnitude = {config['body_force']:.16f} 0.0 0.0 # Force acceleration in m/s^2")
    
    boxes = config['refinement_boxes']
    l_nu = boxes['viscous_length_scale']
    yp1 = boxes['z_extent_level1_yplus']
    yp2 = boxes['z_extent_level2_yplus']
    
    z_ext_l1 = boxes['level1_bottom']['zaxis'][2]
    z_ext_l2 = boxes['level2_bottom']['zaxis'][2]
    
    print(f"\n# Viscous length scale (l_ν = ν/u_τ): {l_nu:.2e} m")
    print(f"# Level 1 extent: y⁺ ≈ {yp1:.1f} → z_extent = {z_ext_l1:.6f} m ({z_ext_l1/config['delta']:.3f}δ)")
    print(f"# Level 2 extent: y⁺ ≈ {yp2:.1f} → z_extent = {z_ext_l2:.6f} m ({z_ext_l2/config['delta']:.3f}δ)")
    
    print(f"\ntagging.labels = l1_bottom l1_top l2_bottom l2_top")
    
    # Level 1 Bottom Refinement
    print(f"\n# Level 1 Refinement Bottom (buffer layer, y⁺ ~ {yp1:.0f})")
    print(f"tagging.l1_bottom.type = GeometryRefinement")
    print(f"tagging.l1_bottom.shapes = l1_bottom")
    print(f"tagging.l1_bottom.max_level = 1")
    print(f"\ntagging.l1_bottom.l1_bottom.type = box")
    print(f"tagging.l1_bottom.l1_bottom.origin = {boxes['level1_bottom']['origin'][0]:.6f} {boxes['level1_bottom']['origin'][1]:.6f} {boxes['level1_bottom']['origin'][2]:.6f}")
    print(f"tagging.l1_bottom.l1_bottom.xaxis = {boxes['level1_bottom']['xaxis'][0]:.6f} {boxes['level1_bottom']['xaxis'][1]:.6f} {boxes['level1_bottom']['xaxis'][2]:.6f}")
    print(f"tagging.l1_bottom.l1_bottom.yaxis = {boxes['level1_bottom']['yaxis'][0]:.6f} {boxes['level1_bottom']['yaxis'][1]:.6f} {boxes['level1_bottom']['yaxis'][2]:.6f}")
    print(f"tagging.l1_bottom.l1_bottom.zaxis = {boxes['level1_bottom']['zaxis'][0]:.6f} {boxes['level1_bottom']['zaxis'][1]:.6f} {boxes['level1_bottom']['zaxis'][2]:.6f}")
    
    # Level 1 Top Refinement
    print(f"\n# Level 1 Refinement Top (buffer layer, y⁺ ~ {yp1:.0f})")
    print(f"tagging.l1_top.type = GeometryRefinement")
    print(f"tagging.l1_top.shapes = l1_top")
    print(f"tagging.l1_top.max_level = 1")
    print(f"\ntagging.l1_top.l1_top.type = box")
    print(f"tagging.l1_top.l1_top.origin = {boxes['level1_top']['origin'][0]:.6f} {boxes['level1_top']['origin'][1]:.6f} {boxes['level1_top']['origin'][2]:.6f}")
    print(f"tagging.l1_top.l1_top.xaxis = {boxes['level1_top']['xaxis'][0]:.6f} {boxes['level1_top']['xaxis'][1]:.6f} {boxes['level1_top']['xaxis'][2]:.6f}")
    print(f"tagging.l1_top.l1_top.yaxis = {boxes['level1_top']['yaxis'][0]:.6f} {boxes['level1_top']['yaxis'][1]:.6f} {boxes['level1_top']['yaxis'][2]:.6f}")
    print(f"tagging.l1_top.l1_top.zaxis = {boxes['level1_top']['zaxis'][0]:.6f} {boxes['level1_top']['zaxis'][1]:.6f} {boxes['level1_top']['zaxis'][2]:.6f}")
    
    # Level 2 Bottom Refinement
    print(f"\n# Level 2 Refinement Bottom (viscous sublayer, y⁺ ~ {yp2:.0f})")
    print(f"tagging.l2_bottom.type = GeometryRefinement")
    print(f"tagging.l2_bottom.shapes = l2_bottom")
    print(f"tagging.l2_bottom.max_level = 2")
    print(f"\ntagging.l2_bottom.l2_bottom.type = box")
    print(f"tagging.l2_bottom.l2_bottom.origin = {boxes['level2_bottom']['origin'][0]:.6f} {boxes['level2_bottom']['origin'][1]:.6f} {boxes['level2_bottom']['origin'][2]:.6f}")
    print(f"tagging.l2_bottom.l2_bottom.xaxis = {boxes['level2_bottom']['xaxis'][0]:.6f} {boxes['level2_bottom']['xaxis'][1]:.6f} {boxes['level2_bottom']['xaxis'][2]:.6f}")
    print(f"tagging.l2_bottom.l2_bottom.yaxis = {boxes['level2_bottom']['yaxis'][0]:.6f} {boxes['level2_bottom']['yaxis'][1]:.6f} {boxes['level2_bottom']['yaxis'][2]:.6f}")
    print(f"tagging.l2_bottom.l2_bottom.zaxis = {boxes['level2_bottom']['zaxis'][0]:.6f} {boxes['level2_bottom']['zaxis'][1]:.6f} {boxes['level2_bottom']['zaxis'][2]:.6f}")
    
    # Level 2 Top Refinement
    print(f"\n# Level 2 Refinement Top (viscous sublayer, y⁺ ~ {yp2:.0f})")
    print(f"tagging.l2_top.type = GeometryRefinement")
    print(f"tagging.l2_top.shapes = l2_top")
    print(f"tagging.l2_top.max_level = 2")
    print(f"\ntagging.l2_top.l2_top.type = box")
    print(f"tagging.l2_top.l2_top.origin = {boxes['level2_top']['origin'][0]:.6f} {boxes['level2_top']['origin'][1]:.6f} {boxes['level2_top']['origin'][2]:.6f}")
    print(f"tagging.l2_top.l2_top.xaxis = {boxes['level2_top']['xaxis'][0]:.6f} {boxes['level2_top']['xaxis'][1]:.6f} {boxes['level2_top']['xaxis'][2]:.6f}")
    print(f"tagging.l2_top.l2_top.yaxis = {boxes['level2_top']['yaxis'][0]:.6f} {boxes['level2_top']['yaxis'][1]:.6f} {boxes['level2_top']['yaxis'][2]:.6f}")
    print(f"tagging.l2_top.l2_top.zaxis = {boxes['level2_top']['zaxis'][0]:.6f} {boxes['level2_top']['zaxis'][1]:.6f} {boxes['level2_top']['zaxis'][2]:.6f}")
    
    print(f"\namr.max_level = 2  # Max AMR level in hierarchy")


def main():
    parser = argparse.ArgumentParser(
        description="Setup domain parameters for turbulent channel flow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python case_setup.py
  python case_setup.py --Nx 512 --IB
  python case_setup.py --delta 0.005 --Nx 256
  python case_setup.py --blocking-factor 8 --Nx 512
        """
    )
    
    parser.add_argument('--delta', type=float, default=0.005,
                        help='Channel half width in meters (default: 0.005)')
    parser.add_argument('--Nx', type=int, default=384,
                        help='Number of cells in x-direction (default: 384)')
    parser.add_argument('--IB', action='store_true',
                        help='Use immersed boundary (extends domain in z)')
    parser.add_argument('--blocking-factor', type=int, default=8,
                        help='AMR blocking factor (default: 8)')
    parser.add_argument('--Re', type=int, default=180,
                        help='Stress Reynolds number (Re_tau) - options: 180, 395, 934 (default: 180)')
    
    args = parser.parse_args()
    
    params = domain_and_flow(
        delta=args.delta,
        Nx=args.Nx,
        IB=args.IB,
        blocking_factor=args.blocking_factor,
        Re=args.Re
    )
    
    print_config(params)


if __name__ == "__main__":
    main()
