"""
Setup domain parameters for turbulent channel flow validation case.
Based on DNS configuration at Re_tau = 180

This script calculates:
1. Domain bounds (geometry.prob_lo/hi)
2. Grid resolution (amr.n_cell) - divisible by blocking factor
3. Cell spacing maintained isotropic in all directions
"""

import argparse
import numpy as np


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


def get_refinement_boxes(delta, prob_lo, prob_hi, n_cell, density, mu, u_tau, IB=False, ref_ratio=2):
    """
    Define AMR refinement boxes aligned with the base grid cells and verify cell spacing.
    """
    # Extract domain dimensions and grid resolution
    Lx = prob_hi[0] - prob_lo[0]
    Ly = prob_hi[1] - prob_lo[1]
    Lz = prob_hi[2] - prob_lo[2]
    Nz = n_cell[2]

    # Viscous Length Scale
    nu = mu / density
    l_nu = nu / u_tau

    # Grid Spacing at each AMR level
    dz_base = Lz / Nz
    dz_l1 = dz_base / ref_ratio
    dz_l2 = dz_l1 / ref_ratio

    # Convert cell spacings to wall units (dz+)
    dz_plus_base = dz_base / l_nu
    dz_plus_l1 = dz_l1 / l_nu
    dz_plus_l2 = dz_l2 / l_nu

    # First cell center off the wall (y+ of first node)
    z1_plus_base = 0.5 * dz_plus_base
    z1_plus_l1 = 0.5 * dz_plus_l1
    z1_plus_l2 = 0.5 * dz_plus_l2

    # Target Box Extents (Where we want refinement to reach)
    yplus_target_l1 = 40.0  # Captures up to y+ ~ 40
    yplus_target_l2 = 5.0   # Captures up to y+ ~ 5

    # Convert target y+ to physical distance
    z_extent_l1 = yplus_target_l1 * l_nu
    z_extent_l2 = yplus_target_l2 * l_nu

    # Align to base grid cell spacing (dz_base) to ensure boxes align with grid
    cells_l1 = max(1, int(np.ceil(z_extent_l1 / dz_base)))
    cells_l2 = max(1, int(np.ceil(z_extent_l2 / dz_base)))
    z_extent_l1 = cells_l1 * dz_base
    z_extent_l2 = cells_l2 * dz_base

    # Limit to half/quarter channel height to avoid overlapping with center
    z_extent_l1 = min(z_extent_l1, 0.5 * delta)
    z_extent_l2 = min(z_extent_l2, 0.25 * delta)

    # Actual physical bounds of the boxes in wall units
    l1_box_extent_plus = z_extent_l1 / l_nu
    l2_box_extent_plus = z_extent_l2 / l_nu

    # Check adequcate resolution for first cell center off the wall
    if z1_plus_l2 >= 1.0:
        print(f"Warning: First cell center off the wall at Level 2 (z1⁺ = {z1_plus_l2:.2f}) > 1.0!")
        print(f"   Consider increasing Nz (base grid) or using a higher ref_ratio.")

    boxes = {
        'level1_bottom': {
            'origin': [prob_lo[0], prob_lo[1], prob_lo[2]],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_l1],
        },
        'level1_top': {
            'origin': [prob_lo[0], prob_lo[1], prob_hi[2] - z_extent_l1],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_l1],
        },
        'level2_bottom': {
            'origin': [prob_lo[0], prob_lo[1], prob_lo[2]],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_l2],
        },
        'level2_top': {
            'origin': [prob_lo[0], prob_lo[1], prob_hi[2] - z_extent_l2],
            'xaxis': [Lx, 0.0, 0.0],
            'yaxis': [0.0, Ly, 0.0],
            'zaxis': [0.0, 0.0, z_extent_l2],
        },
        'viscous_length_scale': l_nu,
        'dz_base': dz_base,
        'dz_l1': dz_l1,
        'dz_l2': dz_l2,
        'dz_plus_base': dz_plus_base,
        'dz_plus_l1': dz_plus_l1,
        'dz_plus_l2': dz_plus_l2,
        'z1_plus_base': z1_plus_base,
        'z1_plus_l1': z1_plus_l1,
        'z1_plus_l2': z1_plus_l2,
        'l1_box_extent_plus': l1_box_extent_plus,
        'l2_box_extent_plus': l2_box_extent_plus,
    }
    return boxes



def domain_and_flow(delta=0.005, Nx=384, IB=False, blocking_factor=4, Re=180):
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
    
    # Recalculate cell spacings after blocking factor adjustment (now slightly different)
    dx = Lx / Nx
    dy = Ly / Ny
    dz = Lz / Nz

    # Check that Lz ~= 2.0δ (no IB) or 2.08δ (with IB)
    if IB:
        assert np.isclose(Lz, Lz_ib, rtol=1e-3), f"Lz={Lz:.6f} m not close to expected {Lz_ib:.6f} m for IB"
    else:
        assert np.isclose(Lz, Lz_no_ib, rtol=1e-3), f"Lz={Lz:.6f} m not close to expected {Lz_no_ib:.6f} m for no IB"
    
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
    n_cell = [Nx, Ny, Nz]
    
    # Get flow characteristics for specified stress Reynolds number
    density = 0.468793
    mu = 3.57816e-5  # Dynamic viscosity (air at 750K)
    dpdx = get_pressure_gradient(Re)
    tau_w = - delta * dpdx # Wall shear stress
    u_tau = np.sqrt(tau_w / density)  # Friction velocity
    u_tau = Re * mu / (density * delta)  # Friction velocity
    Re_tau = u_tau * delta / (mu / density)  # Re_tau calculated
    t_star = delta / u_tau  # Time scale based on friction velocity and channel half-width

    # Body force acceleration: F = |dp/dx| / rho (units: m/s^2)
    body_force = -dpdx / density
    
    # Get refinement boxes (using computed physical parameters)
    refinement_boxes = get_refinement_boxes(delta, prob_lo, prob_hi, n_cell, density, mu, u_tau, IB=IB)
    
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
        'dy': dy,
        'dz': dz,
        'prob_lo': prob_lo,
        'prob_hi': prob_hi,
        'Re': Re_tau,
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
    print(f"  Viscous length scale (l_ν):    {config['refinement_boxes']['viscous_length_scale']:.6e} m")
    print(f"  Density:                       {config['density']} kg/m^3")
    print(f"  Dynamic viscosity:             {config['mu']:.5e} Pa.s")
    
    
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
    print(f"  n_cell.L0 = {config['Nx']} {config['Ny']} {config['Nz']}")
    print(f"  n_cell.L1 = {config['Nx']*2} {config['Ny']*2} {config['Nz']*2}")
    print(f"  n_cell.L2 = {config['Nx']*4} {config['Ny']*4} {config['Nz']*4}")
    print(f"  Cell spacing (L0): dx = {config['dx']:.8f} m, dy = {config['dy']:.8f} m, dz = {config['dz']:.8f} m")
    
    # Grid resolution analysis
    boxes = config['refinement_boxes']
    print(f"  L0: dz⁺ = {boxes['dz_plus_base']:.2f}")
    print(f"  L1: dz⁺ = {boxes['dz_plus_l1']:.2f}")
    print(f"  L2: dz⁺ = {boxes['dz_plus_l2']:.2f}" \
          f"  z1+ = {boxes['z1_plus_l2']:.2f}")
    
    print("\n" + "="*60)
    fname = f"turbulent-flat-re-{config['Re']}"
    if config['IB']:
        fname += "-IB"
    print(f"Parameters for {fname}.inp file:")
    print("="*60)
    print(f"\ntime.stop_time = {20*config['t_star']:.6f}  # ~20 flow-through time")
    print(f"#time.stop_time = {30*config['t_star']:.6f}  # ~30 flow-through time for statistics")

    print("\n#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#" \
          "\n#              GEOMETRY                 #" \
          "\n#.......................................#")
    print(f"geometry.prob_lo = {lo[0]:.1f} {lo[1]:.6f} {lo[2]:.6f}")
    print(f"geometry.prob_hi = {hi[0]:.6f} {hi[1]:.6f} {hi[2]:.6f}")

    print("\n#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#" \
          "\n#               PHYSICS                 #" \
          "\n#.......................................#")
    print(f"BodyForce.magnitude = {config['body_force']:.16f} 0.0 0.0 # Force acceleration in m/s^2")

    print("\n#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#" \
          "\n#       CHANNEL FLOW PARAMETERS         #" \
          "\n#.......................................#")
    print(f"ChannelFlow.re_tau = {config['Re']}")
    
    print("\n#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#" \
          "\n#        ADAPTIVE MESH REFINEMENT       #" \
          "\n#.......................................#")
    # Extents of refinement boxes in wall units
    l1_box_extent_plus = boxes['l1_box_extent_plus']
    l2_box_extent_plus = boxes['l2_box_extent_plus']
    print(f"amr.n_cell = {config['Nx']} {config['Ny']} {config['Nz']}")
    print(f"amr.blocking_factor = {config['blocking_factor']}")
    print(f"amr.max_level = 2")
    print(f"\ntagging.labels = l1 l2")
    
    # Level 1 Refinement (buffer layer)
    print(f"\n# Level 1 Refinement (buffer layer, z⁺ ~ {l1_box_extent_plus:.1f})")
    print(f"tagging.l1.type = GeometryRefinement")
    print(f"tagging.l1.shapes = bottom top")
    print(f"tagging.l1.max_level = 1")
    print(f"\ntagging.l1.bottom.type = box")
    print(f"tagging.l1.bottom.origin = {boxes['level1_bottom']['origin'][0]:.6f} {boxes['level1_bottom']['origin'][1]:.6f} {boxes['level1_bottom']['origin'][2]:.6f}")
    print(f"tagging.l1.bottom.xaxis = {boxes['level1_bottom']['xaxis'][0]:.6f} {boxes['level1_bottom']['xaxis'][1]:.6f} {boxes['level1_bottom']['xaxis'][2]:.6f}")
    print(f"tagging.l1.bottom.yaxis = {boxes['level1_bottom']['yaxis'][0]:.6f} {boxes['level1_bottom']['yaxis'][1]:.6f} {boxes['level1_bottom']['yaxis'][2]:.6f}")
    print(f"tagging.l1.bottom.zaxis = {boxes['level1_bottom']['zaxis'][0]:.6f} {boxes['level1_bottom']['zaxis'][1]:.6f} {boxes['level1_bottom']['zaxis'][2]:.6f}")
    
    print(f"\ntagging.l1.top.type = box")
    print(f"tagging.l1.top.origin = {boxes['level1_top']['origin'][0]:.6f} {boxes['level1_top']['origin'][1]:.6f} {boxes['level1_top']['origin'][2]:.6f}")
    print(f"tagging.l1.top.xaxis = {boxes['level1_top']['xaxis'][0]:.6f} {boxes['level1_top']['xaxis'][1]:.6f} {boxes['level1_top']['xaxis'][2]:.6f}")
    print(f"tagging.l1.top.yaxis = {boxes['level1_top']['yaxis'][0]:.6f} {boxes['level1_top']['yaxis'][1]:.6f} {boxes['level1_top']['yaxis'][2]:.6f}")
    print(f"tagging.l1.top.zaxis = {boxes['level1_top']['zaxis'][0]:.6f} {boxes['level1_top']['zaxis'][1]:.6f} {boxes['level1_top']['zaxis'][2]:.6f}")
    
    # Level 2 Refinement (viscous sublayer)
    print(f"\n# Level 2 Refinement (viscous sublayer, z⁺ ~ {l2_box_extent_plus:.1f})")
    print(f"tagging.l2.type = GeometryRefinement")
    print(f"tagging.l2.shapes = bottom top")
    print(f"tagging.l2.max_level = 2")
    print(f"\ntagging.l2.bottom.type = box")
    print(f"tagging.l2.bottom.origin = {boxes['level2_bottom']['origin'][0]:.6f} {boxes['level2_bottom']['origin'][1]:.6f} {boxes['level2_bottom']['origin'][2]:.6f}")
    print(f"tagging.l2.bottom.xaxis = {boxes['level2_bottom']['xaxis'][0]:.6f} {boxes['level2_bottom']['xaxis'][1]:.6f} {boxes['level2_bottom']['xaxis'][2]:.6f}")
    print(f"tagging.l2.bottom.yaxis = {boxes['level2_bottom']['yaxis'][0]:.6f} {boxes['level2_bottom']['yaxis'][1]:.6f} {boxes['level2_bottom']['yaxis'][2]:.6f}")
    print(f"tagging.l2.bottom.zaxis = {boxes['level2_bottom']['zaxis'][0]:.6f} {boxes['level2_bottom']['zaxis'][1]:.6f} {boxes['level2_bottom']['zaxis'][2]:.6f}")
    
    print(f"\ntagging.l2.top.type = box")
    print(f"tagging.l2.top.origin = {boxes['level2_top']['origin'][0]:.6f} {boxes['level2_top']['origin'][1]:.6f} {boxes['level2_top']['origin'][2]:.6f}")
    print(f"tagging.l2.top.xaxis = {boxes['level2_top']['xaxis'][0]:.6f} {boxes['level2_top']['xaxis'][1]:.6f} {boxes['level2_top']['xaxis'][2]:.6f}")
    print(f"tagging.l2.top.yaxis = {boxes['level2_top']['yaxis'][0]:.6f} {boxes['level2_top']['yaxis'][1]:.6f} {boxes['level2_top']['yaxis'][2]:.6f}")
    print(f"tagging.l2.top.zaxis = {boxes['level2_top']['zaxis'][0]:.6f} {boxes['level2_top']['zaxis'][1]:.6f} {boxes['level2_top']['zaxis'][2]:.6f}")


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
