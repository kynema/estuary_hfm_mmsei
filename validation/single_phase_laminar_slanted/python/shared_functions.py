#!/usr/bin/env python3
"""
Shared utility functions for channel post-processing
"""

import os
import re
import numpy as np


def parse_inp_for_channel_geometry(inp_file, return_dp_dx=False):
    """Extract channel geometry and physics parameters from base-*.inp file
    
    Args:
        inp_file: Path to .inp file
        return_dp_dx: If True, returns dp_dx and segment info; if False, returns compact form
    
    Returns:
        If return_dp_dx=True:
            tuple: (z_s, H, theta_deg, x_hi, z_hi, Umax, rho, mu, x_start, x_end, 
                   z_start, z_end, L_channel, dp_dx)
        If return_dp_dx=False:
            tuple: (z_s, H, theta_deg, x_hi, z_hi, Umax, rho, mu)
        
        Raises RuntimeError if required parameters cannot be found
    """
    
    with open(inp_file, 'r') as f:
        content = f.read()
    
    missing_params = []
    
    # Extract segment start and end points
    start_match = re.search(r'ChannelBuilder\.s1\.segment_start_point\s*=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)', content)
    end_match = re.search(r'ChannelBuilder\.s1\.segment_end_point\s*=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)', content)
    height_match = re.search(r'ChannelBuilder\.s1\.height_start\s*=\s*([\d.\-]+)', content)
    
    if not start_match or not end_match:
        missing_params.append('ChannelBuilder.s1.segment_start_point or segment_end_point')
    else:
        x_start = float(start_match.group(1))
        y_start = float(start_match.group(2))
        z_start = float(start_match.group(3))
        x_end = float(end_match.group(1))
        y_end = float(end_match.group(2))
        z_end = float(end_match.group(3))
        H = float(height_match.group(1)) if height_match else None
        
        # z_s is the z-coordinate where the centerline crosses x=0
        # For a line from (x_start, z_start) to (x_end, z_end):
        # z(x) = z_start + (z_end - z_start)/(x_end - x_start) * (x - x_start)
        # At x=0: z_s = z_start + (z_end - z_start)/(x_end - x_start) * (0 - x_start)
        if x_end != x_start:
            z_s = z_start - (z_end - z_start) * x_start / (x_end - x_start)
        else:
            z_s = z_start
        
        # Calculate theta from segment orientation
        dx = x_end - x_start
        dz = z_end - z_start
        theta_rad = np.arctan2(dz, dx) if dx != 0 else 0.0
        theta_deg = np.rad2deg(theta_rad)
        
        # Channel length along centerline
        L_channel = np.sqrt(dx**2 + dz**2)
        
        # Get x_hi (end x coordinate for segment s1)
        x_hi = x_end
    
    # Extract domain bounds from prob_hi
    prob_hi_match = re.search(r'geometry\.prob_hi\s*=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)', content)
    if prob_hi_match:
        x_hi = float(prob_hi_match.group(1))  # prob_hi[0]
        z_hi = float(prob_hi_match.group(3))  # prob_hi[2]
    else:
        z_hi = None
        missing_params.append('geometry.prob_hi')
    
    # Extract physics parameters
    # Umax from ChannelBuilder.s1.flow_speed
    flow_speed_match = re.search(r'ChannelBuilder\.s1\.flow_speed\s*=\s*([\d.\-eE]+)', content)
    if flow_speed_match:
        Umax = float(flow_speed_match.group(1))
    else:
        Umax = None
        missing_params.append('ChannelBuilder.s1.flow_speed')
    
    # Density from incflo.density
    density_match = re.search(r'incflo\.density\s*=\s*([\d.\-eE]+)', content)
    if density_match:
        rho = float(density_match.group(1))
    else:
        rho = None
        missing_params.append('incflo.density')
    
    # Viscosity from transport.viscosity
    viscosity_match = re.search(r'transport\.viscosity\s*=\s*([\d.\-eE]+)', content)
    if viscosity_match:
        mu = float(viscosity_match.group(1))
    else:
        mu = None
        missing_params.append('transport.viscosity')
    
    if missing_params:
        error_msg = f"Error: Could not parse required parameters from {os.path.basename(inp_file)}:\n"
        error_msg += f"  Missing: {', '.join(missing_params)}"
        raise RuntimeError(error_msg)
    
    # Optionally extract pressure gradient (dp/dx) from BodyForce.magnitude[0]
    if return_dp_dx:
        bodyforce_match = re.search(r'BodyForce\.magnitude\s*=\s*([\d.\-eE]+)\s+([\d.\-eE]+)\s+([\d.\-eE]+)', content)
        if bodyforce_match:
            dp_dx = float(bodyforce_match.group(1))  # First component of BodyForce
        else:
            dp_dx = None
            missing_params.append('BodyForce.magnitude')
        
        if missing_params:
            error_msg = f"Error: Could not parse required parameters from {os.path.basename(inp_file)}:\n"
            error_msg += f"  Missing: {', '.join(missing_params)}"
            raise RuntimeError(error_msg)
        
        return z_s, H, theta_deg, x_hi, z_hi, Umax, rho, mu, x_start, x_end, z_start, z_end, L_channel, dp_dx
    else:
        return z_s, H, theta_deg, x_hi, z_hi, Umax, rho, mu


def radial_distance_from_axis(x_coord, z_coord, z_s, theta):
    """Perpendicular distance from tilted centerline
    
    Args:
        x_coord: x-coordinate(s)
        z_coord: z-coordinate(s) 
        z_s: z-intercept of centerline (where it crosses x=0)
        theta: tilt angle in radians
    
    Returns:
        Perpendicular distance(s) as scalar or array
    """
    return np.abs((z_coord - z_s) * np.cos(theta) - x_coord * np.sin(theta))


def velocity_profile_analytical(r, U_max=1.0, H=1.0):
    """Parabolic Poiseuille profile for velocity perpendicular to centerline
    
    Args:
        r: Perpendicular distance from centerline
        U_max: Maximum velocity at centerline
        H: Channel height (perpendicular)
    
    Returns:
        Velocity as scalar or array, clipped to be non-negative
    """
    return U_max * np.maximum(0.0, 1.0 - (2.0 * r / H)**2)


def axial_velocity(u, w, theta):
    """Project 2D velocity components to direction along tilted flow axis
    
    Args:
        u: Velocity component in x-direction
        w: Velocity component in z-direction
        theta: Tilt angle in radians
    
    Returns:
        Axial velocity component (along tilted flow direction)
    """
    return u * np.cos(theta) + w * np.sin(theta)


def discover_case_variants(rootDir, mode_name, resolutions=None):
    """Auto-detect available drag variants in cases directory
    
    Args:
        rootDir: Root cases directory path
        mode_name: 'flat' or 'slanted'
        resolutions: List of resolutions to check (default: [32, 64, 128, 256, 512])
    
    Returns:
        Sorted list of available variants: ['og', 'temp', 'tf1']
    """
    if resolutions is None:
        resolutions = [32, 64, 128, 256, 512]
    
    available_variants = set()
    possible_variants = ['og', 'temp', 'tf1']
    
    for variant in possible_variants:
        for resolution in resolutions:
            case_name = f'{mode_name}-drag-{variant}-{resolution}'
            case_dir = os.path.join(rootDir, case_name)
            if os.path.exists(case_dir):
                plt_dirs = [d for d in os.listdir(case_dir) if d.startswith('plt')]
                if plt_dirs:
                    available_variants.add(variant)
                    break
    
    # Sort for consistent ordering: og, temp, tf1
    variant_order = ['og', 'temp', 'tf1']
    return [v for v in variant_order if v in available_variants]


def discover_case_variants_with_align(rootDir, mode_name, align='cf', resolutions=None):
    """Auto-detect available drag variants with alignment support for flat cases
    
    Args:
        rootDir: Root cases directory path
        mode_name: 'flat' or 'slanted'
        align: 'cf' or 'cc' (only for flat mode)
        resolutions: List of resolutions to check (default: [32, 64, 128, 256, 512])
    
    Returns:
        Sorted list of available variants: ['og', 'temp', 'tf1']
    """
    if resolutions is None:
        resolutions = [32, 64, 128, 256, 512]
    
    available_variants = set()
    possible_variants = ['og', 'temp', 'tf1']
    
    for variant in possible_variants:
        for resolution in resolutions:
            if mode_name == 'flat':
                case_name = f'{mode_name}-drag-{variant}-{align}-{resolution}'
            else:
                case_name = f'{mode_name}-drag-{variant}-{resolution}'
            case_dir = os.path.join(rootDir, case_name)
            if os.path.exists(case_dir):
                plt_dirs = [d for d in os.listdir(case_dir) if d.startswith('plt')]
                if plt_dirs:
                    available_variants.add(variant)
                    break
    
    # Sort for consistent ordering: og, temp, tf1
    variant_order = ['og', 'temp', 'tf1']
    return [v for v in variant_order if v in available_variants]


def get_analytical_solution_grid(x_arr, y_arr, z_arr, z_s_analytical, theta_rad, maxVelocity, 
                                  channelHeight, is_slanted=False):
    """Compute analytical velocity solution on a 3D grid
    
    Args:
        x_arr: 1D x-coordinate array
        y_arr: 1D y-coordinate array
        z_arr: 1D z-coordinate array
        z_s_analytical: z-intercept of channel centerline (where it crosses x=0)
        theta_rad: Tilt angle in radians
        maxVelocity: Maximum velocity at centerline
        channelHeight: Channel height (perpendicular)
        is_slanted: If True, use looser mask for 3-segment channels (default: False)
    
    Returns:
        tuple: (u_r, r, z_s_analytical) where:
            u_r: 3D array of analytical axial velocity
            r: 3D array of radial distances from centerline
            z_s_analytical: z-intercept of centerline
    """
    X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
    r = radial_distance_from_axis(X, Z, z_s_analytical, theta_rad)
    u_r = velocity_profile_analytical(r, maxVelocity, channelHeight)
    
    # Apply masking based on mode
    if is_slanted:
        # For slanted cases, use looser mask (r <= 1.2*H) to account for 3 segments
        channel_mask = r <= (1.2 * channelHeight)
    else:
        # For flat cases, mask strictly to channel interior (r <= H/2)
        channel_mask = r <= (channelHeight / 2.0)
    u_r[~channel_mask] = 0.0

    return u_r, r, z_s_analytical


def extract_centerline_z_s_from_inp(inp_file):
    """Extract centerline z-intercept from case-specific .inp file
    
    Reads ChannelBuilder.s1.segment_start_point to get z_s (z-coordinate at x=0)
    Format: segment_start_point = x_start y_center z_s
    
    Args:
        inp_file: Path to case-specific .inp file
    
    Returns:
        float: z_s value (z-coordinate of centerline at x=0)
        
    Raises:
        RuntimeError: If required parameter cannot be found
    """
    with open(inp_file, 'r') as f:
        content = f.read()
    
    # Extract segment_start_point: expects format "x y z" (three values)
    pattern = r'ChannelBuilder\.s1\.segment_start_point\s*=\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)'
    match = re.search(pattern, content)
    
    if match:
        x_start = float(match.group(1))
        y_center = float(match.group(2))
        z_s = float(match.group(3))
        return z_s
    
    raise RuntimeError(f'Could not extract ChannelBuilder.s1.segment_start_point from {inp_file}')


def get_case_specific_centerlines(rootDir, mode_name, align, valid_cases, base_z_s):
    """Get resolution-specific centerline z_s for each case
    
    For cc (cell-center) alignment, the centerline z_s is different for each mesh resolution
    because it depends on finding the nearest cell centers to target positions.
    Read each case's .inp file to get the correct z_s for that resolution.
    
    For cf (cell-face) alignment or slanted mode, use the base centerline for all cases.
    
    Args:
        rootDir: Root cases directory path
        mode_name: 'flat' or 'slanted'
        align: 'cf' or 'cc'
        valid_cases: List of case directory names
        base_z_s: Base centerline z_s value (fallback)
    
    Returns:
        list: z_s values, one per case (same length as valid_cases)
        
    Raises:
        RuntimeError: For cc alignment, if required .inp file is not found
    """
    z_s_case_specific = []
    
    if mode_name == 'flat' and align == 'cc':
        print('Reading resolution-specific centerline positions for cc alignment:')
        inp_filename = 'flat.inp'  # Flat cases use flat.inp
        for case_idx, case_name in enumerate(valid_cases):
            case_dir = os.path.join(rootDir, case_name)
            inp_file = os.path.join(case_dir, inp_filename)
            
            if os.path.exists(inp_file):
                try:
                    z_s_res = extract_centerline_z_s_from_inp(inp_file)
                    z_s_case_specific.append(z_s_res)
                    print(f'  {case_name}: z_s = {z_s_res:.6f}')
                except Exception as e:
                    raise RuntimeError(f'Error reading centerline from {inp_file}: {e}')
            else:
                raise RuntimeError(f'Required .inp file not found for cc alignment: {inp_file}')
    else:
        # For flat cf mode or slanted mode, use base centerline for all cases
        z_s_case_specific = [base_z_s] * len(valid_cases)
    
    return z_s_case_specific


def print_effective_height_correction(align, channelHeight, dp_dx_from_inp, rho_from_inp, 
                                      mu_from_inp, nx, dx):
    """Print effective height correction for flat channel with cell-face alignment
    
    For IB aligned to cell face: H_eff = H + h (cell-center offset above/below IB)
    BodyForce.magnitude[0] is -1/rho * dp/dx, so dp/dx = -BodyForce.magnitude[0] * rho
    U_max_eff = BodyForce.magnitude[0] * rho * H_eff^2 / (8*mu)
    
    Args:
        align: Grid alignment mode ('cf' or 'cc')
        channelHeight: Geometric channel height
        dp_dx_from_inp: Pressure gradient coefficient from BodyForce.magnitude[0]
        rho_from_inp: Fluid density
        mu_from_inp: Dynamic viscosity
        nx: List of grid resolutions in x-direction
        dx: List of cell sizes
    """
    if align == 'cf':
        print('\nEffective Height Correction (IB at cell face):')
        print(f'  Geometric height H = {channelHeight:.6f}')
        print(f'  BodyForce.magnitude[0] = {dp_dx_from_inp:.6f} (= -1/rho * dp/dx)')
        print(f'  rho = {rho_from_inp:.6f}, mu = {mu_from_inp:.6f}')
        for case_idx in range(len(nx)):
            h_cell = dx[case_idx]
            H_eff = channelHeight + h_cell
            U_max_eff = dp_dx_from_inp * rho_from_inp * (H_eff ** 2) / (8.0 * mu_from_inp)
            print(f'  nx={nx[case_idx]:3d}: h={h_cell:.6f}, H_eff={H_eff:.6f}, U_max_eff={U_max_eff:.6f}')
