#!/usr/bin/env python3
"""
Slanted Channel Post-Processing: Analytical vs Numerical Solutions

Compare numerical simulations to analytical solution for tilted parabolic pipe flow.
Analyze error convergence with grid refinement and wall boundary errors.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import yt
import warnings
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

warnings.filterwarnings('ignore')

print('Modules loaded')

# ================================================================================
# Helper functions
# ================================================================================
def parse_inp_for_channel_geometry(inp_file):
    """Extract channel geometry and physics parameters from base-*.inp file
    
    Returns:
        tuple: (z_s, H, theta_deg, x_hi, z_hi, Umax, rho, mu)
        Raises RuntimeError if required parameters cannot be found
    """
    import re
    
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
        x_start, y_start, z_start = float(start_match.group(1)), float(start_match.group(2)), float(start_match.group(3))
        x_end, y_end, z_end = float(end_match.group(1)), float(end_match.group(2)), float(end_match.group(3))
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
        
        # Get x_hi (end x coordinate for segment s1)
        x_hi = x_end
    
    # Extract domain bounds from prob_hi
    prob_hi_match = re.search(r'geometry\.prob_hi\s*=\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)', content)
    if prob_hi_match:
        x_hi = float(prob_hi_match.group(1))  # prob_hi[0]
        # y_hi = float(prob_hi_match.group(2))  # prob_hi[1] - not needed
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
    
    return z_s, H, theta_deg, x_hi, z_hi, Umax, rho, mu

# ================================================================================
# Parse command line arguments
# ================================================================================
parser = argparse.ArgumentParser(description='Post-process slanted/flat channel simulations')

mode_group = parser.add_mutually_exclusive_group(required=True)
mode_group.add_argument('--flat', action='store_true', help='Post-process flat channel cases')
mode_group.add_argument('--slanted', action='store_true', help='Post-process slanted channel cases')

parser.add_argument('--align', type=str, default='default', choices=['default', 'cf'], help='Grid alignment when theta=0: default (none) or cf (cell face)')
parser.add_argument('--nx_align', type=int, default=64, help='Reference nx for grid alignment (default: 64)')
parser.add_argument('--H', type=float, default=1.0, help='Channel height (default: 1.0)')
args = parser.parse_args()

# ================================================================================
# Physical parameters
# ================================================================================
# Set theta_deg based on mode
if args.flat:
    theta_deg = 0.0
    mode_name = 'flat'
else:  # slanted
    theta_deg = 45.0  # Default for slanted; will be overridden by parsing .inp
    mode_name = 'slanted'

# Define paths
rootDir = '/Users/dmontgo2/Documents/Kynema/estuary_hfm_mmsei/validation/single_phase_laminar_slanted/cases'
figureDir_base = '/Users/dmontgo2/Documents/Kynema/estuary_hfm_mmsei/validation/single_phase_laminar_slanted/figures'
base_inp = os.path.join(rootDir, 'base-flat-poiseuille.inp' if args.flat else 'base-slanted-poiseuille.inp')

# Figure directory for drag comparison (mode_name-drag-comparison)
figureDir = os.path.join(figureDir_base, f'{mode_name}-drag-comparison')
os.makedirs(figureDir, exist_ok=True)

# Initialize parameters
theta_rad = np.deg2rad(theta_deg)
channelHeight = args.H
maxVelocity = 1.0

# Read channel geometry and physics parameters from base .inp file
if not os.path.exists(base_inp):
    raise RuntimeError(f"Error: Base input file not found: {base_inp}")

try:
    z_s_analytical, H_from_inp, theta_from_inp, x_hi_from_inp, z_hi_from_inp, Umax_from_inp, rho_from_inp, mu_from_inp = parse_inp_for_channel_geometry(base_inp)
    
    print(f'Read channel configuration from {os.path.basename(base_inp)}:')
    print(f'  Channel centerline z_s = {z_s_analytical:.6f}')
    print(f'  Channel height H = {H_from_inp:.6f}')
    print(f'  Channel angle theta = {theta_from_inp:.2f}°')
    print(f'  Domain: x_hi={x_hi_from_inp:.4f}, z_hi={z_hi_from_inp:.4f}')
    print(f'  Physics: Umax={Umax_from_inp:.4f}, rho={rho_from_inp:.4f}, mu={mu_from_inp:.4f}')
    
    # Update parameters from file
    channelHeight = H_from_inp
    theta_deg = theta_from_inp
    theta_rad = np.deg2rad(theta_deg)
    maxVelocity = Umax_from_inp
    
except RuntimeError as e:
    print(e)
    raise

# ================================================================================
# Helper functions for analytical solution
# ================================================================================
def radial_distance_from_axis(x_coord, z_coord, z_s, theta):
    """Perpendicular distance from tilted centerline"""
    return np.abs((z_coord - z_s) * np.cos(theta) - x_coord * np.sin(theta))

def velocity_profile_analytical(r, U_max=1.0, H=1.0):
    """Parabolic profile: u(r) = U_max * (1 - (2r/H)^2)"""
    return U_max * np.maximum(0.0, 1.0 - (2.0 * r / H)**2)

def axial_velocity(u, w, theta):
    """Compute axial velocity along the tilted flow direction"""
    return u * np.cos(theta) + w * np.sin(theta)

# ================================================================================
# Case discovery and loading
# ================================================================================
# Build case specifications to load both og and temp variants
resolutions = [32, 64, 128, 256, 512]
drag_variants = ['og', 'temp']

case_paths_dict = {'og': [], 'temp': []}  # Organize by drag variant
valid_cases_dict = {'og': [], 'temp': []}

for variant in drag_variants:
    for resolution in resolutions:
        case_name = f'{mode_name}-drag-{variant}-{resolution}'
        case_dir = os.path.join(rootDir, case_name)
        if os.path.exists(case_dir):
            plt_dirs = sorted([d for d in os.listdir(case_dir) if d.startswith('plt')])
            if plt_dirs:
                case_paths_dict[variant].append(os.path.join(case_dir, plt_dirs[-1]))
                valid_cases_dict[variant].append(case_name)
                print(f'{case_name}: {plt_dirs[-1]}')

nCases_og = len(valid_cases_dict['og'])
nCases_temp = len(valid_cases_dict['temp'])
print(f'\nFound {nCases_og} OG cases and {nCases_temp} TEMP cases')

print(f'\nAnalytical Solution Parameters:')
print(f'  θ = {theta_deg:.2f}°')
print(f'  H = {channelHeight:.6f}')
print(f'  Channel centerline z_s = {z_s_analytical:.6f}')
print(f'  Alignment: {args.align if np.abs(theta_deg) < 0.01 else "N/A (slanted)"}')

# ================================================================================
# Load datasets using yt for both drag variants
# ================================================================================
# Initialize data structures for each variant
data_by_variant = {}

for variant in drag_variants:
    data_by_variant[variant] = {
        'ds_list': [], 'x': [], 'y': [], 'z': [],
        'nx': [], 'ny': [], 'nz': [],
        'dx': [], 'dy': [], 'dz': [],
        'u': [], 'v': [], 'w': [], 'p': [],
        'domain_bounds': [],
        'ur_exact': [], 'r_dist': [], 'centerline_shift': []
    }

# Load data for each variant
for variant in drag_variants:
    print(f'\nLoading {variant.upper()} variant cases...')
    for case_idx, case_path in enumerate(case_paths_dict[variant]):
        case_name = valid_cases_dict[variant][case_idx]
        print(f'Loading {case_name}...')
        try:
            # Load with yt
            ds_yt = yt.load(case_path)
            data_by_variant[variant]['ds_list'].append(ds_yt)
            
            # Get grid information
            dims = ds_yt.domain_dimensions
            domain_left = ds_yt.domain_left_edge.d
            domain_right = ds_yt.domain_right_edge.d
            data_by_variant[variant]['domain_bounds'].append((domain_left, domain_right))
            
            print(f'  Dimensions: {dims}')
            
            # Use covering_grid to get structured data directly
            cube = ds_yt.covering_grid(level=0, left_edge=ds_yt.domain_left_edge, dims=dims)
            
            # Extract coordinate arrays (automatically cell-centered)
            x_arr = cube["index", "x"].d[:, 0, 0]  # Extract 1D by taking first slice in y,z
            y_arr = cube["index", "y"].d[0, :, 0]  # Extract 1D by taking first slice in x,z
            z_arr = cube["index", "z"].d[0, 0, :]  # Extract 1D by taking first slice in x,y
            
            data_by_variant[variant]['x'].append(x_arr)
            data_by_variant[variant]['y'].append(y_arr)
            data_by_variant[variant]['z'].append(z_arr)
            
            data_by_variant[variant]['nx'].append(dims[0])
            data_by_variant[variant]['ny'].append(dims[1])
            data_by_variant[variant]['nz'].append(dims[2])
            
            data_by_variant[variant]['dx'].append(x_arr[1] - x_arr[0] if len(x_arr) > 1 else 0)
            data_by_variant[variant]['dy'].append(y_arr[1] - y_arr[0] if len(y_arr) > 1 else 0)
            data_by_variant[variant]['dz'].append(z_arr[1] - z_arr[0] if len(z_arr) > 1 else 0)
            
            # Extract 3D field data
            u_data_3d = cube["velocityx"].d.copy()
            v_data_3d = cube["velocityy"].d.copy()
            w_data_3d = cube["velocityz"].d.copy()
            p_data_3d = cube["p"].d.copy()
            
            # Mask numerical data outside analysis region immediately after reading
            # Create radial distance field for masking
            X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
            r_mask = radial_distance_from_axis(X, Z, z_s_analytical, theta_rad)
            
            # Determine threshold based on mode
            if args.slanted:
                mask_outside = r_mask > (1.2 * channelHeight / 2.0)
            else:
                mask_outside = r_mask > (1.2 * channelHeight / 2.0)
            
            # Apply mask - explicitly set masked values to 0
            u_data_3d[mask_outside] = 0.0
            v_data_3d[mask_outside] = 0.0
            w_data_3d[mask_outside] = 0.0
            p_data_3d[mask_outside] = 0.0
            
            data_by_variant[variant]['u'].append(u_data_3d)
            data_by_variant[variant]['v'].append(v_data_3d)
            data_by_variant[variant]['w'].append(w_data_3d)
            data_by_variant[variant]['p'].append(p_data_3d)
            
            print(f'  Grid spacing: dx={data_by_variant[variant]["dx"][-1]:.6f}, dy={data_by_variant[variant]["dy"][-1]:.6f}, dz={data_by_variant[variant]["dz"][-1]:.6f}')
            print(f'  Velocity range: u=[{u_data_3d.min():.4f}, {u_data_3d.max():.4f}]')
            
        except Exception as e:
            print(f'  Error: {e}')
            import traceback
            traceback.print_exc()

# Compute analytical solutions for each variant
def get_analytical_solution_grid(x_arr, y_arr, z_arr):
    X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
    r = radial_distance_from_axis(X, Z, z_s_analytical, theta_rad)
    u_r = velocity_profile_analytical(r, maxVelocity, channelHeight)
    
    # For slanted cases, use looser mask (r <= 1.2*H) to account for 3 segments
    # For flat cases, mask strictly to channel interior (r <= H/2)
    if args.slanted:
        channel_mask = r <= (1.2 * channelHeight)
    else:
        channel_mask = r <= (channelHeight / 2.0)
    u_r[~channel_mask] = 0.0

    return u_r, r, z_s_analytical

for variant in drag_variants:
    for case_idx in range(len(data_by_variant[variant]['x'])):
        ur_a, r, z_s = get_analytical_solution_grid(
            data_by_variant[variant]['x'][case_idx],
            data_by_variant[variant]['y'][case_idx],
            data_by_variant[variant]['z'][case_idx]
        )
        data_by_variant[variant]['ur_exact'].append(ur_a)
        data_by_variant[variant]['r_dist'].append(r)
        data_by_variant[variant]['centerline_shift'].append(z_s)

print('\nAnalytical solutions computed for all variants')

# ================================================================================
# Calculate Errors for Both Variants
# ================================================================================
errors_by_variant = {}

for variant in drag_variants:
    print(f'\nComputing errors for {variant.upper()} variant:')
    errors_by_variant[variant] = {
        'error_max': [],
        'error_l2': [],
        'error_max_wall': [],
        'error_l2_wall': [],
        'nx_array': [],
        'cell_size_array': []
    }
    
    nCases_variant = len(data_by_variant[variant]['x'])
    for case_idx in range(nCases_variant):
        case_name = valid_cases_dict[variant][case_idx]
        j_mid = data_by_variant[variant]['ny'][case_idx] // 2
        
        # Extract 2D slices at midpoint
        u_2d = data_by_variant[variant]['u'][case_idx][:, j_mid, :]
        w_2d = data_by_variant[variant]['w'][case_idx][:, j_mid, :]
        ur_exact_2d = data_by_variant[variant]['ur_exact'][case_idx][:, j_mid, :]
        r_dist_2d = data_by_variant[variant]['r_dist'][case_idx][:, j_mid, :]
        
        # Calculate axial velocity: component along the tilted flow direction
        u_axial_num = axial_velocity(u_2d, w_2d, theta_rad)
        u_axial_ana = ur_exact_2d
        
        # Compute axial velocity errors
        u_axial_err = u_axial_num - u_axial_ana
        
        # For slanted cases, use looser mask (r <= 1.2 * H/2) to account for 3 segments
        # For flat cases, use channel interior only (r <= H/2)
        if args.slanted:
            channel_mask = r_dist_2d <= (1.2 * channelHeight / 2.0)
        else:
            channel_mask = r_dist_2d <= (channelHeight / 2.0)
        
        u_axial_max = np.max(np.abs(u_axial_err[channel_mask])) if np.any(channel_mask) else 0.0
        u_axial_l2 = np.sqrt(np.sum(u_axial_err[channel_mask]**2)) / np.sqrt(np.sum(channel_mask)) if np.any(channel_mask) else 0.0
        
        # Near-wall error (quarter of channel height from wall)
        wall_distance = channelHeight / 4.0
        near_wall_mask = (r_dist_2d > (channelHeight / 2.0 - wall_distance)) & (r_dist_2d <= (channelHeight / 2.0))
        u_axial_max_wall = np.max(np.abs(u_axial_err[near_wall_mask])) if np.any(near_wall_mask) else 0.0
        u_axial_l2_wall = np.sqrt(np.sum(u_axial_err[near_wall_mask]**2)) / np.sqrt(np.sum(near_wall_mask)) if np.any(near_wall_mask) else 0.0
        
        errors_by_variant[variant]['error_max'].append(u_axial_max)
        errors_by_variant[variant]['error_l2'].append(u_axial_l2)
        errors_by_variant[variant]['error_max_wall'].append(u_axial_max_wall)
        errors_by_variant[variant]['error_l2_wall'].append(u_axial_l2_wall)
        errors_by_variant[variant]['nx_array'].append(data_by_variant[variant]['nx'][case_idx])
        errors_by_variant[variant]['cell_size_array'].append(data_by_variant[variant]['dx'][case_idx])
        
        print(f'  {case_name} (y-index {j_mid}):')
        print(f'    Overall - Max: {u_axial_max:.6e}, L2: {u_axial_l2:.6e}')
        print(f'    Wall    - Max: {u_axial_max_wall:.6e}, L2: {u_axial_l2_wall:.6e}')

# Convert to numpy arrays for plotting
for variant in drag_variants:
    errors_by_variant[variant]['error_max'] = np.array(errors_by_variant[variant]['error_max'])
    errors_by_variant[variant]['error_l2'] = np.array(errors_by_variant[variant]['error_l2'])
    errors_by_variant[variant]['error_max_wall'] = np.array(errors_by_variant[variant]['error_max_wall'])
    errors_by_variant[variant]['error_l2_wall'] = np.array(errors_by_variant[variant]['error_l2_wall'])
    errors_by_variant[variant]['cell_size_array'] = np.array(errors_by_variant[variant]['cell_size_array'])
    errors_by_variant[variant]['nx_array'] = np.array(errors_by_variant[variant]['nx_array'])

# ================================================================================
# Error Convergence vs Cell Size (Drag Comparison)
# ================================================================================
fig, ax = plt.subplots(figsize=(10, 7))

colors = {'og': 'C0', 'temp': 'C1'}
markers = {'og': 'o', 'temp': 's'}

for variant in drag_variants:
    cell_size_array = errors_by_variant[variant]['cell_size_array']
    u_axial_error_max = errors_by_variant[variant]['error_max']
    u_axial_error_l2 = errors_by_variant[variant]['error_l2']
    
    ax.loglog(cell_size_array, u_axial_error_max, marker=markers[variant], linestyle='-', linewidth=2, 
              markersize=10, label=f'{variant.upper()} - Max Error', color=colors[variant])
    ax.loglog(cell_size_array, u_axial_error_l2, marker=markers[variant], linestyle='--', linewidth=2, 
              markersize=8, label=f'{variant.upper()} - L2 Error', color=colors[variant])

# Add reference slopes using combined data
all_cell_sizes = np.concatenate([errors_by_variant[v]['cell_size_array'] for v in drag_variants])
all_errors = np.concatenate([errors_by_variant[v]['error_max'] for v in drag_variants])

if len(all_cell_sizes) > 1:
    cell_trend = np.logspace(np.log10(all_cell_sizes.min()), np.log10(all_cell_sizes.max()), 50)
    
    # Scale reference lines to pass through mean of max error data
    mean_error_max = np.mean(all_errors)
    mean_cell_size = np.mean(all_cell_sizes)
    
    # O(h) reference: error ~ C * h
    C_slope1 = 0.85 * mean_error_max / (mean_cell_size ** 1.0)
    ax.loglog(cell_trend, C_slope1 * cell_trend**1.0, ':', alpha=0.6, linewidth=2, color='gray', label='Slope: 1')
    
    # O(h^1/2) reference: error ~ C * h^0.5
    C_slope05 = 0.85 * mean_error_max / (mean_cell_size ** 0.5)
    ax.loglog(cell_trend, C_slope05 * cell_trend**0.5, '-.', alpha=0.6, linewidth=2, color='gray', label='Slope: 0.5')

ax.set_xlabel('Cell Size (h)', fontsize=14)
ax.set_ylabel('Axial Velocity Error', fontsize=14)
ax.set_title(f'Drag Forcing Comparison: {mode_name.title()} Channel Convergence', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(f'{figureDir}/error_convergence.png', dpi=150)
plt.show()

print('\nPost-processing complete!')

print('Post-processing complete!')
