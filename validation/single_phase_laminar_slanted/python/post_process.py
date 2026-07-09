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

drag_group = parser.add_mutually_exclusive_group(required=False)
drag_group.add_argument('--og', action='store_true', help='Use original drag forcing variant')
drag_group.add_argument('--temp', action='store_true', help='Use temporal drag forcing variant (default)')

parser.add_argument('--align', type=str, default='default', choices=['default', 'cf'], help='Grid alignment when theta=0: default (none) or cf (cell face)')
parser.add_argument('--nx_align', type=int, default=64, help='Reference nx for grid alignment (default: 64)')
parser.add_argument('--H', type=float, default=1.0, help='Channel height (default: 1.0)')
args = parser.parse_args()

# Determine drag variant (default to 'temp' if not specified)
drag_variant = 'og' if args.og else 'temp'

# ================================================================================
# Physical parameters
# ================================================================================
# Set theta_deg based on mode
if args.flat:
    theta_deg = 0.0
    case_prefix = f'flat-drag-{drag_variant}'
else:  # slanted
    theta_deg = 45.0  # Default for slanted; will be overridden by parsing .inp
    case_prefix = f'slanted-drag-{drag_variant}'

# Define paths
rootDir = '/Users/dmontgo2/Documents/Kynema/estuary_hfm_mmsei/validation/single_phase_laminar_slanted/cases'
figureDir_base = '/Users/dmontgo2/Documents/Kynema/estuary_hfm_mmsei/validation/single_phase_laminar_slanted/figures'
base_inp = os.path.join(rootDir, 'base-flat-poiseuille.inp' if args.flat else 'base-slanted-poiseuille.inp')

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

# Create figure output directory with drag variant subdirectory
figureDir = os.path.join(figureDir_base, case_prefix)
os.makedirs(figureDir, exist_ok=True)
    

# ================================================================================
# Helper functions for analytical solution
# ================================================================================
def radial_distance_from_axis(x_coord, z_coord, z_s, theta):
    """Perpendicular distance from tilted centerline"""
    return np.abs((z_coord - z_s) * np.cos(theta) - x_coord * np.sin(theta))

def velocity_profile_analytical(r, U_max=1.0, H=1.0):
    """Parabolic profile: u(r) = U_max * (1 - (2r/H)^2)"""
    return U_max * np.maximum(0.0, 1.0 - (2.0 * r / H)**2)

# ================================================================================
# Case discovery and loading
# ================================================================================
# Build case specifications based on mode
if args.flat:
    case_specs = [(f'{case_prefix}-32', 32), (f'{case_prefix}-64', 64), (f'{case_prefix}-128', 128), (f'{case_prefix}-256', 256), (f'{case_prefix}-512', 512)]
else:  # slanted
    case_specs = [(f'{case_prefix}-32', 32), (f'{case_prefix}-64', 64), (f'{case_prefix}-128', 128), (f'{case_prefix}-256', 256), (f'{case_prefix}-512', 512)]

case_paths, case_resolutions, valid_cases = [], [], []
for case_name, resolution in case_specs:
    case_dir = os.path.join(rootDir, case_name)
    if os.path.exists(case_dir):
        plt_dirs = sorted([d for d in os.listdir(case_dir) if d.startswith('plt')])
        if plt_dirs:
            case_paths.append(os.path.join(case_dir, plt_dirs[-1]))
            case_resolutions.append(resolution)
            valid_cases.append(case_name)
            print(f'{case_name}: {plt_dirs[-1]}')

nCases = len(valid_cases)
print(f'\nLoading {nCases} cases')

print(f'\nAnalytical Solution Parameters:')
print(f'  θ = {theta_deg:.2f}°')
print(f'  H = {channelHeight:.6f}')
print(f'  Channel centerline z_s = {z_s_analytical:.6f}')
print(f'  Alignment: {args.align if np.abs(theta_deg) < 0.01 else "N/A (slanted)"}')

# ================================================================================
# Load datasets using yt
# ================================================================================
ds_list, x, y, z, nx, ny, nz, dx, dy, dz, u, v, w, p = [], [], [], [], [], [], [], [], [], [], [], [], [], []
domain_bounds = []  # Store domain boundaries for each case

for case_idx, case_path in enumerate(case_paths):
    print(f'Loading {valid_cases[case_idx]}...')
    try:
        # Load with yt
        ds_yt = yt.load(case_path)
        ds_list.append(ds_yt)
        
        # Get grid information
        dims = ds_yt.domain_dimensions
        domain_left = ds_yt.domain_left_edge.d
        domain_right = ds_yt.domain_right_edge.d
        domain_bounds.append((domain_left, domain_right))
        
        print(f'  Dimensions: {dims}')
        
        # Use covering_grid to get structured data directly
        cube = ds_yt.covering_grid(level=0, left_edge=ds_yt.domain_left_edge, dims=dims)
        
        # Extract coordinate arrays (automatically cell-centered)
        x_arr = cube["index", "x"].d[:, 0, 0]  # Extract 1D by taking first slice in y,z
        y_arr = cube["index", "y"].d[0, :, 0]  # Extract 1D by taking first slice in x,z
        z_arr = cube["index", "z"].d[0, 0, :]  # Extract 1D by taking first slice in x,y
        
        x.append(x_arr)
        y.append(y_arr)
        z.append(z_arr)
        
        nx.append(dims[0])
        ny.append(dims[1])
        nz.append(dims[2])
        
        dx.append(x_arr[1] - x_arr[0] if len(x_arr) > 1 else 0)
        dy.append(y_arr[1] - y_arr[0] if len(y_arr) > 1 else 0)
        dz.append(z_arr[1] - z_arr[0] if len(z_arr) > 1 else 0)
        
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
        
        u.append(u_data_3d)
        v.append(v_data_3d)
        w.append(w_data_3d)
        p.append(p_data_3d)
        
        print(f'  Grid spacing: dx={dx[-1]:.6f}, dy={dy[-1]:.6f}, dz={dz[-1]:.6f}')
        print(f'  Velocity range: u=[{u_data_3d.min():.4f}, {u_data_3d.max():.4f}]\n')
        
    except Exception as e:
        print(f'  Error: {e}\n')
        import traceback
        traceback.print_exc()

nCases = len(ds_list)
print(f'Loaded {nCases} cases\n')

# ================================================================================
# Define Analytical Solution
# ================================================================================
# The axial velocity (along the tilted flow direction) for parabolic pipe flow is:
# u(r) = U_max * (1 - (2r/H)^2)
# where the perpendicular distance from the tilted centerline is:
# r = |x*sin(theta) - (z-z_s)*cos(theta)|
# 
# The centerline is shifted vertically by z_s in the z-direction:
# z_s = H/(2*cos(theta)) + sigma
# 
# and rotated at an angle theta about the point (0, 0, sigma) in x-z plane.

def get_analytical_solution_grid(case_idx):
    X, Y, Z = np.meshgrid(x[case_idx], y[case_idx], z[case_idx], indexing='ij')
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

def axial_velocity(u, w, theta):
    """Compute axial velocity along the tilted flow direction"""
    return u * np.cos(theta) + w * np.sin(theta)

# Compute analytical solutions
ur_exact, r_dist, centerline_shift = [], [], []
for case_idx in range(nCases):
    ur_a, r, z_s = get_analytical_solution_grid(case_idx)
    ur_exact.append(ur_a)
    r_dist.append(r)
    centerline_shift.append(z_s)

print('Analytical solutions computed')

# ================================================================================
# Compute domain bounds for all plots
# ================================================================================
x_min_all = min([bounds[0][0] for bounds in domain_bounds])
x_max_all = max([bounds[1][0] for bounds in domain_bounds])
z_min_all = min([bounds[0][2] for bounds in domain_bounds])
z_max_all = max([bounds[1][2] for bounds in domain_bounds])

# ================================================================================
# Plots of Loaded Data
# ================================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()

for case_idx in range(4):
    ax = axes[case_idx]
    
    if case_idx < nCases:
        # Get middle y-index
        mid_y_idx = ny[case_idx] // 2
        
        # Extract 2D slice at middle y
        x_2d = x[case_idx]
        z_2d = z[case_idx]
        u_2d = u[case_idx][:, mid_y_idx, :]
        u_r = axial_velocity(u_2d, w[case_idx][:, mid_y_idx, :], theta_rad)
        
        # Create 2D coordinate meshes for pcolormesh
        X_mesh, Z_mesh = np.meshgrid(x_2d, z_2d, indexing='ij')
        
        # Plot with "RdYlBu_r" colormap using pcolormesh
        mesh = ax.pcolormesh(X_mesh, Z_mesh, u_r, cmap="RdYlBu_r", shading="auto")
        cbar = plt.colorbar(mesh, ax=ax, shrink=0.5)
        cbar.set_label('$u$ (m/s)', fontsize=11)
        
        # Set axis limits to true domain bounds
        ax.set_xlim(x_min_all, x_max_all)
        ax.set_ylim(z_min_all, z_max_all)
        
        ax.set_xlabel('x (m)', fontsize=12)
        ax.set_ylabel('z', fontsize=12)
        ax.set_title(f'{valid_cases[case_idx]}', fontsize=13, fontweight='bold')
        ax.set_aspect('equal')
    else:
        ax.axis('off')

plt.tight_layout()
plt.savefig(f'{figureDir}/loaded_velocity_visualization.png', dpi=150, bbox_inches='tight')
plt.show()

# ================================================================================
# Plot Geometry Setup
# ================================================================================
fig, ax = plt.subplots(figsize=(16, 10))
fontSize = 22
tol = 0.05

# Domain outline
ax.plot([x_min_all, x_max_all, x_max_all, x_min_all, x_min_all], [z_min_all, z_min_all, z_max_all, z_max_all, z_min_all], 'k-', linewidth=5)

# Channel geometry
z_s = centerline_shift[0]
x_centerline = np.linspace(x_min_all, x_max_all, 100)
z_centerline = z_s + x_centerline * np.tan(theta_rad)
ax.plot(x_centerline, z_centerline, 'b:', linewidth=4, label='Centerline')

# Channel walls and immersed boundary regions
perp_offset = channelHeight / 2.0
x_boundary_full = np.linspace(x_min_all, x_max_all, 100)

# Lower channel boundary
z_lower = z_s + x_boundary_full * np.tan(theta_rad) - perp_offset * np.cos(theta_rad)
# Upper channel boundary
z_upper = z_s + x_boundary_full * np.tan(theta_rad) + perp_offset * np.cos(theta_rad)

# Draw immersed boundary walls in deeppink
ax.plot(x_boundary_full, z_lower, color='deeppink', linewidth=5, label='Immersed boundary')
ax.plot(x_boundary_full, z_upper, color='deeppink', linewidth=5)

# Fill regions outside the channel with grey (gainsboro)
# Lower solid region - from domain bottom to lower boundary
lower_x = np.concatenate([x_boundary_full, [x_max_all, x_min_all]])
lower_z = np.concatenate([z_lower, [z_min_all, z_min_all]])
lower_polygon = plt.matplotlib.patches.Polygon(list(zip(lower_x, lower_z)), 
                                 facecolor='gainsboro', edgecolor='none', alpha=1.0, zorder=0)
ax.add_patch(lower_polygon)

# Upper solid region - from upper boundary to domain top
upper_x = np.concatenate([x_boundary_full, [x_max_all, x_min_all]])
upper_z = np.concatenate([z_upper, [z_max_all, z_max_all]])
upper_polygon = plt.matplotlib.patches.Polygon(list(zip(upper_x, upper_z)), 
                                 facecolor='gainsboro', edgecolor='none', alpha=1.0, zorder=0)
ax.add_patch(upper_polygon)

# Flow direction arrow (paleturquoise)
arrow_x_start = x_min_all + 0.15 * (x_max_all - x_min_all)
arrow_z_start = z_s + arrow_x_start * np.tan(theta_rad)
arrow_x_end = arrow_x_start + 0.25 * (x_max_all - x_min_all)
arrow_z_end = z_s + arrow_x_end * np.tan(theta_rad)

arrow = plt.matplotlib.patches.FancyArrowPatch((arrow_x_start, arrow_z_start), (arrow_x_end, arrow_z_end),
                                 mutation_scale=200, facecolor='paleturquoise', edgecolor='paleturquoise', linewidth=2.5, zorder=2)
ax.add_patch(arrow)

# Add "flow" label (rotated at angle theta, shifted more to the right)
arrow_mid_x = (arrow_x_start + arrow_x_end) / 2
arrow_mid_z = z_s + arrow_mid_x * np.tan(theta_rad)
ax.text(arrow_mid_x, arrow_mid_z, 'flow', fontsize=fontSize - 2, fontweight='bold', rotation=theta_deg, 
        verticalalignment='center', horizontalalignment='center', zorder=3)

# Set axis limits
ax.set_xlim(x_min_all, x_max_all)
ax.set_ylim(z_min_all, z_max_all)

# Axis labels and title
ax.set_xlabel('x (m)', fontsize=fontSize, fontweight='normal')
ax.set_ylabel('z (m)', fontsize=fontSize, fontweight='normal')
geometry_title = 'Flat Channel Geometry' if args.flat else 'Slanted Channel Geometry'
ax.set_title(geometry_title, fontsize=fontSize + 2, fontweight='normal')
ax.tick_params(axis='both', which='major', labelsize=fontSize)

# Legend and grid
ax.legend(fontsize=fontSize - 2, loc='upper left')
ax.grid(True, alpha=0.2, zorder=1)
ax.set_aspect('equal')

# Adjust layout and save
plt.tight_layout()
plt.savefig(f'{figureDir}/geometry_setup.png', dpi=150)
plt.show()
print('Geometry visualization complete\n')

# ================================================================================
# Calculate Errors
# ================================================================================
u_axial_error_max = []
u_axial_error_l2 = []
u_axial_error_max_wall = []
u_axial_error_l2_wall = []

for case_idx in range(nCases):
    # Get midpoint slice in y to avoid periodic boundary effects
    j_mid = ny[case_idx] // 2
    
    # Extract 2D slices at midpoint
    u_2d = u[case_idx][:, j_mid, :]
    w_2d = w[case_idx][:, j_mid, :]
    ur_exact_2d = ur_exact[case_idx][:, j_mid, :]
    r_dist_2d = r_dist[case_idx][:, j_mid, :]
    
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
    
    u_axial_error_max.append(u_axial_max)
    u_axial_error_l2.append(u_axial_l2)
    u_axial_error_max_wall.append(u_axial_max_wall)
    u_axial_error_l2_wall.append(u_axial_l2_wall)
    
    print(f'{valid_cases[case_idx]} (midpoint slice at y-index {j_mid}):')
    print(f'  Overall - Max: {u_axial_max:.6e}, L2: {u_axial_l2:.6e}')
    print(f'  Wall    - Max: {u_axial_max_wall:.6e}, L2: {u_axial_l2_wall:.6e}')

u_axial_error_max = np.array(u_axial_error_max)
u_axial_error_l2 = np.array(u_axial_error_l2)
u_axial_error_max_wall = np.array(u_axial_error_max_wall)
u_axial_error_l2_wall = np.array(u_axial_error_l2_wall)
nx_array = np.array(nx)
nz_array = np.array(nz)
dx_array = np.array(dx)
cell_size_array = dx_array  # Isotropic cell size
grid_count = nx_array * nz_array
print()

# ================================================================================
# Error Convergence vs Cell Size
# ================================================================================
fig, ax = plt.subplots(figsize=(10, 7))

ax.loglog(cell_size_array, u_axial_error_max, 'o-', linewidth=2, markersize=10, label='Max Error')
ax.loglog(cell_size_array, u_axial_error_l2, 's-', linewidth=2, markersize=8, label='L2 Error')
if len(cell_size_array) > 1:
    # Trend lines
    coeffs_max = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_max), 1)
    cell_trend = np.logspace(np.log10(cell_size_array.min()), np.log10(cell_size_array.max()), 50)
    
    # Scale reference lines to pass through mean of max error data
    mean_error_max = np.mean(u_axial_error_max)
    mean_cell_size = np.mean(cell_size_array)
    
    # O(h) reference: error ~ C * h
    C_slope1 = 0.85 * mean_error_max / (mean_cell_size ** 1.0)
    ax.loglog(cell_trend, C_slope1 * cell_trend**1.0, ':', alpha=0.6, linewidth=2, color='gray', label='Slope: 1')
    
    # O(h^1/2) reference: error ~ C * h^0.5
    C_slope05 = 0.85 * mean_error_max / (mean_cell_size ** 0.5)
    ax.loglog(cell_trend, C_slope05 * cell_trend**0.5, '-.', alpha=0.6, linewidth=2, color='gray', label='Slope: 0.5')

ax.set_xlabel('Cell Size (h)', fontsize=14)
ax.set_ylabel('Axial Velocity Error', fontsize=14)
ax.set_title('Convergence vs Cell Size (Axial Velocity)', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(f'{figureDir}/error_convergence.png', dpi=150)
plt.show()

# ================================================================================
# Error Field Visualization
# ================================================================================
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
axes = axes.flatten()

# Find maximum error across all cases for consistent colorbar scaling
max_error_all = 0.0
for case_idx in range(nCases):
    j_mid = ny[case_idx] // 2
    
    # Extract 2D slices at midpoint
    u_2d = u[case_idx][:, j_mid, :]
    w_2d = w[case_idx][:, j_mid, :]
    ur_exact_2d = ur_exact[case_idx][:, j_mid, :]
    
    # Calculate axial velocity errors
    u_axial_num = axial_velocity(u_2d, w_2d, theta_rad)
    u_axial_ana = ur_exact_2d
    u_axial_err = u_axial_num - u_axial_ana
    
    max_error_all = max(max_error_all, np.abs(u_axial_err).max())

print(f'Maximum axial velocity error across all cases: {max_error_all:.6e}\n')

# Create a shared ScalarMappable to ensure all colorbars use the same scale
sm = ScalarMappable(cmap='RdYlBu_r', norm=Normalize(vmin=0, vmax=max_error_all))

for case_idx in range(4):
    ax = axes[case_idx]
    
    if case_idx < nCases:
        # Get axial velocity error field for this case
        j_mid = ny[case_idx] // 2
        
        # Extract 2D slices at midpoint
        u_2d = u[case_idx][:, j_mid, :]
        w_2d = w[case_idx][:, j_mid, :]
        ur_exact_2d = ur_exact[case_idx][:, j_mid, :]
        
        # Calculate axial velocity errors
        u_axial_num = axial_velocity(u_2d, w_2d, theta_rad)
        u_axial_ana = ur_exact_2d
        u_axial_err = u_axial_num - u_axial_ana
        
        # Create meshgrid for x and z
        X_mesh, Z_mesh = np.meshgrid(x[case_idx], z[case_idx], indexing='ij')
        
        # Plot contour with error magnitude using shared normalizer
        contour = ax.contourf(X_mesh, Z_mesh, np.abs(u_axial_err), levels=20, cmap='RdYlBu_r', norm=sm.norm)
        cbar = plt.colorbar(sm, ax=ax, shrink=0.4)
        cbar.set_label('|Error|', fontsize=11)
        
        # Set axis limits to true domain bounds
        ax.set_xlim(x_min_all, x_max_all)
        ax.set_ylim(z_min_all, z_max_all)
        
        ax.set_xlabel('x (m)', fontsize=12)
        ax.set_ylabel('z', fontsize=12)
        ax.set_title(f'{valid_cases[case_idx]} (Nx={nx[case_idx]}, Nz={nz[case_idx]})', fontsize=13, fontweight='bold')
        ax.set_aspect('equal')
    else:
        # Leave empty subplots blank
        ax.axis('off')

plt.tight_layout()
plt.savefig(f'{figureDir}/error_field_visualization.png', dpi=150, bbox_inches='tight')
plt.show()
print('Axial velocity error field visualization complete (all colorbars use same scale)\n')

# ================================================================================
# Axial Velocity - Wall Error Analysis
# ================================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Max Error vs Cell Size
ax = axes[0]
ax.loglog(cell_size_array, u_axial_error_max, 'o-', linewidth=2, markersize=10, label='Overall')
ax.loglog(cell_size_array, u_axial_error_max_wall, 's-', linewidth=2, markersize=8, label='Wall region')
if len(cell_size_array) > 1:
    coeffs_max_overall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_max), 1)
    coeffs_max_wall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_max_wall), 1)
    cell_trend = np.logspace(np.log10(cell_size_array.min()), np.log10(cell_size_array.max()), 50)
    ax.loglog(cell_trend, 0.85 * np.exp(coeffs_max_overall[1]) * cell_trend**coeffs_max_overall[0], '--', alpha=0.5, linewidth=1.5, color='C0', label=f'Overall Trend (slope={coeffs_max_overall[0]:.2f})')
    ax.loglog(cell_trend, 0.85 * np.exp(coeffs_max_wall[1]) * cell_trend**coeffs_max_wall[0], '--', alpha=0.5, linewidth=1.5, color='C1', label=f'Wall Trend (slope={coeffs_max_wall[0]:.2f})')
ax.set_xlabel('Cell Size (h)', fontsize=12)
ax.set_ylabel('Max Error', fontsize=12)
ax.set_title('Max Norm Error vs Cell Size', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, which='both')

# L2 Error vs Cell Size
ax = axes[1]
ax.loglog(cell_size_array, u_axial_error_l2, 'o-', linewidth=2, markersize=10, label='Overall')
ax.loglog(cell_size_array, u_axial_error_l2_wall, 's-', linewidth=2, markersize=8, label='Wall region')
if len(cell_size_array) > 1:
    coeffs_l2_overall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_l2), 1)
    coeffs_l2_wall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_l2_wall), 1)
    cell_trend = np.logspace(np.log10(cell_size_array.min()), np.log10(cell_size_array.max()), 50)
    ax.loglog(cell_trend, 0.85 * np.exp(coeffs_l2_overall[1]) * cell_trend**coeffs_l2_overall[0], '--', alpha=0.5, linewidth=1.5, color='C0', label=f'Overall Trend (slope={coeffs_l2_overall[0]:.2f})')
    ax.loglog(cell_trend, 0.85 * np.exp(coeffs_l2_wall[1]) * cell_trend**coeffs_l2_wall[0], '--', alpha=0.5, linewidth=1.5, color='C1', label=f'Wall Trend (slope={coeffs_l2_wall[0]:.2f})')
ax.set_xlabel('Cell Size (h)', fontsize=12)
ax.set_ylabel('L2 Error', fontsize=12)
ax.set_title('L2 Error vs Cell Size', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(f'{figureDir}/wall_error_analysis.png', dpi=150)
plt.show()
print('Wall error analysis complete\n')

# Print trend slopes
print('='*90)
print('CONVERGENCE RATE ANALYSIS (Slopes vs Cell Size)')
print('='*90)
if len(cell_size_array) > 1:
    coeffs_max_overall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_max), 1)
    coeffs_max_wall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_max_wall), 1)
    coeffs_l2_overall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_l2), 1)
    coeffs_l2_wall = np.polyfit(np.log(cell_size_array), np.log(u_axial_error_l2_wall), 1)
    
    print(f'{"Error Type":<25} {"Region":<20} {"Convergence Rate (p)":<20}')
    print('='*90)
    print(f'{"Max Error":<25} {"Overall":<20} {coeffs_max_overall[0]:<20.2f}')
    print(f'{"Max Error":<25} {"Wall":<20} {coeffs_max_wall[0]:<20.2f}')
    print(f'{"L2 Error":<25} {"Overall":<20} {coeffs_l2_overall[0]:<20.2f}')
    print(f'{"L2 Error":<25} {"Wall":<20} {coeffs_l2_wall[0]:<20.2f}')
    print('='*90)

# Error Summary
print('='*90)
print('AXIAL VELOCITY ERROR SUMMARY')
print('='*90)
print(f'{"Case":<20} {"Grid":<15} {"Max Error (Overall)":<20} {"L2 Error (Overall)":<20} {"Max Error (Wall)":<20} {"L2 Error (Wall)":<20}')
print('='*90)
for i in range(nCases):
    print(f'{valid_cases[i]:<20} {nx[i]} x {nz[i]:<8} {u_axial_error_max[i]:<20.6e} {u_axial_error_l2[i]:<20.6e} {u_axial_error_max_wall[i]:<20.6e} {u_axial_error_l2_wall[i]:<20.6e}')
print('='*90)

# ================================================================================
# Velocity Profile Comparison
# ================================================================================
# Plot the axial velocity at three locations: x_lo (start), x_mid (middle), and x_hi (end)

case_finest = nCases - 1
j_mid = ny[case_finest] // 2

fig, axes = plt.subplots(3, 2, figsize=(14, 12))

# Define x-locations to plot: x_hi * (0.25, 0.5, 0.75)
x_hi = x[case_finest][-1]  # Last x coordinate
x_positions = [x_hi * frac for frac in [0.25, 0.5, 0.75]]

# Find indices closest to these positions
i_locations = [np.argmin(np.abs(x[case_finest] - xp)) for xp in x_positions]
x_labels = ['$0.25 x_{hi}$', '$0.5 x_{hi}$', '$0.75 x_{hi}$']

# Extract 2D slices at midpoint
u_2d = u[case_finest][:, j_mid, :]
w_2d = w[case_finest][:, j_mid, :]
ur_exact_2d = ur_exact[case_finest][:, j_mid, :]

# Compute axial velocities
u_axial_num = axial_velocity(u_2d, w_2d, theta_rad)
u_axial_ana = ur_exact_2d

# First pass: find max error for consistent x-axis scaling on right column
max_error_all = 0.0
for i_loc in i_locations:
    velocity_error = u_axial_num[i_loc, :] - u_axial_ana[i_loc, :]
    max_error_all = max(max_error_all, np.abs(velocity_error).max())

for row_idx, (i_loc, x_label) in enumerate(zip(i_locations, x_labels)):
    x_val = x[case_finest][i_loc]
    
    # Left subplot: velocity profile
    ax = axes[row_idx, 0]
    ax.plot(u_axial_ana[i_loc, :], z[case_finest], 'o', color='tab:red', markersize=5, label='Analytical')
    ax.plot(u_axial_num[i_loc, :], z[case_finest], '-', color='tab:blue', linewidth=3, label='Numerical')
    ax.set_xlabel('Axial Velocity (m/s)', fontsize=12)
    ax.set_ylabel('z (m)', fontsize=12)
    ax.set_title(f'Axial Velocity Profile ({x_label}: x={x_val:.3f})', fontsize=12)
    # Reverse legend order so Numerical appears first
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Right subplot: error
    ax = axes[row_idx, 1]
    velocity_error = u_axial_num[i_loc, :] - u_axial_ana[i_loc, :]
    ax.plot(velocity_error, z[case_finest], '-', color='tab:green', linewidth=3, markersize=6)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Error in Axial Velocity (m/s)', fontsize=12)
    ax.set_ylabel('z (m)', fontsize=12)
    ax.set_title(f'Axial Velocity Error ({x_label}: x={x_val:.3f})', fontsize=12)
    ax.set_xlim(-max_error_all, max_error_all)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{figureDir}/velocity_profile_comparison.png', dpi=150)
plt.show()

print('Post-processing complete!')
