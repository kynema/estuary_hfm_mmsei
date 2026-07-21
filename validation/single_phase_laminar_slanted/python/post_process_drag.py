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
from shared_functions import (
    parse_inp_for_channel_geometry,
    radial_distance_from_axis,
    axial_velocity,
    discover_case_variants_with_align,
    get_analytical_solution_grid,
    print_effective_height_correction,
    get_case_specific_centerlines
)

warnings.filterwarnings('ignore')

print('Modules loaded')

# ================================================================================
# Helper function to read time from Header file
# ================================================================================
def get_time_from_header(plt_dir):
    """Read simulation time from plt###/Header file.
    
    Header format varies: time is at line (num_variables + 4).
    - Line 2: num_variables
    - Line (num_variables + 4): time value
    Rounds to nearest integer to handle floating point noise.
    """
    header_file = os.path.join(plt_dir, 'Header')
    try:
        with open(header_file, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 2:
                # Read number of variables from line 2 (0-indexed line 1)
                try:
                    num_variables = int(lines[1].strip())
                except:
                    return None
                
                # Time is at line (num_variables + 4) (0-indexed: num_variables + 3)
                time_line_idx = num_variables + 3
                if len(lines) > time_line_idx:
                    time_line = lines[time_line_idx].strip()
                    # Parse time value - could be "time = 100.0" or just "100.008612413216781"
                    if '=' in time_line:
                        time_value = float(time_line.split('=')[1].strip())
                    else:
                        time_value = float(time_line)
                    # Round to nearest integer to handle floating point noise
                    return round(time_value)
    except:
        pass
    return None

# ================================================================================
# Parse command line arguments
# ================================================================================
parser = argparse.ArgumentParser(description='Post-process slanted/flat channel simulations')

mode_group = parser.add_mutually_exclusive_group(required=True)
mode_group.add_argument('--flat', action='store_true', help='Post-process flat channel cases')
mode_group.add_argument('--slanted', action='store_true', help='Post-process slanted channel cases')

parser.add_argument('--align', type=str, default='cf', choices=['cf', 'cc'], help='Grid alignment for flat cases: cf (cell face, default) or cc (cell center)')
parser.add_argument('--rel', action='store_true', help='Plot errors relative to maximum velocity (default: absolute errors)')
parser.add_argument('--drag', type=str, default='og,temp,tf1', help='Drag models to analyze, comma-separated (default: og,temp,tf1)')
parser.add_argument('--time', type=str, default=None, help='Simulation time(s) to analyze, comma-separated (e.g., 100,200,300). Max 3 times. Default: use last plot file')
args = parser.parse_args()

# Parse --time argument: can be single time or comma-separated list (max 4)
times_to_analyze = []
if args.time is not None:
    times_to_analyze = [int(float(t.strip())) for t in args.time.split(',')]
    if len(times_to_analyze) > 3:
        raise ValueError("Maximum 3 times can be specified")
    times_to_analyze = sorted(times_to_analyze)  # Sort for consistent line styling
    print(f'Analyzing {len(times_to_analyze)} time(s): {times_to_analyze}')
else:
    times_to_analyze = [None]  # None means use latest plot file

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
file_dir = os.path.dirname(os.path.abspath(__file__))
rootDir = os.path.join(os.path.dirname(file_dir), 'cases')
figureDir_base = os.path.join(os.path.dirname(file_dir), 'figures')
base_inp = os.path.join(rootDir, 'base-flat-aligned-cf.inp' if args.flat else 'base-slanted.inp')

# Figure directory for drag comparison (mode_name-drag-comparison)
if args.flat:
    figureDir = os.path.join(figureDir_base, f'{mode_name}-drag-comparison-{args.align}')
else:
    figureDir = os.path.join(figureDir_base, f'{mode_name}-drag-comparison')
os.makedirs(figureDir, exist_ok=True)

# Parse --drag argument
drag_variants = [v.strip() for v in args.drag.split(',')]
valid_variants = {'og', 'temp', 'tf1'}
for variant in drag_variants:
    if variant not in valid_variants:
        raise ValueError(f"Invalid drag variant: {variant}. Must be one of: og, temp, tf1")

# Verify at least one case exists for the selected variants
print('Verifying availability of selected drag variants...')
available_variants = set()
resolutions = [32, 64, 128, 256, 512, 1024, 2048]

for variant in drag_variants:
    for resolution in resolutions:
        if mode_name == 'flat':
            case_name = f'{mode_name}-drag-{variant}-{args.align}-{resolution}'
        else:
            case_name = f'{mode_name}-drag-{variant}-{resolution}'
        case_dir = os.path.join(rootDir, case_name)
        if os.path.exists(case_dir):
            plt_dirs = [d for d in os.listdir(case_dir) if d.startswith('plt') and 'old' not in d]
            if plt_dirs:
                available_variants.add(variant)
                break

# Warn if requested variants not found
missing = set(drag_variants) - available_variants
if missing:
    print(f'Warning: No cases found for drag variant(s): {", ".join(missing)}')

# Keep only variants that were found
drag_variants = [v for v in drag_variants if v in available_variants]
if not drag_variants:
    raise RuntimeError(f"No cases found for any of the specified drag variants: {args.drag}")

print(f'Using drag variants: {drag_variants}\n')

# Initialize parameters
theta_rad = np.deg2rad(theta_deg)
channelHeight = 500.0
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
# Case discovery and loading
# ================================================================================
# ================================================================================
# Build case specifications to load selected variants and times
# ================================================================================
resolutions = [32, 64, 128, 256, 512, 1024, 2048]

# First, discover all available cases (independent of time)
available_cases = {variant: [] for variant in drag_variants}
for variant in drag_variants:
    for resolution in resolutions:
        if mode_name == 'flat':
            case_name = f'{mode_name}-drag-{variant}-{args.align}-{resolution}'
        else:
            case_name = f'{mode_name}-drag-{variant}-{resolution}'
        case_dir = os.path.join(rootDir, case_name)
        if os.path.exists(case_dir):
            plt_dirs = sorted([d for d in os.listdir(case_dir) if d.startswith('plt') and 'old' not in d])
            if plt_dirs:
                available_cases[variant].append((case_name, case_dir, plt_dirs))

# Nested dictionaries: case_paths_dict[variant][time_idx][case_idx]
case_paths_dict = {variant: {time_idx: [] for time_idx in range(len(times_to_analyze))} for variant in drag_variants}
valid_cases_dict = {variant: {time_idx: [] for time_idx in range(len(times_to_analyze))} for variant in drag_variants}

# Line styles for different times
line_styles = ['-', '-.', ':']
case_times = [(times_to_analyze[i], line_styles[i]) for i in range(len(times_to_analyze))]

# Now, for each available case, find the best plot file at each requested time
for variant in drag_variants:
    for case_name, case_dir, plt_dirs in available_cases[variant]:
        print(f'{case_name}:')
        # For each requested time, find the closest plot file
        for time_idx, (target_time, linestyle) in enumerate(case_times):
            if target_time is None:
                # Use last plot file
                case_paths_dict[variant][time_idx].append(os.path.join(case_dir, plt_dirs[-1]))
                valid_cases_dict[variant][time_idx].append(case_name)
                print(f'  {plt_dirs[-1]} (latest)')
            else:
                # Find plot file exactly matching requested time (within 1s tolerance)
                best_plt = None
                best_time_diff = float('inf')
                best_time_val = None
                for plt_dir in plt_dirs:
                    full_plt_path = os.path.join(case_dir, plt_dir)
                    file_time = get_time_from_header(full_plt_path)
                    if file_time is not None:
                        time_diff = abs(file_time - target_time)
                        if time_diff < best_time_diff:
                            best_time_diff = time_diff
                            best_plt = full_plt_path
                            best_time_val = file_time
                # Only use if time matches within 10 second tolerance
                if best_plt is not None and best_time_diff <= 10:
                    case_paths_dict[variant][time_idx].append(best_plt)
                    valid_cases_dict[variant][time_idx].append(case_name)
                    print(f'  t={best_time_val}s (target={target_time}s)')

# Report found cases
for variant in drag_variants:
    nCases = len(valid_cases_dict[variant][0])
    print(f'\nFound {nCases} {variant.upper()} cases')

print()
if len(times_to_analyze) > 1 or times_to_analyze[0] is not None:
    print(f'Analyzing times: {[t for t, _ in case_times if t is not None]}')

print(f'\nAnalytical Solution Parameters:')
print(f'  θ = {theta_deg:.2f}°')
print(f'  H = {channelHeight:.6f}')
print(f'  Channel centerline z_s = {z_s_analytical:.6f}')
print(f'  Alignment: {args.align if np.abs(theta_deg) < 0.01 else "N/A (slanted)"}')

# ================================================================================
# Read Resolution-Specific Centerline Positions (for --align cc mode)
# ================================================================================
# Note: Centerlines are the same for all times, so use first time's case list
z_s_case_specific_dict = {}
for variant in drag_variants:
    # Extract case names from first time index (they're the same for all times)
    valid_cases_flat = valid_cases_dict[variant][0]
    z_s_case_specific_dict[variant] = get_case_specific_centerlines(
        rootDir, mode_name, args.align, valid_cases_flat, z_s_analytical
    )
if mode_name == 'flat' and args.align == 'cc':
    print()

# ================================================================================
# Load datasets using yt for both drag variants and times
# ================================================================================
# Initialize data structures: data_by_variant[variant][time_idx][resolution_idx]
data_by_variant = {}

for variant in drag_variants:
    data_by_variant[variant] = {}
    for time_idx in range(len(times_to_analyze)):
        data_by_variant[variant][time_idx] = {
            'ds_list': [], 'x': [], 'y': [], 'z': [],
            'nx': [], 'ny': [], 'nz': [],
            'dx': [], 'dy': [], 'dz': [],
            'u': [], 'v': [], 'w': [], 'p': [],
            'domain_bounds': [],
            'ur_exact': [], 'r_dist': [], 'centerline_shift': []
        }

# Load data for each variant, time, and case
for variant in drag_variants:
    for time_idx, (target_time, linestyle) in enumerate(case_times):
        print(f'\nLoading {variant.upper()} variant at time_idx={time_idx}...')
        for case_idx, case_path in enumerate(case_paths_dict[variant][time_idx]):
            case_name = valid_cases_dict[variant][time_idx][case_idx]
            print(f'Loading {case_name}...')
            try:
                # Load with yt
                ds_yt = yt.load(case_path)
                data_by_variant[variant][time_idx]['ds_list'].append(ds_yt)
                
                # Get grid information
                dims = ds_yt.domain_dimensions
                domain_left = ds_yt.domain_left_edge.d
                domain_right = ds_yt.domain_right_edge.d
                data_by_variant[variant][time_idx]['domain_bounds'].append((domain_left, domain_right))
                
                print(f'  Dimensions: {dims}')
                
                # Use covering_grid to get structured data directly
                cube = ds_yt.covering_grid(level=0, left_edge=ds_yt.domain_left_edge, dims=dims)
                
                # Extract coordinate arrays (automatically cell-centered)
                x_arr = cube["index", "x"].d[:, 0, 0]
                y_arr = cube["index", "y"].d[0, :, 0]
                z_arr = cube["index", "z"].d[0, 0, :]
                
                data_by_variant[variant][time_idx]['x'].append(x_arr)
                data_by_variant[variant][time_idx]['y'].append(y_arr)
                data_by_variant[variant][time_idx]['z'].append(z_arr)
                
                data_by_variant[variant][time_idx]['nx'].append(dims[0])
                data_by_variant[variant][time_idx]['ny'].append(dims[1])
                data_by_variant[variant][time_idx]['nz'].append(dims[2])
                
                data_by_variant[variant][time_idx]['dx'].append(x_arr[1] - x_arr[0] if len(x_arr) > 1 else 0)
                data_by_variant[variant][time_idx]['dy'].append(y_arr[1] - y_arr[0] if len(y_arr) > 1 else 0)
                data_by_variant[variant][time_idx]['dz'].append(z_arr[1] - z_arr[0] if len(z_arr) > 1 else 0)
                
                # Extract 3D field data
                u_data_3d = cube["velocityx"].d.copy()
                v_data_3d = cube["velocityy"].d.copy()
                w_data_3d = cube["velocityz"].d.copy()
                p_data_3d = cube["p"].d.copy()
                
                # Mask numerical data outside analysis region immediately after reading
                X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
                r_mask = radial_distance_from_axis(X, Z, z_s_analytical, theta_rad)
                
                if args.slanted:
                    mask_outside = r_mask > (1.2 * channelHeight / 2.0)
                else:
                    mask_outside = r_mask > (1.2 * channelHeight / 2.0)
                
                u_data_3d[mask_outside] = 0.0
                v_data_3d[mask_outside] = 0.0
                w_data_3d[mask_outside] = 0.0
                p_data_3d[mask_outside] = 0.0
                
                data_by_variant[variant][time_idx]['u'].append(u_data_3d)
                data_by_variant[variant][time_idx]['v'].append(v_data_3d)
                data_by_variant[variant][time_idx]['w'].append(w_data_3d)
                data_by_variant[variant][time_idx]['p'].append(p_data_3d)
                
                print(f'  Grid spacing: dx={data_by_variant[variant][time_idx]["dx"][-1]:.6f}')
                print(f'  Velocity range: u=[{u_data_3d.min():.4f}, {u_data_3d.max():.4f}]')
                
            except Exception as e:
                print(f'  Error: {e}')
                import traceback
                traceback.print_exc()

# Compute analytical solutions for each variant and time
for variant in drag_variants:
    for time_idx in range(len(times_to_analyze)):
        for case_idx in range(len(data_by_variant[variant][time_idx]['x'])):
            ur_a, r, z_s = get_analytical_solution_grid(
                data_by_variant[variant][time_idx]['x'][case_idx],
                data_by_variant[variant][time_idx]['y'][case_idx],
                data_by_variant[variant][time_idx]['z'][case_idx],
                z_s_case_specific_dict[variant][case_idx], theta_rad, maxVelocity, channelHeight, 
                is_slanted=args.slanted
            )
            data_by_variant[variant][time_idx]['ur_exact'].append(ur_a)
            data_by_variant[variant][time_idx]['r_dist'].append(r)
            data_by_variant[variant][time_idx]['centerline_shift'].append(z_s)

print('\nAnalytical solutions computed for all variants and times')

# ================================================================================
# Calculate Errors for Both Variants and All Times
# ================================================================================
errors_by_variant = {}

for variant in drag_variants:
    errors_by_variant[variant] = {}
    for time_idx, (target_time, linestyle) in enumerate(case_times):
        print(f'\nComputing errors for {variant.upper()} variant (time_idx={time_idx}):')
        errors_by_variant[variant][time_idx] = {
            'error_max': [],
            'error_l2': [],
            'nx_array': [],
            'cell_size_array': []
        }
        
        nCases_variant = len(data_by_variant[variant][time_idx]['x'])
        for case_idx in range(nCases_variant):
            j_mid = data_by_variant[variant][time_idx]['ny'][case_idx] // 2
            
            # Extract 2D slices at midpoint
            u_2d = data_by_variant[variant][time_idx]['u'][case_idx][:, j_mid, :]
            w_2d = data_by_variant[variant][time_idx]['w'][case_idx][:, j_mid, :]
            ur_exact_2d = data_by_variant[variant][time_idx]['ur_exact'][case_idx][:, j_mid, :]
            r_dist_2d = data_by_variant[variant][time_idx]['r_dist'][case_idx][:, j_mid, :]
            
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
            
            # Normalize to relative errors if requested
            if args.rel:
                u_axial_max /= maxVelocity
                u_axial_l2 /= maxVelocity
            
            errors_by_variant[variant][time_idx]['error_max'].append(u_axial_max)
            errors_by_variant[variant][time_idx]['error_l2'].append(u_axial_l2)
            errors_by_variant[variant][time_idx]['nx_array'].append(data_by_variant[variant][time_idx]['nx'][case_idx])
            errors_by_variant[variant][time_idx]['cell_size_array'].append(data_by_variant[variant][time_idx]['dx'][case_idx])

# Convert to numpy arrays for plotting
for variant in drag_variants:
    for time_idx in range(len(times_to_analyze)):
        errors_by_variant[variant][time_idx]['error_max'] = np.array(errors_by_variant[variant][time_idx]['error_max'])
        errors_by_variant[variant][time_idx]['error_l2'] = np.array(errors_by_variant[variant][time_idx]['error_l2'])
        errors_by_variant[variant][time_idx]['cell_size_array'] = np.array(errors_by_variant[variant][time_idx]['cell_size_array'])
        errors_by_variant[variant][time_idx]['nx_array'] = np.array(errors_by_variant[variant][time_idx]['nx_array'])

# ================================================================================
# Error Convergence vs Cell Size (Drag Comparison with Multiple Times)
# ================================================================================
# Define styling for all possible variants
all_colors = {'og': 'C0', 'temp': 'C1', 'tf1': 'C2'}
all_markers = {'og': 'o', 'temp': 's', 'tf1': '^'}
all_labels = {'og': 'OG Limiter', 'temp': 'Temporal Limiter', 'tf1': 'Temporal Implementation'}

# Select only the colors/markers/labels for variants being compared
colors = {v: all_colors[v] for v in drag_variants}
markers = {v: all_markers[v] for v in drag_variants}
labels = {v: all_labels[v] for v in drag_variants}

fig, ax = plt.subplots(figsize=(10, 7))

# Collect all cell sizes and errors for reference line scaling
all_cell_sizes = []
all_errors = []

for variant in drag_variants:
    for time_idx, (target_time, linestyle) in enumerate(case_times):
        cell_size_array = errors_by_variant[variant][time_idx]['cell_size_array']
        u_axial_error_max = errors_by_variant[variant][time_idx]['error_max']
        
        # Build legend label
        if target_time is not None:
            legend_label = f'{labels[variant]} (t={target_time}s)'
        else:
            legend_label = f'{labels[variant]}'
        
        ax.loglog(cell_size_array, u_axial_error_max, marker=markers[variant], linestyle=linestyle, 
                  linewidth=2.5, markersize=10, label=legend_label, color=colors[variant])
        
        all_cell_sizes.extend(cell_size_array)
        all_errors.extend(u_axial_error_max)

# Add reference slopes using combined data
all_cell_sizes = np.array(all_cell_sizes)
all_errors = np.array(all_errors)

if len(all_cell_sizes) > 1:
    cell_trend = np.logspace(np.log10(all_cell_sizes.min()), np.log10(all_cell_sizes.max()), 50)
    
    # Scale reference lines to pass through mean of max error data
    mean_error_max = np.mean(all_errors)
    mean_cell_size = np.mean(all_cell_sizes)
    
    # O(h^1/2) reference: error ~ C * h^0.5
    C_slope05 = 1.1 * mean_error_max / (mean_cell_size ** 0.5)
    ax.loglog(cell_trend, C_slope05 * cell_trend**0.5, ':', linewidth=2, color='black', label='Slope: 0.5')

    # O(h) reference: error ~ C * h
    C_slope1 = 0.6 * mean_error_max / (mean_cell_size ** 1.0)
    ax.loglog(cell_trend, C_slope1 * cell_trend**1.0, '-.', linewidth=2, color='black', label='Slope: 1')

ax.set_xlabel('Cell Size (h)', fontsize=14)
ylab = 'Relative Error' if args.rel else 'Absolute Error'
ax.set_ylabel(ylab, fontsize=14)
ax.set_xticks([1e0, 4e0, 1e1, 4e1])
ax.set_xticklabels(['$10^0$', '$4 \\times 10^0$', '$10^1$', '$4 \\times 10^1$'])
#ax.set_yticks([4e0, 1e1, 4e1])
#ax.set_yticklabels(['$4 \\times 10^0$', '$10^1$', '$4 \\times 10^1$'])
ax.set_title(f'Drag Forcing Comparison: {mode_name.title()} Channel Convergence', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
fname = 'error_convergence'
if args.rel:    
    fname += '_relative'
fname += '.png'
plt.savefig(f'{figureDir}/{fname}', dpi=150)
plt.show()

print('\nPost-processing complete!')
