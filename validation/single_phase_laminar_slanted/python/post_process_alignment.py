#!/usr/bin/env python3
"""
Alignment Comparison Post-Processing: Cell-Face vs Cell-Center

Compare numerical error convergence between cell-face (cf) and cell-center (cc)
grid alignments for flat channel cases with specified drag forcing variants.
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
    velocity_profile_analytical,
    axial_velocity,
    get_analytical_solution_grid,
    get_case_specific_centerlines
)

warnings.filterwarnings('ignore')

print('Modules loaded')

# ================================================================================
# Parse command line arguments
# ================================================================================
parser = argparse.ArgumentParser(description='Compare cell alignment effects (cf vs cc)')

parser.add_argument('--drag', type=str, default='og', 
                   help='Drag models to compare (comma-separated: og,temp,tf1). Default: og')
parser.add_argument('--no_ib', action='store_true', help='Include flat channel without IB cases for comparison')
parser.add_argument('--rel', action='store_true', help='Plot errors relative to maximum velocity (default: absolute errors)')

args = parser.parse_args()

# Parse drag models
drag_models = [d.strip() for d in args.drag.split(',')]
print(f'Comparing alignments for drag models: {drag_models}\n')

# ================================================================================
# Physical parameters
# ================================================================================
theta_deg = 0.0
theta_rad = np.deg2rad(theta_deg)
channelHeight = 1.0
maxVelocity = 1.0

# Define paths
file_dir = os.path.dirname(os.path.abspath(__file__))
rootDir = os.path.join(os.path.dirname(file_dir), 'cases')
figureDir_base = os.path.join(os.path.dirname(file_dir), 'figures')

# Read base configuration from cell-face reference
base_inp_cf = os.path.join(rootDir, 'base-flat-aligned-cf.inp')
if not os.path.exists(base_inp_cf):
    raise RuntimeError(f"Error: Base input file not found: {base_inp_cf}")

# For no-ib, use the no-ib base file
if args.no_ib:
    base_inp_no_ib = os.path.join(rootDir, 'base-flat-no-ib.inp')
    if not os.path.exists(base_inp_no_ib):
        raise RuntimeError(f"Error: Base input file not found: {base_inp_no_ib}")

try:
    z_s_analytical, H_from_inp, theta_from_inp, x_hi_from_inp, z_hi_from_inp, Umax_from_inp, rho_from_inp, mu_from_inp, x_start_cline, x_end_cline, z_start_cline, z_end_cline, L_channel, dp_dx_from_inp = parse_inp_for_channel_geometry(base_inp_cf, return_dp_dx=True)
    
    print(f'Read channel configuration from base-flat-aligned-cf.inp:')
    print(f'  Channel centerline z_s = {z_s_analytical:.6f}')
    print(f'  Channel height H = {H_from_inp:.6f}')
    print(f'  Physics: Umax={Umax_from_inp:.4f}, rho={rho_from_inp:.4f}, mu={mu_from_inp:.4f}\n')
    
    channelHeight = H_from_inp
    maxVelocity = Umax_from_inp
    
except RuntimeError as e:
    print(e)
    raise

# ================================================================================
# Case discovery and loading for each alignment
# ================================================================================
data_by_alignment = {}
resolutions = [32, 64, 128, 256, 512, 1024, 2048]

alignments_to_load = ['cf', 'cc']
if args.no_ib:
    alignments_to_load.append('no-ib')

for align in alignments_to_load:
    data_by_alignment[align] = {}
    
    for drag_model in drag_models:
        print(f'\n{"="*80}')
        print(f'Loading {align.upper()} alignment, {drag_model} drag model')
        print(f'{"="*80}')
        
        # Case naming differs for no-ib: flat-no-ib-{nx} instead of flat-drag-{variant}-{align}-{nx}
        if align == 'no-ib':
            case_prefix = f'flat-no-ib'
        else:
            case_prefix = f'flat-drag-{drag_model}-{align}'
        case_specs = [(f'{case_prefix}-{nx}', nx) for nx in resolutions]
        
        case_paths, valid_cases, case_resolutions = [], [], []
        for case_name, resolution in case_specs:
            case_dir = os.path.join(rootDir, case_name)
            if os.path.exists(case_dir):
                plt_dirs = sorted([d for d in os.listdir(case_dir) if d.startswith('plt')])
                if plt_dirs:
                    case_paths.append(os.path.join(case_dir, plt_dirs[-1]))
                    valid_cases.append(case_name)
                    case_resolutions.append(resolution)
                    print(f'  {case_name}: {plt_dirs[-1]}')
        
        nCases = len(valid_cases)
        if nCases == 0:
            print(f'  WARNING: No cases found for {align.upper()} alignment, {drag_model} drag model')
            continue
            
        print(f'\nLoading {nCases} cases\n')
        
        # Load datasets and compute errors
        u_axial_error_max = []
        u_axial_error_l2 = []
        nx_array = []
        dx_array = []
        
        # Get base z_s for this alignment
        if align == 'no-ib':
            base_inp_align = os.path.join(rootDir, 'base-flat-no-ib.inp')
        else:
            base_inp_align = os.path.join(rootDir, f'base-flat-aligned-{align}.inp')
        z_s_align, _, _, _, _, _, _, _ = parse_inp_for_channel_geometry(base_inp_align, return_dp_dx=False)
        
        # Get case-specific centerlines (only for cc alignment; cf and no-ib use constant z_s)
        if align == 'cc':
            z_s_cases = get_case_specific_centerlines(rootDir, 'flat', align, valid_cases, z_s_align)
        else:
            z_s_cases = [z_s_align] * len(valid_cases)
        
        for case_idx, case_path in enumerate(case_paths):
            print(f'  Loading {valid_cases[case_idx]}...')
            try:
                ds_yt = yt.load(case_path)
                
                # Get grid information
                dims = ds_yt.domain_dimensions
                cube = ds_yt.covering_grid(level=0, left_edge=ds_yt.domain_left_edge, dims=dims)
                
                # Extract coordinate arrays
                x_arr = cube["index", "x"].d[:, 0, 0]
                y_arr = cube["index", "y"].d[0, :, 0]
                z_arr = cube["index", "z"].d[0, 0, :]
                
                nx = dims[0]
                ny = dims[1]
                nz = dims[2]
                dx = x_arr[1] - x_arr[0] if len(x_arr) > 1 else 0
                
                nx_array.append(nx)
                dx_array.append(dx)
                
                # Extract velocity data
                u_data_3d = cube["velocityx"].d.copy()
                w_data_3d = cube["velocityz"].d.copy()
                
                # Mask data outside channel
                X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing='ij')
                r_mask = radial_distance_from_axis(X, Z, z_s_cases[case_idx], theta_rad)
                mask_outside = r_mask > (1.2 * channelHeight / 2.0)
                u_data_3d[mask_outside] = 0.0
                w_data_3d[mask_outside] = 0.0
                
                # Compute analytical solution
                ur_exact, r_dist, _ = get_analytical_solution_grid(x_arr, y_arr, z_arr, 
                                                                   z_s_cases[case_idx], theta_rad, 
                                                                   maxVelocity, channelHeight, 
                                                                   is_slanted=False)
                
                # Calculate errors at midpoint slice
                j_mid = ny // 2
                u_2d = u_data_3d[:, j_mid, :]
                w_2d = w_data_3d[:, j_mid, :]
                ur_exact_2d = ur_exact[:, j_mid, :]
                r_dist_2d = r_dist[:, j_mid, :]
                
                u_axial_num = axial_velocity(u_2d, w_2d, theta_rad)
                u_axial_ana = ur_exact_2d
                u_axial_err = u_axial_num - u_axial_ana
                
                # Overall channel error
                channel_mask = r_dist_2d <= (channelHeight / 2.0)
                u_axial_max = np.max(np.abs(u_axial_err[channel_mask])) if np.any(channel_mask) else 0.0
                u_axial_l2 = np.sqrt(np.sum(u_axial_err[channel_mask]**2)) / np.sqrt(np.sum(channel_mask)) if np.any(channel_mask) else 0.0
                
                # Normalize to relative errors if requested
                if args.rel:
                    u_axial_max /= maxVelocity
                    u_axial_l2 /= maxVelocity
                
                u_axial_error_max.append(u_axial_max)
                u_axial_error_l2.append(u_axial_l2)
                
                print(f'    nx={nx}: Max={u_axial_max:.6e}, L2={u_axial_l2:.6e}')
                
            except Exception as e:
                print(f'    Error: {e}')
                import traceback
                traceback.print_exc()
        
        # Store results
        if len(u_axial_error_max) > 0:
            data_by_alignment[align][drag_model] = {
                'nx': np.array(nx_array),
                'dx': np.array(dx_array),
                'error_max': np.array(u_axial_error_max),
                'error_l2': np.array(u_axial_error_l2)
            }

print(f'\n\nData loading complete')

# Check if any data was collected
total_data_points = sum(len(data_by_alignment[align]) for align in alignments_to_load)
if total_data_points == 0:
    print('ERROR: No cases found to analyze.')
    if args.no_ib:
        print('Make sure cases have been generated with --no_ib:')
        print('  python generate_cases.py --flat --no_ib --nx 32,64,128,256,512')
    else:
        print('Make sure flat cases have been generated:')
        print('  python generate_cases.py --flat --align cf --nx 32,64,128,256,512')
        print('  python generate_cases.py --flat --align cc --nx 32,64,128,256,512')
    exit(1)

# ================================================================================
# Convergence comparison plot (single axis)
# ================================================================================
fig, ax = plt.subplots(1, 1, figsize=(10, 7))

# Keep color consistent within each drag model across alignments
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
drag_color_map = {drag_model: color_cycle[idx % len(color_cycle)] for idx, drag_model in enumerate(drag_models)}

# Collect all points for reference-line scaling
all_dx_data = []
all_error_data = []

for drag_model in drag_models:
    color = drag_color_map[drag_model]

    # Cell Centers: ':o' (dotted with circle markers)
    if drag_model in data_by_alignment['cc']:
        data_cc = data_by_alignment['cc'][drag_model]
        ax.loglog(
            data_cc['dx'],
            data_cc['error_max'],
            ':o',
            linewidth=2.5,
            markersize=8,
            color=color,
            label=f'Cell Centers ({drag_model})',
            alpha=0.9
        )
        all_dx_data.append(data_cc['dx'])
        all_error_data.append(data_cc['error_max'])
    
    # Cell Faces: '--s' (dashed with square markers)
    if drag_model in data_by_alignment['cf']:
        data_cf = data_by_alignment['cf'][drag_model]
        ax.loglog(
            data_cf['dx'],
            data_cf['error_max'],
            '--s',
            linewidth=2.5,
            markersize=8,
            color=color,
            label=f'Cell Faces ({drag_model})',
            alpha=0.9
        )
        all_dx_data.append(data_cf['dx'])
        all_error_data.append(data_cf['error_max'])

# Plot No-IB cases separately (outside drag_model loop, since no-ib has no drag variants)
if args.no_ib and 'no-ib' in data_by_alignment and len(drag_models) > 0:
    # Get a drag_model's data for no-ib (they're all the same, so use first one)
    first_drag = list(data_by_alignment['no-ib'].keys())[0] if data_by_alignment['no-ib'] else None
    if first_drag:
        data_no_ib = data_by_alignment['no-ib'][first_drag]
        ax.loglog(
            data_no_ib['dx'],
            data_no_ib['error_max'],
            '-^',
            linewidth=2.5,
            markersize=8,
            color='purple',
            label='No IB',
            alpha=0.9
        )
        all_dx_data.append(data_no_ib['dx'])
        all_error_data.append(data_no_ib['error_max'])

# Add O(h) and O(h^0.5) reference lines once
if len(all_dx_data) > 0:
    combined_dx = np.concatenate(all_dx_data)
    combined_error = np.concatenate(all_error_data)

    if len(combined_dx) > 1:
        cell_trend = np.logspace(np.log10(combined_dx.min()), np.log10(combined_dx.max()), 50)
        mean_error = np.mean(combined_error)
        mean_cell_size = np.mean(combined_dx)

        C_slope05 = 0.2 * mean_error / (mean_cell_size ** 0.5)
        ax.loglog(cell_trend, C_slope05 * cell_trend**0.5, ':', linewidth=2, color='black', label='Slope: 0.5')

        C_slope1 = 0.1 * mean_error / (mean_cell_size ** 1.0)
        ax.loglog(cell_trend, C_slope1 * cell_trend**1.0, '-.', linewidth=2, color='black', label='Slope: 1')
        
        # Add O(h^2) reference if no-ib cases are included
        if args.no_ib:
            C_slope2 = 0.01 * mean_error / (mean_cell_size ** 2.0)
            ax.loglog(cell_trend, C_slope2 * cell_trend**2.0, '--', linewidth=2, color='black', label='Slope: 2')

ax.set_xlabel('Cell Size (h)', fontsize=12)
ylab = 'Relative Error' if args.rel else 'Absolute Error'
ax.set_ylabel(ylab, fontsize=12)
ax.set_title('Alignment Comparison Across Drag Variants', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
figureDir = os.path.join(figureDir_base, 'alignment-comparison')
os.makedirs(figureDir, exist_ok=True)
fname = 'alignment_convergence_comparison'
if args.rel:
    fname += '_relative'
fname += '.png'
plt.savefig(f'{figureDir}/{fname}', dpi=150, bbox_inches='tight')
plt.show()

print('\nAlignment comparison plots saved to figures/alignment-comparison/')
print('Post-processing complete!')
