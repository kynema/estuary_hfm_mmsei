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
resolutions = [32, 64, 128, 256, 512]

for align in ['cf', 'cc']:
    data_by_alignment[align] = {}
    
    for drag_model in drag_models:
        print(f'\n{"="*80}')
        print(f'Loading {align.upper()} alignment, {drag_model} drag model')
        print(f'{"="*80}')
        
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
        u_axial_error_max_wall = []
        u_axial_error_l2_wall = []
        nx_array = []
        dx_array = []
        
        # Get base z_s for this alignment
        base_inp_align = os.path.join(rootDir, f'base-flat-aligned-{align}.inp')
        z_s_align, _, _, _, _, _, _, _ = parse_inp_for_channel_geometry(base_inp_align, return_dp_dx=False)
        
        # Get case-specific centerlines for cc alignment
        z_s_cases = get_case_specific_centerlines(rootDir, 'flat', align, valid_cases, z_s_align)
        
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
                
                # Near-wall error (quarter of channel height from wall)
                wall_distance = channelHeight / 4.0
                near_wall_mask = (r_dist_2d > (channelHeight / 2.0 - wall_distance)) & (r_dist_2d <= (channelHeight / 2.0))
                u_axial_max_wall = np.max(np.abs(u_axial_err[near_wall_mask])) if np.any(near_wall_mask) else 0.0
                u_axial_l2_wall = np.sqrt(np.sum(u_axial_err[near_wall_mask]**2)) / np.sqrt(np.sum(near_wall_mask)) if np.any(near_wall_mask) else 0.0
                
                u_axial_error_max.append(u_axial_max)
                u_axial_error_l2.append(u_axial_l2)
                u_axial_error_max_wall.append(u_axial_max_wall)
                u_axial_error_l2_wall.append(u_axial_l2_wall)
                
                print(f'    nx={nx}: Overall Max={u_axial_max:.6e}, L2={u_axial_l2:.6e}  |  Wall Max={u_axial_max_wall:.6e}, L2={u_axial_l2_wall:.6e}')
                
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
                'error_l2': np.array(u_axial_error_l2),
                'error_max_wall': np.array(u_axial_error_max_wall),
                'error_l2_wall': np.array(u_axial_error_l2_wall)
            }

print(f'\n\nData loading complete')

# ================================================================================
# Convergence comparison plots
# ================================================================================
fig, axes = plt.subplots(1, len(drag_models), figsize=(7*len(drag_models), 6))
if len(drag_models) == 1:
    axes = [axes]

for col_idx, drag_model in enumerate(drag_models):
    ax = axes[col_idx]
    
    # Collect all data for this drag model to compute reference lines
    all_dx_data = []
    all_error_data = []
    
    for align in ['cf', 'cc']:
        if drag_model in data_by_alignment[align]:
            data = data_by_alignment[align][drag_model]
            align_label = 'Cell Centers' if align == 'cc' else 'Cell Faces'
            ax.loglog(data['dx'], data['error_max'], 'o-', linewidth=2.5, 
                     markersize=10, label=f'{align_label} ({drag_model})', alpha=0.85)
            all_dx_data.append(data['dx'])
            all_error_data.append(data['error_max'])
    
    # Add reference lines if we have data
    if len(all_dx_data) > 0:
        # Combine all data for reference line scaling
        combined_dx = np.concatenate(all_dx_data)
        combined_error = np.concatenate(all_error_data)
        
        if len(combined_dx) > 1:
            cell_trend = np.logspace(np.log10(combined_dx.min()), np.log10(combined_dx.max()), 50)
            
            # Scale reference lines to pass through mean of error data
            mean_error = np.mean(combined_error)
            mean_cell_size = np.mean(combined_dx)
            
            # O(h) reference: error ~ C * h
            C_slope1 = 0.85 * mean_error / (mean_cell_size ** 1.0)
            ax.loglog(cell_trend, C_slope1 * cell_trend**1.0, ':', alpha=0.6, linewidth=2, color='gray', label='Slope: 1')
            
            # O(h^0.5) reference: error ~ C * h^0.5
            C_slope05 = 1.25 * mean_error / (mean_cell_size ** 0.5)
            ax.loglog(cell_trend, C_slope05 * cell_trend**0.5, '-.', alpha=0.6, linewidth=2, color='gray', label='Slope: 0.5')
    
    ax.set_xlabel('Cell Size (h)', fontsize=12)
    ax.set_ylabel('Maximum Error', fontsize=12)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
figureDir = os.path.join(figureDir_base, 'alignment-comparison')
os.makedirs(figureDir, exist_ok=True)
plt.savefig(f'{figureDir}/alignment_convergence_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

print('\nAlignment comparison plots saved to figures/alignment-comparison/')
print('Post-processing complete!')
