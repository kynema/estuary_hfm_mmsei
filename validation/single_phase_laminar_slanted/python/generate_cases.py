#!/usr/bin/env python3
"""
Generate simulation case directories with .inp files for flat and slanted configurations.

Usage:
    python generate_cases.py --flat --nx 64,128,256,512
    python generate_cases.py --slanted --nx 64,128,256,512
"""

import os
import argparse
import math
from pathlib import Path
from base_setup import calculate_channel_config



def calculate_channel_center_cc(nz, z_hi=1000.0, H=500.0):
    """Calculate channel centerline z-position aligned to cell centers."""
    dz = z_hi / nz
    
    # Cell center coordinates in z
    z_centers = [(i + 0.5) * dz for i in range(nz)]
    
    # Find cell centers closest to target boundaries
    z_lo_target = (z_hi - H)/2
    z_hi_target = z_lo_target + H
    z_lo = min(z_centers, key=lambda z: abs(z - z_lo_target))
    z_hi = min(z_centers, key=lambda z: abs(z - z_hi_target))
    z_center = 0.5 * (z_hi + z_lo)
    
    return z_center


def generate_inp_content(config, template_path, drag_variant='temp', mode_name='flat', nx=32, align='cf', ny=4, nz=32, x_hi=1000.0, y_hi=50.0, z_hi=1000.0, H=500.0, no_ib=False):
    """
    Generate minimal .inp file that includes template and overrides n_cell.
    
    drag_variant: 'temp' (default), 'og', 'tf1'
    mode_name: 'flat' or 'slanted'
    align: 'cf' (cell face) or 'cc' (cell center) for flat cases
    nx: resolution (for computing mac_proj smoothing parameters)
    no_ib: True for flat cases without immersed boundary
    """
    import math
    
    # Make path relative to cases directory
    output = f"FILE = ../{template_path.name}\n"
    n_cell = config['n_cell']
    output += f"amr.n_cell = {n_cell[0]} {n_cell[1]} {n_cell[2]}\n"
    
    # Add drag forcing settings based on variant (skip for no_ib cases)
    if not no_ib:
        if drag_variant == 'og':
            output += "DragForcing.terrain_use_original_limiter = true\n"
        elif drag_variant == 'temp':
            output += "DragForcing.terrain_use_original_limiter = false\n"
            output += "DragForcing.terrain_use_temporal_limiter = true\n"
        elif drag_variant == 'tf1':
            output += "DragForcing.terrain_use_original_limiter = false\n"
            output += "DragForcing.bc_forcing_time_factor = 1.\n"
    
    # For flat cases with cell-center alignment, override ChannelBuilder parameters
    # Skip ChannelBuilder parameters for no-ib cases
    if mode_name == 'flat' and align == 'cc' and not no_ib:
        z_center = calculate_channel_center_cc(nz, z_hi)
        output += f"ChannelBuilder.s1.segment_start_point = 0.0 62.5 {z_center:10.10f}\n"
        output += f"ChannelBuilder.s1.segment_end_point = {x_hi:10.1f} 62.5 {z_center:10.10f}"
    
    # Add mac_proj smoothing parameters for slanted cases with nx >= 128
    if mode_name == 'slanted' and nx >= 128:
        num_smooth = 4 * (1 + math.log2(nx / 128.0))
        num_smooth = int(round(num_smooth))
        output += f"mac_proj.num_pre_smooth = {num_smooth}\n"
        output += f"mac_proj.num_post_smooth = {num_smooth}"
    
    return output


def create_case_directory(case_path, inp_content, mode_name):
    """
    Create a case directory with minimal .inp file.
    """
    
    # Create directory
    os.makedirs(case_path, exist_ok=True)
    
    # Write .inp file
    inp_path = os.path.join(case_path, f"{mode_name}.inp")
    with open(inp_path, 'w') as f:
        f.write(inp_content)
    
    print(f"Created: {inp_path}")
    return inp_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate case directories with configured .inp files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Flat channel with temporal drag limiter (default)
  python generate_cases.py --flat --nx 32,64,128,256,512
  
  # Slanted channel with original drag limiter
  python generate_cases.py --slanted --og --nx 32,64,128,256,512
  
  # Flat channel with original drag limiter
  python generate_cases.py --flat --og --nx 32,64,128,256,512
        """
    )
    
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--flat', action='store_true', help='Generate flat channel cases')
    mode_group.add_argument('--slanted', action='store_true', help='Generate slanted channel cases')
    
    drag_group = parser.add_mutually_exclusive_group()
    drag_group.add_argument('--og', action='store_true', help='Use original drag limiter variant')
    drag_group.add_argument('--temp', action='store_true', help='Use temporal drag limiter variant (default)')
    drag_group.add_argument('--tf1', action='store_true', help='Use temporal implementation with time factor 1')
    
    parser.add_argument('--align', type=str, default='cf', choices=['cf', 'cc'], 
                       help='Channel alignment for flat cases: cf (cell face, default) or cc (cell center)')
    parser.add_argument('--nx', type=str, default=None, 
                       help='Comma-separated list of nx values (default: 32,64,128,256,512 for both modes)')
    parser.add_argument('--dry_run', action='store_true', help='Print what would be created without creating')
    parser.add_argument('--no_ib', action='store_true', help='Generate flat channel without IB')
    
    args = parser.parse_args()
    
    # Validate --no_ib can only be used with --flat
    if args.no_ib and not args.flat:
        print("Error: --no_ib can only be used with --flat")
        return 1
    
    # Set mode-specific defaults for nx if not provided
    if args.nx is None:
        if args.flat:
            args.nx = '32,64,128,256,512'
        else:  # slanted
            args.nx = '32,64,128,256,512'
    
    # Get paths relative to script location
    script_dir = Path(__file__).parent
    cases_dir = script_dir.parent / 'cases'
    
    # Parse nx values
    nx_values = [int(x.strip()) for x in args.nx.split(',')]
    
    # Determine drag variant
    if args.og:
        drag_variant = 'og'
        drag_suffix = 'drag-og'
    elif args.temp:
        drag_variant = 'temp'
        drag_suffix = 'drag-temp'
    elif args.tf1:
        drag_variant = 'tf1'
        drag_suffix = 'drag-tf1'
    else:
        # Default to temporal variant
        drag_variant = 'temp'
        drag_suffix = 'drag-temp'
    
    # Determine mode and template
    if args.flat:
        mode_name = 'flat'
        align_suffix = f"-{args.align}"
        case_prefix = f'flat-{drag_suffix}{align_suffix}'
        # Use alignment-specific template
        if args.align == 'cc':
            template_file = 'base-flat-aligned-cc.inp'
        else:
            template_file = 'base-flat-aligned-cf.inp'
        if args.no_ib:
            case_prefix = "flat-no-ib"
            template_file = 'base-flat-no-ib.inp'
    else:
        mode_name = 'slanted'
        case_prefix = f'slanted-{drag_suffix}'
        template_file = 'base-slanted.inp'
        args.align = None  # Alignment only applies to flat cases
    
    template_path = cases_dir / template_file
    
    if not template_path.exists():
        print(f"Error: Template file not found: {template_path}")
        return 1
    
    print(f"Using template: {template_path}")
    print(f"Cases directory: {cases_dir}")
    print(f"Mode: {mode_name}")
    if args.no_ib:
        print(f"Immersed boundary: DISABLED")
    print()
    
    # Set domain parameters (matching case_setup.py main())
    x_hi = 1000.0
    y_hi = 50.0
    z_hi = 1000.0
    H = 500.0
    
    # Generate cases for each nx
    for nx in nx_values:
        print(f"Generating {case_prefix}-{nx}...")
        
        # For flat cases, ny is always 4; for slanted, use same
        ny = 4
        nz = nx
        
        config = calculate_channel_config(
            flat_mode=args.flat,
            slanted_mode=args.slanted,
            x_hi=x_hi,
            y_hi=y_hi,
            z_hi=z_hi,
            H=H,
            nx=nx,
            align=args.align if args.flat else 'default',
            nx_align=32,
        )
        
        # For no-ib cases, override n_cell: nx, ny, nx/2
        if args.no_ib:
            nz_no_ib = nx // 2
            config['n_cell'] = (nx, ny, nz_no_ib)
        
        case_name = f"{case_prefix}-{nx}"
        case_path = cases_dir / case_name
        
        # Generate input content
        inp_content = generate_inp_content(config, template_path, drag_variant=drag_variant, mode_name=mode_name, nx=nx, align=args.align if args.flat else 'cf', ny=ny, nz=nz, x_hi=x_hi, y_hi=y_hi, z_hi=z_hi, H=H, no_ib=args.no_ib)
        
        if args.dry_run:
            print(f"  Would create: {case_path}")
            print(f"  Input file content:")
            print("  " + "\n  ".join(inp_content.split("\n")))
        else:
            create_case_directory(str(case_path), inp_content, mode_name)
        print(f"  n_cell: {config['n_cell']}")
        print(f"  BodyForce: ({config['body_force_x']:.6f}, {config['body_force_y']:.6f}, {config['body_force_z']:.6f})")
        print(f"  Drag variant: {drag_variant}")
        if args.flat and args.align == 'cc':
            z_center = calculate_channel_center_cc(nz, z_hi)
            print(f"  Channel center (cc): z={z_center:10.6f}")
        if mode_name == 'slanted' and nx >= 128:
            num_smooth = 4 * (1 + math.log2(nx / 128.0))
            num_smooth = int(round(num_smooth))
            print(f"  MAC projection pre smoothers: {num_smooth}")
            print(f"  MAC projection post smoothers: {num_smooth}")
        print()
    
    return 0


if __name__ == "__main__":
    exit(main())