#!/usr/bin/env python3
"""
Generate simulation case directories with .inp files for flat and slanted configurations.

Usage:
    python generate_cases.py --flat --nx 64,128,256,512
    python generate_cases.py --slanted --nx 64,128,256,512
"""

import os
import argparse
from pathlib import Path
from base_setup import calculate_channel_config



def generate_inp_content(config, template_path, drag_variant='temp'):
    """
    Generate minimal .inp file that includes template and overrides n_cell.
    
    drag_variant: 'temp' (default), 'og'
    """
    # Make path relative to cases directory
    output = f"FILE = ../{template_path.name}\n"
    n_cell = config['n_cell']
    output += f"amr.n_cell = {n_cell[0]} {n_cell[1]} {n_cell[2]}\n"
    
    # Add drag forcing settings based on variant
    if drag_variant == 'og':
        output += "DragForcing.use_original_drag_limiter = true\n"
        output += "DragForcing.use_temporal_drag_limiter = false"
    elif drag_variant == 'temp':
        output += "DragForcing.use_original_drag_limiter = false\n"
        output += "DragForcing.use_temporal_drag_limiter = true"
    
    return output


def create_case_directory(case_path, config, template_path, mode_name, drag_variant='temp'):
    """
    Create a case directory with minimal .inp file that references template.
    """
    
    # Create directory
    os.makedirs(case_path, exist_ok=True)
    
    # Generate .inp content (just references template and overrides n_cell)
    inp_content = generate_inp_content(config, template_path, drag_variant=drag_variant)
    
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
    
    parser.add_argument('--nx', type=str, default=None, 
                       help='Comma-separated list of nx values (default: 32,64,128,256,512 for both modes)')
    parser.add_argument('--dry_run', action='store_true', help='Print what would be created without creating')
    
    args = parser.parse_args()
    
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
    else:
        # Default to temporal variant
        drag_variant = 'temp'
        drag_suffix = 'drag-temp'
    
    # Determine mode and template
    if args.flat:
        mode_name = 'flat'
        case_prefix = f'flat-{drag_suffix}'
        template_file = 'base-flat-poiseuille.inp'
    else:
        mode_name = 'slanted'
        case_prefix = f'slanted-{drag_suffix}'
        template_file = 'base-slanted-poiseuille.inp'
    
    template_path = cases_dir / template_file
    
    if not template_path.exists():
        print(f"Error: Template file not found: {template_path}")
        return 1
    
    print(f"Using template: {template_path}")
    print(f"Cases directory: {cases_dir}")
    print(f"Mode: {mode_name}")
    print()
    
    # Set domain parameters (matching case_setup.py main())
    x_hi = 1000.0
    y_hi = 50.0
    z_hi = 1000.0
    H = 500.0
    
    # Generate cases for each nx
    for nx in nx_values:
        print(f"Generating {case_prefix}-{nx}...")
        
        config = calculate_channel_config(
            flat_mode=args.flat,
            slanted_mode=args.slanted,
            x_hi=x_hi,
            y_hi=y_hi,
            z_hi=z_hi,
            H=H,
            nx=nx,
            align='cf' if args.flat else 'default',
            nx_align=32,
        )
        
        case_name = f"{case_prefix}-{nx}"
        case_path = cases_dir / case_name
        
        if args.dry_run:
            print(f"  Would create: {case_path}")
            print(f"  n_cell: {config['n_cell']}")
            print(f"  BodyForce: ({config['body_force_x']:.6f}, {config['body_force_y']:.6f}, {config['body_force_z']:.6f})")
            print(f"  Drag variant: {drag_variant}")
        else:
            create_case_directory(str(case_path), config, template_path, mode_name, drag_variant=drag_variant)
            print(f"  n_cell: {config['n_cell']}")
            print(f"  BodyForce: ({config['body_force_x']:.6f}, {config['body_force_y']:.6f}, {config['body_force_z']:.6f})")
            print(f"  Drag variant: {drag_variant}")
        print()
    
    return 0


if __name__ == "__main__":
    exit(main())