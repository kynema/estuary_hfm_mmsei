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
from case_setup import calculate_channel_config



def generate_inp_content(config, template_path):
    """
    Generate minimal .inp file that includes template and overrides n_cell.
    """
    # Make path relative to cases directory
    output = f"FILE = ../{template_path.name}\n"
    n_cell = config['n_cell']
    output += f"amr.n_cell = {n_cell[0]} {n_cell[1]} {n_cell[2]}"
    return output


def create_case_directory(case_path, config, template_path, mode_name):
    """
    Create a case directory with minimal .inp file that references template.
    """
    
    # Create directory
    os.makedirs(case_path, exist_ok=True)
    
    # Generate .inp content (just references template and overrides n_cell)
    inp_content = generate_inp_content(config, template_path)
    
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
  # Flat channel with multiple mesh refinements
  python generate_cases.py --flat --nx 64,128,256,512
  
  # Slanted channel with multiple mesh refinements
  python generate_cases.py --slanted --nx 64,128,256,512
        """
    )
    
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument('--flat', action='store_true', help='Generate flat channel cases')
    mode_group.add_argument('--slanted', action='store_true', help='Generate slanted channel cases')
    
    parser.add_argument('--nx', type=str, default=None, 
                       help='Comma-separated list of nx values (default: 64,128,256,512 for flat, 32,64,128,256 for slanted)')
    parser.add_argument('--dry_run', action='store_true', help='Print what would be created without creating')
    
    args = parser.parse_args()
    
    # Set mode-specific defaults for nx if not provided
    if args.nx is None:
        if args.flat:
            args.nx = '64,128,256,512'
        else:  # slanted
            args.nx = '32,64,128,256'
    
    # Get paths relative to script location
    script_dir = Path(__file__).parent
    cases_dir = script_dir.parent / 'cases'
    
    # Parse nx values
    nx_values = [int(x.strip()) for x in args.nx.split(',')]
    
    # Determine mode and template
    if args.flat:
        mode_name = 'flat'
        case_prefix = 'flat-aligned'
        template_file = 'base-flat-poiseuille.inp'
    else:
        mode_name = 'slanted'
        case_prefix = 'slanted'
        template_file = 'base-slanted-poiseuille.inp'
    
    template_path = cases_dir / template_file
    
    if not template_path.exists():
        print(f"Error: Template file not found: {template_path}")
        return 1
    
    print(f"Using template: {template_path}")
    print(f"Cases directory: {cases_dir}")
    print(f"Mode: {mode_name}")
    print()
    
    # Generate cases for each nx
    for nx in nx_values:
        print(f"Generating {case_prefix}-{nx}...")
        
        config = calculate_channel_config(
            flat_mode=args.flat,
            slanted_mode=args.slanted,
            nx=nx,
            align='cf' if args.flat else 'default',
            nx_align=64,
        )
        
        case_name = f"{case_prefix}-{nx}"
        case_path = cases_dir / case_name
        
        if args.dry_run:
            print(f"  Would create: {case_path}")
            print(f"  n_cell: {config['n_cell']}")
            print(f"  BodyForce: ({config['body_force_x']:.6f}, {config['body_force_y']:.6f}, {config['body_force_z']:.6f})")
        else:
            create_case_directory(str(case_path), config, template_path, mode_name)
            print(f"  n_cell: {config['n_cell']}")
            print(f"  BodyForce: ({config['body_force_x']:.6f}, {config['body_force_y']:.6f}, {config['body_force_z']:.6f})")
        print()
    
    return 0


if __name__ == "__main__":
    exit(main())