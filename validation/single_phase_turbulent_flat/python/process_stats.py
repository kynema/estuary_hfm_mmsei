import numpy as np
import pandas as pd
import sys
import argparse
import glob
from pathlib import Path

# Add path to kynema-sgf tools
kynema_sgf_tools = Path("/Users/dmontgo2/Documents/Kynema/kynema-sgf/tools")
if kynema_sgf_tools.exists():
    sys.path.insert(0, str(kynema_sgf_tools))

# Import the official Kynema utility class
from amrex_particle import AmrexParticleFile


def find_sampling_folder(Re: int, IB: bool = False) -> Path:
    """
    Find the latest sampling folder for the given Reynolds number and IB variant.
    
    Parameters
    ----------
    Re : int
        Stress Reynolds number
    IB : bool
        If True, look for ReTau{Re}_IB/ directory
    
    Returns
    -------
    Path
        Path to the latest completed sampling folder (before sampling09000)
    """
    case_dir = Path(__file__).parent.parent / "cases"
    
    if IB:
        sampling_pattern = case_dir / f"ReTau{Re}_IB" / "post_processing" / "sampling*"
    else:
        sampling_pattern = case_dir / f"ReTau{Re}" / "post_processing" / "sampling*"
    
    folders = sorted(glob.glob(str(sampling_pattern)))
    
    if not folders:
        case_variant = f"ReTau{Re}_IB" if IB else f"ReTau{Re}"
        raise FileNotFoundError(
            f"No sampling folders found in {case_dir / case_variant / 'post_processing'}/"
        )
    
    # Filter to only actual directories
    folders = [f for f in folders if Path(f).is_dir()]
    
    if not folders:
        raise FileNotFoundError("No sampling folders found")
    
    # Return the latest folder that has sampling_info.yaml
    # (in case the latest is still being written to by the running simulation)
    for folder_path in reversed(folders):
        info_file = Path(folder_path) / "sampling_info.yaml"
        if info_file.exists():
            return Path(folder_path)
    
    # If no folder has the info file yet, warn and return latest anyway
    print(f"Warning: No completed sampling folder found (latest is being written to)",
          file=sys.stderr)
    return Path(folders[-1])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Process Kynema particle sampling statistics"
    )
    parser.add_argument("--Re", type=int, default=180,
                        help="Stress Reynolds number (default: 180)")
    parser.add_argument("--IB", action="store_true",
                        help="Use immersed boundary variant")
    args = parser.parse_args()

    # Find the sampling folder
    try:
        sampling_folder = find_sampling_folder(args.Re, IB=args.IB)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    
    print(f"Processing sampling data from: {sampling_folder}")
    
    # Load particle data using AmrexParticleFile
    # The particles are stored in sampling_folder/particles/
    particles_folder = sampling_folder / "particles"
    
    if not particles_folder.exists():
        print(f"Error: particles folder not found: {particles_folder}", file=sys.stderr)
        return 1
    
    try:
        # Initialize AmrexParticleFile with the particles directory
        pfile = AmrexParticleFile(str(particles_folder))
        
        # Try parsing the Header file directly
        pfile.parse_header()
        
        # Look for particle count
        num_particles = pfile.num_particles if hasattr(pfile, 'num_particles') else None
        if num_particles is not None:
            print(f"  Found {num_particles} particles")
        
        # Load all binary data
        pfile.load_binary_data()
    except FileNotFoundError as e:
        print(f"Error: Could not find required files in {particles_folder}", file=sys.stderr)
        print(f"  Details: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error loading particle file: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    
    # Extract the dataframe
    try:
        df = pfile()
        print(f"  Loaded {len(df)} particle records")
        print(f"  Columns: {df.columns.tolist()}")
    except Exception as e:
        print(f"Error: Could not extract particle data: {e}", file=sys.stderr)
        return 1
    
    # Load sampling metadata to map particles to z-planes
    import yaml
    sampling_info_file = sampling_folder / "sampling_info.yaml"
    if not sampling_info_file.exists():
        print(f"Error: sampling_info.yaml not found at {sampling_info_file}", file=sys.stderr)
        return 1
    
    try:
        with open(sampling_info_file, 'r') as f:
            sampling_info = yaml.safe_load(f)
        
        # Extract offsets (z positions) from samplers
        samplers = sampling_info.get('samplers', [])
        if not samplers:
            print(f"Error: No samplers found in {sampling_info_file}", file=sys.stderr)
            return 1
        
        # For now, assume single sampler (channel_stats)
        sampler = samplers[0]
        offsets = sampler.get('offsets', [])
        
        # Get the sampling origin for coordinate mapping
        # Offsets are relative to origin, need to convert to absolute coordinates
        sampling_origin = sampler.get('origin', [0.0, 0.0, -0.005])
        sampling_origin_z = sampling_origin[2] if isinstance(sampling_origin, (list, tuple)) else -0.005
        
    except yaml.YAMLError as e:
        print(f"Error parsing {sampling_info_file}: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error processing sampling metadata: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    
    # Map particles to z-planes based on their z-coordinates
    # The zco column contains the actual z-coordinate in domain coords
    # The offsets from .inp are relative to sampling.origin, so convert to absolute
    offsets_array = np.array(offsets)
    absolute_offsets = offsets_array + sampling_origin_z
    
    print(f"  Offsets (relative): {offsets_array.min():.6f} to {offsets_array.max():.6f}")
    print(f"  Offsets (absolute): {absolute_offsets.min():.6f} to {absolute_offsets.max():.6f}")
    print(f"  Particle zco range: {df['zco'].min():.6f} to {df['zco'].max():.6f}")
    
    # For each particle, find the nearest offset plane (in absolute coordinates)
    def find_nearest_plane(z_coord):
        idx = np.argmin(np.abs(absolute_offsets - z_coord))
        return absolute_offsets[idx]
    
    df['z'] = df['zco'].apply(find_nearest_plane)
    matched_planes = df['z'].nunique()
    print(f"  Matched to {matched_planes} unique z-planes (expected {len(offsets_array)})")
    
    # Rename velocity columns to match kynema-sgf conventions
    # velocityx = u (flow), velocityy = v_sgf (span), velocityz = w_sgf (wall-normal)
    df = df.rename(columns={
        'velocityx': 'u',
        'velocityy': 'v_sgf',
        'velocityz': 'w_sgf'
    })
    
    # Cleanly sort the dataset by Z height
    df = df.sort_values(by="z")

    # Group by each unique Z-coordinate layer to compute profiles
    profiles = df.groupby("z").agg(
        u_mean=('u', 'mean'),
        u_var=('u', 'var'),
        v_sgf_mean=('v_sgf', 'mean'),
        v_sgf_var=('v_sgf', 'var'),
        w_sgf_mean=('w_sgf', 'mean'),
        w_sgf_var=('w_sgf', 'var')
    ).reset_index()

    # Extract profiles as clean NumPy arrays
    z_coords = profiles['z'].to_numpy()
    u_mean = profiles['u_mean'].to_numpy()
    u_variance = profiles['u_var'].to_numpy()
    v_sgf_mean = profiles['v_sgf_mean'].to_numpy()
    v_sgf_variance = profiles['v_sgf_var'].to_numpy()
    w_sgf_mean = profiles['w_sgf_mean'].to_numpy()
    w_sgf_variance = profiles['w_sgf_var'].to_numpy()

    # Save data into a compact ASCII columns file
    final_data = np.column_stack((
        z_coords, 
        u_mean, u_variance,
        v_sgf_mean, v_sgf_variance,
        w_sgf_mean, w_sgf_variance
    ))
    header = "z_coord   u_mean   u_var   v_sgf_mean   v_sgf_var   w_sgf_mean   w_sgf_var"
    
    # Save to the case directory (ReTau{Re}/ or ReTau{Re}_IB/)
    case_dir = Path(__file__).parent.parent / "cases"
    if args.IB:
        output_file = case_dir / f"ReTau{args.Re}_IB" / "particle_stats.txt"
    else:
        output_file = case_dir / f"ReTau{args.Re}" / "particle_stats.txt"
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_file, final_data, header=header)
    print(f"  Saved particle statistics to {output_file}")
    print(f"    z range: {z_coords.min():.6f} to {z_coords.max():.6f}")
    print(f"    {len(z_coords)} z-planes sampled")

    print("Successfully calculated statistics directly from Kynema binary sampling!")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
