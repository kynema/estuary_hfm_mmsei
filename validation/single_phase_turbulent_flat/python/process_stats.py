import numpy as np
import pandas as pd
import sys
import argparse
import glob
import os
import concurrent.futures
from pathlib import Path
from data import build_case_dir_name

# Add path to FVCOM-dataExtraction tools
estuary_hfm_mmsei = Path(__file__).parent.parent.parent.parent
fvcom_tools = estuary_hfm_mmsei / "tools" / "FVCOM-dataExtraction"
if fvcom_tools.exists():
    sys.path.insert(0, str(fvcom_tools))

# Import the official Kynema utility class
from amrex_particle import AmrexParticleFile


def find_sampling_folder(Re: int, DNS: bool = False, LES: bool = False,
                         IB: bool = False, drag: str = None) -> Path:
    """
    Find the latest sampling folder for the given Reynolds number and variant.
    
    Parameters
    ----------
    Re : int
        Stress Reynolds number
    DNS : bool
        If True, look for DNS variant
    LES : bool
        If True, look for LES variant
    IB : bool
        If True, look for immersed boundary variant
    drag : str
        Drag model: 'og' or 'tf1' (only used with IB)
    
    Returns
    -------
    Path
        Path to the latest completed sampling folder (before sampling09000)
    """
    case_dir = Path(__file__).parent.parent / "cases"
    case_name = build_case_dir_name(Re, DNS=DNS, LES=LES, IB=IB, drag=drag)
    sampling_pattern = case_dir / case_name / "post_processing" / "sampling*"
    
    folders = sorted(glob.glob(str(sampling_pattern)))
    
    if not folders:
        raise FileNotFoundError(
            f"No sampling folders found in {case_dir / case_name / 'post_processing'}/"
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


def _read_binary_file_task(task):
    """Read all grid chunks from a single AMReX binary data file.

    Each DATA_XXXXX file is independent, so this function is designed to run
    concurrently with other invocations on different files.

    Parameters
    ----------
    task : tuple
        (fname, chunks, nints_file, nreals, num_ints) where
        chunks is a list of (npts, offset) per grid chunk in this file.

    Returns
    -------
    list of (pidxs, int_block, real_block) — one entry per chunk
    """
    fname, chunks, nints_file, nreals, num_ints = task
    results = []
    with open(fname, 'rb') as fh:
        for npts, offset in chunks:
            fh.seek(offset)
            ivals = np.fromfile(fh, dtype=np.int32, count=nints_file * npts)
            rvals = np.fromfile(fh, dtype=float, count=nreals * npts)
            # Reshape and extract vectorised — avoids Python-level loops
            ivals_2d = ivals.reshape(npts, nints_file)
            rvals_2d = rvals.reshape(npts, nreals)
            pidxs = ivals_2d[:, 2].copy()                   # particle IDs
            int_block = ivals_2d[:, 2:2 + num_ints].copy()  # [pidx, int_var...]
            results.append((pidxs, int_block, rvals_2d.copy()))
    return results


def load_particle_data_parallel(pfile, nprocs=1):
    """Load AMReX particle binary data in parallel across data files.

    Drop-in replacement for AmrexParticleFile.load_binary_data() + pfile().
    Every DATA_XXXXX binary file is independent, so reads are dispatched to a
    thread pool (disk I/O releases the GIL, enabling true concurrency).

    Parameters
    ----------
    pfile : AmrexParticleFile
        Particle file object after parse_header() has been called.
    nprocs : int
        Number of parallel worker threads.

    Returns
    -------
    pd.DataFrame
        Same format as pfile().
    """
    num_ints   = pfile.num_ints
    nints_file = num_ints + 2          # binary layout: cpu, level, pidx, int_vars...
    nreals     = pfile.ndim + pfile.num_reals

    # Group grid entries by unique binary file (each (lev, idx) → one file)
    file_chunks: dict = {}
    for lev, ginfo in enumerate(pfile.grid_info):
        for idx, npts, offset in ginfo:
            if npts < 1:
                continue
            fname = str(pfile.bin_file_name(lev, idx))
            file_chunks.setdefault(fname, []).append((npts, offset))

    tasks = [
        (fname, chunks, nints_file, nreals, num_ints)
        for fname, chunks in file_chunks.items()
    ]

    print(f"  Reading {len(tasks)} binary data file(s) with {nprocs} worker(s)...")

    idata = np.empty((pfile.num_particles, num_ints), dtype=int)
    rdata = np.empty((pfile.num_particles, nreals),   dtype=float)

    # ThreadPoolExecutor: disk I/O releases the GIL → true parallel reads
    with concurrent.futures.ThreadPoolExecutor(max_workers=nprocs) as executor:
        for file_result in executor.map(_read_binary_file_task, tasks):
            for pidxs, int_block, real_block in file_result:
                idata[pidxs, :] = int_block
                rdata[pidxs, :] = real_block

    idict = {key: idata[:, i] for i, key in enumerate(pfile.int_var_names)}
    rvar_names = ["xco", "yco", "zco"] + pfile.real_var_names
    rdict = {key: rdata[:, i] for i, key in enumerate(rvar_names)}
    idict.update(rdict)
    return pd.DataFrame(idict, index=idata[:, 0])


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Process Kynema particle sampling statistics"
    )
    parser.add_argument("--Re", type=int, default=180,
                        help="Stress Reynolds number (default: 180)")
    parser.add_argument("--DNS", action="store_true",
                        help="Use DNS variant")
    parser.add_argument("--LES", action="store_true",
                        help="Use LES variant")
    parser.add_argument("--IB", action="store_true",
                        help="Use immersed boundary variant")
    parser.add_argument("--drag", type=str, choices=['og', 'tf1'],
                        help="Drag model for IB cases: 'og' or 'tf1'")
    parser.add_argument("--nprocs", type=int,
                        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)),
                        help="Number of parallel workers for reading binary data "
                             "(default: $SLURM_CPUS_PER_TASK or 1)")
    args = parser.parse_args()

    # Validate options
    if not args.DNS and not args.LES:
        parser.error("Either --DNS or --LES must be specified")
    
    if args.DNS and args.LES:
        parser.error("Cannot specify both --DNS and --LES")
    
    if args.IB and not args.LES:
        parser.error("--IB can only be used with --LES")
    
    if args.IB and not args.drag:
        parser.error("--drag option required when using --IB")
    
    if args.drag and not args.IB:
        parser.error("--drag can only be used with --IB")

    # Find the sampling folder
    try:
        sampling_folder = find_sampling_folder(args.Re, DNS=args.DNS, LES=args.LES,
                                              IB=args.IB, drag=args.drag)
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
        pfile = AmrexParticleFile(str(particles_folder))
        pfile.parse_header()
        num_particles = pfile.num_particles if hasattr(pfile, 'num_particles') else None
        if num_particles is not None:
            print(f"  Found {num_particles} particles")
        df = load_particle_data_parallel(pfile, nprocs=args.nprocs)
        print(f"  Loaded {len(df)} particle records")
        print(f"  Columns: {df.columns.tolist()}")
    except FileNotFoundError as e:
        print(f"Error: Could not find required files in {particles_folder}", file=sys.stderr)
        print(f"  Details: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error loading particle file: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
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
        # Offsets are relative to origin, need to convert to absolute coordinates.
        # NOTE: sampling_info.yaml written by the solver does NOT include an
        # 'origin' key, so this always falls back to the default below. The
        # channel half-width (and therefore sampling.channel_stats.origin,
        # which is always z=-delta) must match the convention in
        # case_setup.py's main() (delta=0.01 m for Re>180, 0.005 m for
        # Re=180) or higher-Re cases silently get plane-matched against the
        # wrong z-origin, corrupting particle_stats.txt.
        delta = 0.01 if args.Re > 180 else 0.005
        sampling_origin = sampler.get('origin', [0.0, 0.0, -delta])
        sampling_origin_z = sampling_origin[2] if isinstance(sampling_origin, (list, tuple)) else -delta
        
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
    offsets_array = np.array(offsets, dtype=float)
    absolute_offsets = offsets_array + float(sampling_origin_z)
    
    print(f"  Offsets (relative): {offsets_array.min():.6f} to {offsets_array.max():.6f}")
    print(f"  Offsets (absolute): {absolute_offsets.min():.6f} to {absolute_offsets.max():.6f}")
    print(f"  Particle zco range: {df['zco'].min():.6f} to {df['zco'].max():.6f}")
    
    # For each particle, find the nearest offset plane — O(N log P) via searchsorted,
    # avoids the O(N×P) broadcast that would require ~2 TiB for large datasets.
    zco_vals = df['zco'].to_numpy()
    sorted_offsets = np.sort(absolute_offsets)
    # searchsorted returns insertion index; nearest neighbour is idx-1 or idx
    ins = np.searchsorted(sorted_offsets, zco_vals)
    ins = np.clip(ins, 1, len(sorted_offsets) - 1)
    left  = sorted_offsets[ins - 1]
    right = sorted_offsets[ins]
    nearest_idx = np.where(np.abs(zco_vals - left) <= np.abs(zco_vals - right),
                           ins - 1, ins)
    df['z'] = sorted_offsets[nearest_idx]
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
    
    # Save to the case directory
    case_dir = Path(__file__).parent.parent / "cases"
    case_name = build_case_dir_name(args.Re, DNS=args.DNS, LES=args.LES,
                                    IB=args.IB, drag=args.drag)
    output_file = case_dir / case_name / "particle_stats.txt"
    
    output_file.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_file, final_data, header=header)
    print(f"  Saved particle statistics to {output_file}")
    print(f"    z range: {z_coords.min():.6f} to {z_coords.max():.6f}")
    print(f"    {len(z_coords)} z-planes sampled")

    print("Successfully calculated statistics directly from Kynema binary sampling!")
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
