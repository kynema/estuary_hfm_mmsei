"""
Data extraction helpers for the Kynema-SGF Rosario simulation line-sampling
output (Sea Spider "SS" and Stablemoor "STBM" LineSamplers).

This module loads the AMReX native particle data written out by Kynema-SGF
into line_sampling##### folders under a post_processing directory, one folder
per output time:
  - <folder>/sampling_info.yaml: simulation time (s) and the list of
    LineSamplers (index/label/type), e.g. {0: STBM, 1: SS}
  - <folder>/particles/: AMReX particle binary data (Header + Level_*/DATA_*),
    read via AmrexParticleFile from ../../FVCOM-dataExtraction/amrex_particle.py

Kynema-SGF defines z=0 at the water surface (zco negative below the surface),
whereas the Rosario ADCP data (see adcp_data_extraction.py) defines z=0 at the
seafloor. To align the two datasets, the extracted "z_m" (height above
seafloor) is computed as zco + waterdepth, where waterdepth is the local
still-water depth taken from the ADCP data for each unit (STBM uses the fixed
surveyed depth used in plot_adcp_data.py; SS uses the mean measured
"waterdepth" scalar from the SS dataset).

The functions here are intended to be called from other scripts (e.g.
plot_sgf_data.py) rather than run directly. Each loader returns a dict shaped
like the one returned by adcp_data_extraction.load_adcp_data():
  - "time": the full (unfiltered) time vector as a pandas datetime64 array
  - "scalars": DataFrame of per-timestep scalar quantities (time, sim_time_s)
  - "profiles": DataFrame of per-timestep, per-height quantities, including
    "z_m" (height above seafloor), "u"/"v"/"w", "speed", and "vof"
  - "waves": always None (no wave data available from the simulation)
  - "title_range": a human-readable string describing the loaded date range
"""

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import adcp_data_extraction

# Path to the FVCOM-dataExtraction tools directory, which holds
# amrex_particle.py (the AmrexParticleFile reader).
FVCOM_TOOLS_DIR = Path(__file__).parent.parent.parent / "FVCOM-dataExtraction"
if FVCOM_TOOLS_DIR.exists():
    sys.path.insert(0, str(FVCOM_TOOLS_DIR))

from amrex_particle import AmrexParticleFile

# set_id -> label is read per-folder from sampling_info.yaml, but the labels
# we know how to handle are fixed to these two units.
_KNOWN_UNITS = {"ss": "SS", "stbm": "STBM"}


def _discover_line_sampling_folders(root_dir):
    """Return a numerically sorted list of line_sampling##### folders under
    root_dir that contain a sampling_info.yaml file.

    Folders are sorted by the numeric suffix in their name (not lexically,
    since that suffix is not zero-padded to a fixed width and would
    otherwise sort e.g. "line_sampling12292" after "line_sampling122229").
    Folders are then checked, in that numeric (chronological) order. If a
    PermissionError is encountered while checking a folder (e.g. a newer
    folder still being written by a running simulation and its permissions
    haven't been set yet), a warning is printed and any remaining, later
    folders are skipped, since those are the most likely to be incomplete.
    """
    pattern = str(Path(root_dir) / "line_sampling*")
    folders = [f for f in glob.glob(pattern) if Path(f).is_dir()]
    folders = sorted(folders, key=lambda f: int(Path(f).name.removeprefix("line_sampling")))
    accessible_folders = []
    for f in folders:
        try:
            if (Path(f) / "sampling_info.yaml").exists():
                accessible_folders.append(f)
        except PermissionError as e:
            print(f"\nWarning: permission denied checking {f}: {e}")
            print("Stopping further folder discovery; remaining folders may still be incomplete.\n")
            break
    folders = accessible_folders
    if not folders:
        raise FileNotFoundError(f"No line_sampling folders found in {root_dir}")
    return [Path(f) for f in folders]


def _load_one_timestep(folder):
    """Load one line_sampling##### folder into a DataFrame of particles,
    tagged with the simulation time and the sampler label (e.g. "SS",
    "STBM") resolved via that folder's sampling_info.yaml.
    """
    with open(folder / "sampling_info.yaml", "r") as fh:
        info = yaml.safe_load(fh)

    sim_time_s = float(info["time"])
    set_id_to_label = {s["index"]: s["label"] for s in info["samplers"]}

    pfile = AmrexParticleFile(str(folder / "particles"))
    df = pfile()

    df = df.copy()
    df["label"] = df["set_id"].map(set_id_to_label)
    df["sim_time_s"] = sim_time_s
    return df


def _load_all_timesteps(root_dir):
    """Load and concatenate every line_sampling##### folder under root_dir
    into a single combined particle DataFrame.

    Folders are processed in sorted (chronological) order. If a
    PermissionError is encountered (e.g. a newer folder is still being
    written by a running simulation and its permissions haven't been set
    yet), a warning is printed and any remaining, later folders are skipped
    rather than raising, since those are the most likely to be incomplete.
    """
    folders = _discover_line_sampling_folders(root_dir)
    frames = []
    for folder in folders:
        try:
            frames.append(_load_one_timestep(folder))
        except PermissionError as e:
            print(f"\nWarning: permission denied reading {folder}: {e}")
            print("Stopping further folder reads; remaining folders may still be incomplete.\n")
            break
    if not frames:
        raise FileNotFoundError(f"No readable line_sampling folders found in {root_dir}")
    return pd.concat(frames, ignore_index=True)


def _ss_waterdepth():
    """Mean measured waterdepth (m) from the SS ADCP dataset."""
    ss_data = adcp_data_extraction.load_ss_data()
    return ss_data["scalars"]["waterdepth"].mean()


def _seafloor_z_shift(unit):
    """Return the waterdepth (m) used to shift SGF zco (z=0 at surface) to
    z_m (height above seafloor), for the given unit ("SS" or "STBM")."""
    if unit == "STBM":
        return adcp_data_extraction.STBM_WATERDEPTH
    elif unit == "SS":
        return _ss_waterdepth()
    else:
        raise ValueError(f"Unknown unit '{unit}'; expected 'SS' or 'STBM'.")


def _build_result(combined_df, unit, file_prefix, plot_label,
                   start_date=None, start_time=None, vof_threshold=0.5):
    """Filter combined_df to rows for the given unit and build a result dict
    shaped like adcp_data_extraction.load_adcp_data()'s return value.

    Rows with vof < vof_threshold (i.e. mostly air, above the free surface)
    are excluded from "profiles" since velocity there is not physically
    meaningful for comparison to the ADCP data. Pass vof_threshold=None to
    disable this filtering.
    """
    df = combined_df[combined_df["label"] == unit].copy()
    if df.empty:
        raise ValueError(f"No SGF data found for unit '{unit}'.")

    z_shift = _seafloor_z_shift(unit)
    df["z_m"] = df["zco"] + z_shift
    df = df.rename(columns={
        "velocityx": "u",
        "velocityy": "v",
        "velocityz": "w",
    })
    df["speed"] = np.hypot(df["u"], df["v"])

    if start_date is not None and start_time is not None:
        origin = pd.Timestamp(f"{start_date} {start_time}")
        df["time"] = origin + pd.to_timedelta(df["sim_time_s"], unit="s")
    else:
        df["time"] = df["sim_time_s"]

    if vof_threshold is not None:
        profiles = df[df["vof"] >= vof_threshold].copy()
    else:
        profiles = df.copy()
    profiles = profiles[["time", "sim_time_s", "z_m", "u", "v", "w", "speed", "vof"]]

    scalars = df[["time", "sim_time_s"]].drop_duplicates().reset_index(drop=True)

    time = scalars["time"].to_numpy()

    if start_date is not None and start_time is not None:
        title_range = (
            f"{scalars['time'].min():%Y-%m-%d %H:%M:%S} to "
            f"{scalars['time'].max():%Y-%m-%d %H:%M:%S}"
        )
    else:
        title_range = (
            f"t={scalars['sim_time_s'].min():.1f}s to "
            f"t={scalars['sim_time_s'].max():.1f}s"
        )

    return {
        "unit": unit,
        "file_prefix": file_prefix,
        "plot_label": plot_label,
        "time": time,
        "scalars": scalars,
        "profiles": profiles,
        "waves": None,
        "title_range": title_range,
    }


def load_ss_sgf_data(root_dir, start_date=None, start_time=None,
                      vof_threshold=0.5):
    """Load the SGF SS (Sea Spider) LineSampler data across all output
    times found under root_dir."""
    combined = _load_all_timesteps(root_dir)
    return _build_result(combined, "SS", "SGF_SS", "SGF SS",
                          start_date=start_date, start_time=start_time,
                          vof_threshold=vof_threshold)


def load_stbm_sgf_data(root_dir, start_date=None, start_time=None,
                        vof_threshold=0.5):
    """Load the SGF STBM (Stablemoor) LineSampler data across all output
    times found under root_dir."""
    combined = _load_all_timesteps(root_dir)
    return _build_result(combined, "STBM", "SGF_STBM", "SGF STBM",
                          start_date=start_date, start_time=start_time,
                          vof_threshold=vof_threshold)


def load_sgf_data(sgf_unit, root_dir, start_date=None,
                   start_time=None, vof_threshold=0.5):
    """Dispatch to load_ss_sgf_data or load_stbm_sgf_data based on sgf_unit.

    sgf_unit is case-insensitive and accepts "ss"/"SS" or "stbm"/"STBM".
    root_dir is the post_processing directory containing the
    line_sampling##### folders (provided by the caller, e.g. plot_sgf_data.py,
    since this location is expected to move in the future).
    start_date/start_time (e.g. "2024-10-22"/"22:00:00") give the real-world
    UTC datetime corresponding to simulation time=0, used to convert the
    "time" column from raw simulation seconds to real datetimes so it can be
    compared against ADCP timestamps. If omitted, "time" stays as raw
    simulation seconds.
    """
    unit = sgf_unit.strip().lower()
    if unit not in _KNOWN_UNITS:
        raise ValueError(f"Unknown sgf_unit '{sgf_unit}'; expected 'ss' or 'stbm'.")

    if unit == "ss":
        return load_ss_sgf_data(root_dir=root_dir, start_date=start_date,
                                 start_time=start_time, vof_threshold=vof_threshold)
    else:
        return load_stbm_sgf_data(root_dir=root_dir, start_date=start_date,
                                   start_time=start_time, vof_threshold=vof_threshold)
