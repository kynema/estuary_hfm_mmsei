"""
Data extraction helpers for the Rosario ADCP datasets.

This module loads and filters the two Rosario ADCP datasets produced by the
Matlab scripts in ../Matlab:
  - Sea Spider (SS): a seafloor tripod with an uplooking Sig250, stored in
    Matlab/data/SeaSpider_Sig250_Rosario.mat
  - Stablemoor (STBM): a mid-water mooring with up/down-looking Sig500s,
    stored in Matlab/data/STBM_Sig500s_Rosario.mat

The functions here are intended to be called from other scripts (e.g.
plot_adcp_data.py) rather than run directly. Each loader returns:
  - "time": the full (unfiltered) time vector as a pandas datetime64 array
  - "scalars": DataFrame of per-timestep scalar quantities (time, hubspeed,
    hubTI, and waterdepth for SS)
  - "profiles": DataFrame of per-timestep, per-height quantities, always
    including a "z_m" (height above seafloor) and "speed" column
  - "waves": DataFrame of surface wave statistics (STBM only, may be None)
  - "title_range": a human-readable string describing the requested date range
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio

MATLAB_DATA_DIR = Path(__file__).parent.parent / "Matlab" / "data"
SS_MAT_PATH = MATLAB_DATA_DIR / "SeaSpider_Sig250_Rosario.mat"
STBM_MAT_PATH = MATLAB_DATA_DIR / "STBM_Sig500s_Rosario.mat"

DATENUM_EPOCH_OFFSET = 719529  # days from 0000-01-01 to 1970-01-01


def matlab_datenum_to_datetime(datenum):
    """Convert MATLAB datenum (days since 0000-01-01) to pandas datetime64."""
    return pd.to_datetime(np.asarray(datenum) - DATENUM_EPOCH_OFFSET, unit="D")


def filter_dates(df, time_col, start, stop):
    mask = pd.Series(True, index=df.index)
    if start is not None:
        mask &= df[time_col] >= pd.Timestamp(start)
    if stop is not None:
        mask &= df[time_col] <= pd.Timestamp(stop)
    return df[mask]


def _date_bounds(start_date, stop_date, start_time):
    start = f"{start_date} {start_time}" if start_date is not None and start_time is not None else start_date
    stop = f"{stop_date} {start_time}" if stop_date is not None and start_time is not None else stop_date
    return start, stop


def _title_range(time, start_date, stop_date):
    display_start = start_date if start_date is not None else time.min().strftime("%Y-%m-%d %H:%M:%S")
    display_stop = stop_date if stop_date is not None else time.max().strftime("%Y-%m-%d %H:%M:%S")
    return f"{display_start} to {display_stop}"


def load_ss_data(start_date=None, stop_date=None, start_time=None, mat_path=SS_MAT_PATH):
    """Load and filter the Sea Spider (SS) dataset.

    Returns a dict with keys "time", "scalars", "profiles", "waves" (always
    None for SS), and "title_range". "profiles" includes u, v, w, and a
    derived "speed" column (horizontal speed, hypot(u, v)).
    """
    d = sio.loadmat(mat_path, simplify_cells=True)

    time = matlab_datenum_to_datetime(d["t"])
    z = np.asarray(d["z"])  # fixed 1-D depth grid, meters above seafloor
    n_time, n_z = d["u"].shape

    scalars = pd.DataFrame(
        {
            "time": time,
            "hubspeed": d["hubspeed"],
            "hubTI": d["hubTI"],
            "waterdepth": d["waterdepth"],
        }
    )

    profiles = pd.DataFrame(
        {
            "time": np.repeat(time, n_z),
            "z_m": np.tile(z, n_time),
            "u": d["u"].ravel(),
            "v": d["v"].ravel(),
            "w": d["w"].ravel(),
        }
    )
    profiles["speed"] = np.hypot(profiles["u"], profiles["v"])

    start, stop = _date_bounds(start_date, stop_date, start_time)
    scalars = filter_dates(scalars, "time", start, stop)
    profiles = filter_dates(profiles, "time", start, stop)

    if scalars.empty or profiles.empty:
        raise ValueError("No SS data found in the requested date range.")

    return {
        "unit": "SS",
        "file_prefix": "SeaSpider_Sig250_Rosario",
        "plot_label": "Sea Spider",
        "time": time,
        "scalars": scalars,
        "profiles": profiles,
        "waves": None,
        "title_range": _title_range(time, start_date, stop_date),
    }


def load_stbm_data(start_date=None, stop_date=None, start_time=None, mat_path=STBM_MAT_PATH):
    """Load and filter the Stablemoor (STBM) dataset.

    Returns a dict with keys "time", "scalars", "profiles", "waves" (None if
    the 'Waves' struct is not present in the .mat file), and "title_range".
    "profiles" uses the fixed-height-above-seabed gridded data
    (fixedz_profiles in the original Matlab/RosarioStablemoorAnalysis.m) and
    includes "speed", "power_density", and "principal_axis" columns.
    """
    d = sio.loadmat(mat_path, simplify_cells=True)

    time = matlab_datenum_to_datetime(d["t"])
    fixedz = np.asarray(d["fixedz"], dtype=float)  # fixed grid, meters above seafloor
    n_time, n_fixedz = d["speed_fixedz"].shape

    scalars = pd.DataFrame(
        {
            "time": time,
            "hubspeed": d["hubspeed"],
            "hubTI": d["hubTI"],
        }
    )

    profiles = pd.DataFrame(
        {
            "time": np.repeat(time, n_fixedz),
            "z_m": np.tile(fixedz, n_time),
            "speed": d["speed_fixedz"].ravel(),
            "power_density": d["power_fixedz"].ravel(),
            "principal_axis": d["axis_fixedz"].ravel(),
        }
    )

    start, stop = _date_bounds(start_date, stop_date, start_time)
    scalars = filter_dates(scalars, "time", start, stop)
    profiles = filter_dates(profiles, "time", start, stop)

    if scalars.empty or profiles.empty:
        raise ValueError("No STBM data found in the requested date range.")

    waves = None
    if "Waves" in d:
        waves = pd.DataFrame(d["Waves"])
        waves["time"] = matlab_datenum_to_datetime(waves["time"])
        waves = filter_dates(waves, "time", start, stop)
        if waves.empty:
            waves = None

    return {
        "unit": "STBM",
        "file_prefix": "STBM_Sig500s_Rosario",
        "plot_label": "STBM",
        "time": time,
        "scalars": scalars,
        "profiles": profiles,
        "waves": waves,
        "title_range": _title_range(time, start_date, stop_date),
    }


def load_adcp_data(adcp_unit, start_date=None, stop_date=None, start_time=None):
    """Dispatch to load_ss_data or load_stbm_data based on adcp_unit.

    adcp_unit is case-insensitive and accepts "ss"/"SS" or "stbm"/"STBM".
    """
    unit = adcp_unit.strip().lower()
    if unit == "ss":
        return load_ss_data(start_date=start_date, stop_date=stop_date, start_time=start_time)
    elif unit == "stbm":
        return load_stbm_data(start_date=start_date, stop_date=stop_date, start_time=start_time)
    else:
        raise ValueError(f"Unknown adcp_unit '{adcp_unit}'; expected 'ss' or 'stbm'.")


def hub_height_above_seafloor(data, hubdepth, waterdepth):
    """Compute the hub height above the seafloor for a loaded dataset.

    For SS, uses the mean of the measured "waterdepth" scalar; for STBM (no
    measured waterdepth), falls back to the given fixed waterdepth.
    """
    if data["unit"] == "SS":
        return data["scalars"]["waterdepth"].mean() - hubdepth
    return waterdepth - hubdepth


def reflection_zone():
    """Return (start, end) heights above seafloor (m) of the fixed
    surface-reflection zone, matching the area([...],[94 94],83) shading in
    the original Matlab scripts. This zone should be excluded from plotted
    average speed profiles.
    """
    return 83, 94
