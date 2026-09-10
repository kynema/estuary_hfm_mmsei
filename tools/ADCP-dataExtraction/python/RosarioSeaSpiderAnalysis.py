"""
Python port of a subset of Matlab/RosarioSeaSpiderAnalysis.m

Produces two plots for a user-selectable date range from the
SeaSpider_Sig250_Rosario.mat data:
  1. Average speed profile vs. height above seafloor (z=0 at seafloor),
     matching the first panel of SeaSpider_Sig250_Rosario_profiles.png.
     A fixed surface-reflection zone (z=83 to 94 m above seafloor, matching
     the area([...],[94 94],83) shading in Matlab/RosarioSeaSpiderAnalysis.m)
     is shaded red and excluded from the plotted profile line.
  2. Time series of speed at a fixed height above the seafloor (default
     3 m) and at hub height.

Edit start_date, stop_date, and start_time below to narrow the analysis window.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio

MAT_PATH = Path(__file__).parent.parent / "Matlab" / "data" / "SeaSpider_Sig250_Rosario.mat"
OUTPUT_DIR = Path(__file__).parent/ "figures"

DATENUM_EPOCH_OFFSET = 719529  # days from 0000-01-01 to 1970-01-01

# ---- user-editable date range (None, or unset, = use full record) ----
# start_date = "yyyy-mm-dd"
# stop_date = "yyyy-mm-dd"
# start_time = "hh:mm:ss"  # UTC, applied to both start_date and stop_date

# Use full data set
start_date = None
stop_date = None
start_time = None

# Single day (from figure 20 in Thomson et al. 2025)
#start_date = "2024-11-02"
#stop_date = "2024-11-03"
#start_time = "00:00:00" 

# Candidate 1
#start_date = "2024-10-22"
#stop_date = "2024-10-27"
#start_time = "22:00:00"

# Candidate 1 (1 day)
start_date = "2024-10-22"
stop_date = "2024-10-23"
start_time = "22:00:00"

# Candidate 2
#start_date = "2024-11-25"
#stop_date = "2024-12-02"
#start_time = "15:00:00" 

# ---- turbine params (Matlab/RosarioSeaSpiderAnalysis.m) ----
R = 13.5  # rotor radius (m)
hubdepth = 3.5 + R  # hub depth below surface (m)
waterdepth = 90  # m, surveyed depth (not derivable from STBM pressure)

# height above seafloor to report a time series for, in addition to hub depth
near_bottom_height = 3.0  # m

# Specify heights above seafloor for time series plots. 
# Include "avg" to plot the column-averaged speed 
#time_series_heights = [waterdepth-hubdepth, near_bottom_height]  # heights above seafloor for time series
time_series_heights = [waterdepth-hubdepth, near_bottom_height, "avg"]  # include column average


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


def main():
    d = sio.loadmat(MAT_PATH, simplify_cells=True)

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

    start = f"{start_date} {start_time}" if start_date is not None and start_time is not None else start_date
    stop = f"{stop_date} {start_time}" if stop_date is not None and start_time is not None else stop_date
    scalars = filter_dates(scalars, "time", start, stop)
    profiles = filter_dates(profiles, "time", start, stop)

    if scalars.empty or profiles.empty:
        raise ValueError("No data found in the requested date range.")

    display_start = start_date if start_date is not None else time.min().strftime("%Y-%m-%d %H:%M:%S")
    display_stop = stop_date if stop_date is not None else time.max().strftime("%Y-%m-%d %H:%M:%S")
    title_range = f"{display_start} to {display_stop}"

    profiles["speed"] = np.hypot(profiles["u"], profiles["v"])

    # ---- plot 1: average speed profile vs height above seafloor ----
    avg_profile = profiles.groupby("z_m")["speed"].mean().sort_index()

    # Matlab/RosarioSeaSpiderAnalysis.m shades a fixed surface-reflection zone
    # from z=83 to z=94 m above seafloor (area([0 1.1],[94 94],83) in fh8).
    reflection_zone_start = 83
    reflection_zone_end = 94

    mean_waterdepth = scalars["waterdepth"].mean()
    hub_z = mean_waterdepth - hubdepth  # height of hub above seafloor

    # exclude bins affected by surface reflections from the plotted profile line
    plotted_profile = avg_profile[avg_profile.index < reflection_zone_start]

    fig, ax = plt.subplots(figsize=(3, 8))
    ax.axhspan(hub_z - R, hub_z + R, color="0.85", zorder=0)
    ax.axhspan(
        reflection_zone_start,
        reflection_zone_end,
        color="red",
        alpha=0.2,
        zorder=0,
    )
    ax.plot(plotted_profile.values, plotted_profile.index, "b-", linewidth=3)
    ax.axhline(hub_z, color="k", linestyle="--", linewidth=1)
    ax.axhline(hub_z + R, color="k", linestyle=":", linewidth=1)
    ax.axhline(hub_z - R, color="k", linestyle=":", linewidth=1)
    ax.set_xlabel("avg speed [m/s]")
    ax.set_ylabel("height above seafloor [m]")
    ax.set_xlim(left=0)
    ax.set_title(f"Sea Spider avg speed profile\n{title_range}")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "SeaSpider_Sig250_Rosario_profile.png", dpi=150)

    # ---- plot 2: time series of speed at each height in time_series_heights ----
    available_z = profiles["z_m"].unique()

    fig, ax = plt.subplots(figsize=(10, 4))
    for h in time_series_heights:
        if h == "avg":
            column_avg = profiles.groupby("time")["speed"].mean().sort_index()
            ax.plot(column_avg.index, column_avg.values, label="column average", linewidth=3)
        else:
            nearest_z = available_z[np.argmin(np.abs(available_z - h))]
            series = profiles[profiles["z_m"] == nearest_z].sort_values("time")
            ax.plot(series["time"], series["speed"], label=f"z = {nearest_z:.1f} m above seafloor", linewidth=3)
    ax.set_ylabel("speed [m/s]")
    ax.set_title(f"Sea Spider speed time series\n{title_range}")
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "SeaSpider_Sig250_Rosario_timeseries.png", dpi=150)

    print(f"Wrote SeaSpider_Sig250_Rosario_profile.png and SeaSpider_Sig250_Rosario_timeseries.png to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()
