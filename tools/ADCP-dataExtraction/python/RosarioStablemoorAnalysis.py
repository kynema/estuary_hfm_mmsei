"""
Python port of a subset of Matlab/RosarioStablemoorAnalysis.m

Produces plots for a user-selectable date range from the 
STBM_Sig500s_Rosario.mat data:
  1. Average speed profile vs. height above seafloor (z=0 at seafloor),
     using the fixed-height-above-seabed gridded data (fixedz_profiles),
     matching the first panel of STBM_Sig500s_Rosario_profiles.png. A fixed
     surface-reflection zone (z=83 to 94 m above seafloor, matching the
     area([...],[94 94],83) shading in Matlab/RosarioStablemoorAnalysis.m)
     is shaded red and excluded from the plotted profile line.
  2. Time series of speed at a fixed height above the seafloor (default
     3 m) and at hub height.
  3. Surface wave time series (significant wave height, peak wave period),
     matching the fh21 plot in Matlab/RosarioStablemoorAnalysis.m. Uses the
     'Waves' struct array in STBM_Sig500s_Rosario.mat, which is populated
     from the externally preprocessed STBM_Sig500up_preprocessed.mat file
     (see the "surface waves" section of the Matlab script). If 'Waves' is
     not present in the .mat file, this plot is skipped with a warning.
  4. Wave height vs. peak period scatter, matching the fh22 plot.

Edit start_date, stop_date, and start_time below to narrow the analysis window.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio

MAT_PATH = Path(__file__).parent.parent / "Matlab" / "data" / "STBM_Sig500s_Rosario.mat"
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

# ---- turbine / site params (Matlab/RosarioStablemoorAnalysis.m) ----
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
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    d = sio.loadmat(MAT_PATH, simplify_cells=True)

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

    fixedz_profiles = pd.DataFrame(
        {
            "time": np.repeat(time, n_fixedz),
            "z_m": np.tile(fixedz, n_time),
            "speed": d["speed_fixedz"].ravel(),
            "power_density": d["power_fixedz"].ravel(),
            "principal_axis": d["axis_fixedz"].ravel(),
        }
    )

    start = f"{start_date} {start_time}" if start_date is not None and start_time is not None else start_date
    stop = f"{stop_date} {start_time}" if stop_date is not None and start_time is not None else stop_date
    scalars = filter_dates(scalars, "time", start, stop)
    fixedz_profiles = filter_dates(fixedz_profiles, "time", start, stop)

    if scalars.empty or fixedz_profiles.empty:
        raise ValueError("No data found in the requested date range.")

    display_start = start_date if start_date is not None else time.min().strftime("%Y-%m-%d %H:%M:%S")
    display_stop = stop_date if stop_date is not None else time.max().strftime("%Y-%m-%d %H:%M:%S")
    title_range = f"{display_start} to {display_stop}"

    # ---- plot 1: average speed profile vs height above seafloor ----
    avg_profile = fixedz_profiles.groupby("z_m")["speed"].mean().sort_index()

    # Matlab/RosarioStablemoorAnalysis.m shades a fixed surface-reflection zone
    # from z=83 to z=94 m above seafloor (area([0 1.5],[94 94],83) in fh8).
    reflection_zone_start = 83
    reflection_zone_end = 94

    hub_z = waterdepth - hubdepth  # height of hub above seafloor

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
    ax.set_ylim(0, 100)
    ax.set_title(f"STBM avg speed profile\n{title_range}")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "STBM_Sig500s_Rosario_profile.png", dpi=150)

    # ---- plot 2: time series of speed at each height in time_series_heights ----
    available_z = fixedz_profiles["z_m"].unique()

    fig, ax = plt.subplots(figsize=(10, 4))
    for h in time_series_heights:
        if h == "avg":
            column_avg = fixedz_profiles.groupby("time")["speed"].mean().sort_index()
            ax.plot(column_avg.index, column_avg.values, label="column average", linewidth=3)
        else:
            nearest_z = available_z[np.argmin(np.abs(available_z - h))]
            series = fixedz_profiles[fixedz_profiles["z_m"] == nearest_z].sort_values("time")
            ax.plot(series["time"], series["speed"], label=f"z = {nearest_z:.1f} m above seafloor", linewidth=3)
    ax.set_ylabel("speed [m/s]")
    ax.set_title(f"STBM speed time series\n{title_range}")
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "STBM_Sig500s_Rosario_timeseries.png", dpi=150)

    written = ["STBM_Sig500s_Rosario_profile.png", "STBM_Sig500s_Rosario_timeseries.png"]

    # ---- plots 3 & 4: surface waves (only if 'Waves' is present) ----
    if "Waves" in d:
        waves = pd.DataFrame(d["Waves"])
        waves["time"] = matlab_datenum_to_datetime(waves["time"])
        waves = filter_dates(waves, "time", start, stop)

        if waves.empty:
            print("No wave data found in the requested date range; skipping wave plots.")
        else:
            # plot 3: wave height / peak period time series
            fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
            axes[0].plot(waves["time"], waves["sigwaveheight"], "+", markersize=10, markeredgewidth=2.5)
            axes[0].set_ylim(0, 2)
            axes[0].set_ylabel("$H_s$ [m]")
            axes[1].plot(waves["time"], waves["peakwaveperiod"], "+", markersize=10, markeredgewidth=2.5)
            axes[1].set_ylim(0, 10)
            axes[1].set_ylabel("$T_p$ [s]")
            fig.autofmt_xdate()
            fig.suptitle(f"STBM surface waves\n{title_range}")
            fig.tight_layout()
            fig.savefig(OUTPUT_DIR / "STBM_Sig500s_Rosario_waves_timeseries.png", dpi=150)
            written.append("STBM_Sig500s_Rosario_waves_timeseries.png")

            # plot 3b: significant wave height only
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(waves["time"], waves["sigwaveheight"], "+", markersize=10, markeredgewidth=2.5)
            ax.set_ylim(0, 2)
            ax.set_ylabel("$H_s$ [m]")
            ax.set_title(f"STBM significant wave height\n{title_range}")
            fig.autofmt_xdate()
            fig.tight_layout()
            fig.savefig(OUTPUT_DIR / "STBM_Sig500s_Rosario_waveheight.png", dpi=150)
            written.append("STBM_Sig500s_Rosario_waveheight.png")

            # plot 4: wave height vs peak period
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.hexbin(waves["peakwaveperiod"], waves["sigwaveheight"], gridsize=40, mincnt=1)
            ax.set_xlabel("Wave T_p [s]")
            ax.set_ylabel("Wave H_s [m]")
            ax.set_title(f"STBM wave height vs. period\n{title_range}")
            fig.tight_layout()
            fig.savefig(OUTPUT_DIR / "STBM_Sig500s_Rosario_waves_hist.png", dpi=150)
            written.append("STBM_Sig500s_Rosario_waves_hist.png")
    else:
        print("'Waves' not found in the .mat file; skipping wave plots.")

    print(f"Wrote {', '.join(written)} to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()
