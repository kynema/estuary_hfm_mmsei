"""
Plot Rosario ADCP data (Sea Spider or Stablemoor) for a user-selectable date
range and ADCP unit.

Produces:
  1. Average speed profile vs. height above seafloor (z=0 at seafloor). A
     fixed surface-reflection zone (z=83 to 94 m above seafloor, matching the
     area([...],[94 94],83) shading in the original Matlab scripts) is shaded
     red and excluded from the plotted profile line.
  2. Time series of speed at each height in time_series_heights (heights
     above the seafloor, plus optionally the column average).
  3. For STBM only, if wave data is present: surface wave time series
     (significant wave height, peak wave period) and a wave height vs. peak
     period scatter, matching the fh21/fh22 plots in
     Matlab/RosarioStablemoorAnalysis.m.

Edit adcp_unit, start_date, stop_date, and start_time below to select the
ADCP unit and narrow the analysis window.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import adcp_data_extraction

OUTPUT_DIR = Path(__file__).parent / "figures"

# ---- user-editable ADCP unit ----
# "ss" or "SS" for Sea Spider, "stbm" or "STBM" for Stablemoor
adcp_unit = "stbm"

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

# ---- turbine / site params (Matlab/RosarioSeaSpiderAnalysis.m, Matlab/RosarioStablemoorAnalysis.m) ----
R = 13.5  # rotor radius (m)
hubdepth = 3.5 + R  # hub depth below surface (m)
waterdepth = adcp_data_extraction.STBM_WATERDEPTH  # m, surveyed depth (not derivable from STBM pressure)

# height above seafloor to report a time series for, in addition to hub depth
near_bottom_height = 3.0  # m

# Specify heights above seafloor for time series plots.
# Include "avg" to plot the column-averaged speed
#time_series_heights = [waterdepth-hubdepth, near_bottom_height]  # heights above seafloor for time series
time_series_heights = [waterdepth - hubdepth, near_bottom_height, "avg"]  # include column average


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    data = adcp_data_extraction.load_adcp_data(adcp_unit, start_date=start_date, stop_date=stop_date, start_time=start_time)

    # ---- plot 1: average speed profile vs height above seafloor ----
    avg_profile = data["profiles"].groupby("z_m")["speed"].mean().sort_index()

    hub_z = adcp_data_extraction.hub_height_above_seafloor(data, hubdepth, waterdepth)
    reflection_zone_start, reflection_zone_end = adcp_data_extraction.reflection_zone()

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
    if data["unit"] == "STBM":
        ax.set_ylim(0, 100)
    ax.set_title(f"{data['plot_label']} avg speed profile\n{data['title_range']}")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_profile.png", dpi=150)

    # ---- plot 2: time series of speed at each height in time_series_heights ----
    available_z = data["profiles"]["z_m"].unique()

    fig, ax = plt.subplots(figsize=(10, 4))
    for h in time_series_heights:
        if h == "avg":
            column_avg = data["profiles"].groupby("time")["speed"].mean().sort_index()
            ax.plot(column_avg.index, column_avg.values, label="column average", linewidth=3)
        else:
            nearest_z = available_z[np.argmin(np.abs(available_z - h))]
            series = data["profiles"][data["profiles"]["z_m"] == nearest_z].sort_values("time")
            ax.plot(series["time"], series["speed"], label=f"z = {nearest_z:.1f} m above seafloor", linewidth=3)
    ax.set_ylabel("speed [m/s]")
    ax.set_title(f"{data['plot_label']} speed time series\n{data['title_range']}")
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_timeseries.png", dpi=150)

    written = [f"{data['file_prefix']}_profile.png", f"{data['file_prefix']}_timeseries.png"]

    # ---- plots 3 & 4: surface waves (STBM only, if 'Waves' is present) ----
    if data["unit"] == "STBM":
        if data["waves"] is not None:
            # plot 3: wave height / peak period time series
            fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
            axes[0].plot(data["waves"]["time"], data["waves"]["sigwaveheight"], "+", markersize=10, markeredgewidth=2.5)
            axes[0].set_ylim(0, 2)
            axes[0].set_ylabel("$H_s$ [m]")
            axes[1].plot(data["waves"]["time"], data["waves"]["peakwaveperiod"], "+", markersize=10, markeredgewidth=2.5)
            axes[1].set_ylim(0, 10)
            axes[1].set_ylabel("$T_p$ [s]")
            fig.autofmt_xdate()
            fig.suptitle(f"{data['plot_label']} surface waves\n{data['title_range']}")
            fig.tight_layout()
            fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_waves_timeseries.png", dpi=150)
            written.append(f"{data['file_prefix']}_waves_timeseries.png")

            # plot 3b: significant wave height only
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(data["waves"]["time"], data["waves"]["sigwaveheight"], "+", markersize=10, markeredgewidth=2.5)
            ax.set_ylim(0, 2)
            ax.set_ylabel("$H_s$ [m]")
            ax.set_title(f"{data['plot_label']} significant wave height\n{data['title_range']}")
            fig.autofmt_xdate()
            fig.tight_layout()
            fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_waveheight.png", dpi=150)
            written.append(f"{data['file_prefix']}_waveheight.png")

            # plot 4: wave height vs peak period
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.hexbin(data["waves"]["peakwaveperiod"], data["waves"]["sigwaveheight"], gridsize=40, mincnt=1)
            ax.set_xlabel("Wave T_p [s]")
            ax.set_ylabel("Wave H_s [m]")
            ax.set_title(f"{data['plot_label']} wave height vs. period\n{data['title_range']}")
            fig.tight_layout()
            fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_waves_hist.png", dpi=150)
            written.append(f"{data['file_prefix']}_waves_hist.png")
        else:
            print("No wave data available (missing 'Waves' struct or no data in range); skipping wave plots.")

    print(f"Wrote {', '.join(written)} to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()
