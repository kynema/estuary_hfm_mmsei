"""
Plot the FVCOM vertical velocity profile at a Rosario deployment location
(Sea Spider or Stablemoor), in the same style as plot_sgf_data.py and
plot_adcp_data.py.

This is the script version of the "### Plot temporal evolution of flow in
deployment location" cell in
../../FVCOM-dataExtraction/FVCOM-dataExtraction.ipynb: it loads the reduced
FVCOM numpy arrays, finds the FVCOM grid column nearest the deployment
location, and converts the sigma-layer data to heights above the seafloor so
it can be plotted (and optionally overlaid with measured ADCP data) exactly
like the SGF line-sampling data.

Produces:
  1. Average speed profile vs. height above seafloor (z=0 at the seafloor).
  2. Time series of speed at each height in time_series_heights (plus
     optionally the column average).

Edit fvcom_unit, REDUCED_DATA_DIR and the time-window settings below. FVCOM
time is a Modified Julian Day (days since 1858-11-17 00:00 UTC), so the
window is specified the same way as in the notebook.

FVCOM here is a 2015 hindcast while the ADCP measurements are from 2024, so
absolute timestamps cannot be compared. If plot_adcp_data is True, the FVCOM
window start is mapped onto adcp_start_date/adcp_start_time and the ADCP
record over the matching elapsed interval is overlaid as 'o' markers.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import adcp_data_extraction
import fvcom_data_extraction

OUTPUT_DIR = Path(__file__).parent / "figures"

# ---- FVCOM deployment location ----
# "ss" for Sea Spider, "stbm" for Stablemoor
fvcom_unit = "ss"

# ---- data location ----
# Directory holding the reduced FVCOM arrays (PS_*.npy)
case_name = "PS_20150331_20150501"
REDUCED_DATA_DIR = Path(
    "/projects/hfm/churchfield/rosario/FVCOM/WA_puget_sound/v1.0.0/00_raw/"
    "Puget_Sound_corrected/20150331_20150501/reduced"
)

# ---- FVCOM time window (Modified Julian Day) ----
fvcom_start_mjd = 5.7023e4 + 31.0 + 28.0 + 31.0 + (5 - 1) + (11.0 + 30.0) / 24.0
duration_hours = 24

# ---- ADCP overlay ----
# The comparison workflow applies an 8-hour FVCOM-to-UTC offset. Therefore,
# this raw-MJD window is 38 hours after the validated comparison reference:
# 2015-04-05 11:00 UTC maps to ADCP 2024-10-02 22:00 UTC.
adcp_start_date = "2024-10-04"
adcp_start_time = "12:00:00"
plot_adcp_data = True  # Set to False to skip plotting ADCP data

# Specify heights above seafloor for time series plots.
# Include "avg" to plot the column-averaged speed
time_series_heights = [15, 73, "avg"]


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    use_real_time = plot_adcp_data and adcp_start_date is not None and adcp_start_time is not None
    data = fvcom_data_extraction.load_fvcom_column(
        fvcom_unit,
        REDUCED_DATA_DIR,
        fvcom_start_mjd,
        duration_hours,
        start_date=adcp_start_date if use_real_time else None,
        start_time=adcp_start_time if use_real_time else None,
    )

    print(f"Loaded {data['file_prefix']} data from {REDUCED_DATA_DIR}")

    adcp_data = None
    if plot_adcp_data:
        if not use_real_time:
            print("plot_adcp_data is True but adcp_start_date/adcp_start_time are not set; skipping ADCP overlay.")
        else:
            stop_timestamp = data["scalars"]["time"].max()
            print(f"Loading ADCP data for {fvcom_unit}...")
            print(f"   ADCP data time range: {adcp_start_date} {adcp_start_time} to {stop_timestamp:%Y-%m-%d %H:%M:%S}")
            adcp_data = adcp_data_extraction.load_adcp_data(
                fvcom_unit,
                start_date=adcp_start_date,
                stop_date=stop_timestamp.strftime("%Y-%m-%d"),
                start_time=adcp_start_time,
                stop_time=stop_timestamp.strftime("%H:%M:%S"),
            )

    # ---- average speed profile vs height above seafloor ----
    # The sigma layers move with the tide, so average per sigma level and use
    # each level's time-averaged height.
    by_level = data["profiles"].groupby("sigma_level")
    avg_profile = by_level["speed"].mean()
    avg_z = by_level["z_m"].mean()
    order = np.argsort(avg_z.values)

    fig, ax = plt.subplots(figsize=(4.5, 8))
    ax.plot(avg_profile.values[order], avg_z.values[order], "b-", linewidth=3, label=data["plot_label"])
    if adcp_data is not None:
        adcp_avg_profile = adcp_data["profiles"].groupby("z_m")["speed"].mean().sort_index()
        ax.plot(adcp_avg_profile.values, adcp_avg_profile.index, "o", markerfacecolor="none", label="ADCP")
        ax.legend()
    ax.set_xlabel("avg speed [m/s]")
    ax.set_ylabel("height above seafloor [m]")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title(f"{data['plot_label']} avg speed profile\n{data['title_range']}")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_profile_{case_name}.png", dpi=150)

    # ---- time series of speed at each height in time_series_heights ----
    adcp_available_z = adcp_data["profiles"]["z_m"].unique() if adcp_data is not None else None

    fig, ax = plt.subplots(figsize=(10, 4))
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, height in enumerate(time_series_heights):
        color = color_cycle[i % len(color_cycle)]
        if height == "avg":
            column_avg = data["profiles"].groupby("time")["speed"].mean().sort_index()
            ax.plot(column_avg.index, column_avg.values, label="column average", linewidth=3, color=color)
            if adcp_data is not None:
                adcp_column_avg = adcp_data["profiles"].groupby("time")["speed"].mean().sort_index()
                ax.plot(adcp_column_avg.index, adcp_column_avg.values, "o", markerfacecolor="none",
                        label="ADCP column average", color=color)
        else:
            level = int(avg_z.index[np.argmin(np.abs(avg_z.values - height))])
            series = data["profiles"][data["profiles"]["sigma_level"] == level].sort_values("time")
            ax.plot(series["time"], series["speed"], linewidth=3, color=color,
                    label=f"z = {avg_z.loc[level]:.1f} m above seafloor")
            if adcp_data is not None:
                adcp_nearest_z = adcp_available_z[np.argmin(np.abs(adcp_available_z - height))]
                adcp_series = adcp_data["profiles"][adcp_data["profiles"]["z_m"] == adcp_nearest_z].sort_values("time")
                ax.plot(adcp_series["time"], adcp_series["speed"], "o", markerfacecolor="none",
                        label=f"ADCP z = {adcp_nearest_z:.1f} m", color=color)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("speed [m/s]")
    ax.set_xlabel("time" if adcp_data is not None else "time since window start [s]")
    ax.set_title(f"{data['plot_label']} speed time series\n{data['title_range']}")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 0.90, 1))
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_timeseries_{case_name}.png", dpi=150)

    print(f"Wrote {data['file_prefix']}_profile_{case_name}.png and "
          f"{data['file_prefix']}_timeseries_{case_name}.png to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()