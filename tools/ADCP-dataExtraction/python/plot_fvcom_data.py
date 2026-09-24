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

Edit fvcom_unit, REDUCED_DATA_DIR and the time-window settings below. The
FVCOM window can be specified either with native Modified Julian Day (MJD)
values or with UTC start/end dates and times, but never both at once.

FVCOM here is a 2015 hindcast while the ADCP measurements are from 2024, so
absolute timestamps cannot be compared. If plot_adcp_data is True, the FVCOM
window start is mapped onto adcp_start_date/adcp_start_time and the ADCP
record over the matching elapsed interval is overlaid as 'o' markers.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

# ---- FVCOM time window ----
# Choose exactly one configuration:
#   - MJD: set fvcom_start_mjd and either fvcom_end_mjd or duration_hours.
#   - UTC: set start_date, start_time, end_date, and end_time, and set both
#          MJD values to None. All UTC values are interpreted as UTC.
#
# FVCOM MJD is days since 1858-11-17 00:00 UTC. The UTC configuration applies
# the FVCOM-to-UTC offset used by fvcom_data_extraction automatically.
fvcom_start_mjd = None
#fvcom_start_mjd = 5.7023e4 + 31.0 + 28.0 + 31.0 + (5 - 1) + (11.0 + 30.0) / 24.0
fvcom_end_mjd = None
duration_hours = 24  # Used only with MJD when fvcom_end_mjd is None.

# UTC
#start_date = None  # e.g. "2015-04-07"
#start_time = None  # e.g. "01:00:00"
#end_date = None    # e.g. "2015-04-08"
#end_time = None    # e.g. "01:00:00"
start_date = "2015-04-07"
start_time = "01:00:00"
end_date = "2015-04-07" 
end_time = "13:00:00"    

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

# ---- Plot font sizes ----
X_LABEL_FONT_SIZE = 16
Y_LABEL_FONT_SIZE = 16
TICK_LABEL_FONT_SIZE = 14
TITLE_FONT_SIZE = 14


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fvcom_start_mjd_resolved, duration_hours_resolved, utc_start = (
        fvcom_data_extraction.resolve_fvcom_window(
            start_mjd=fvcom_start_mjd,
            end_mjd=fvcom_end_mjd,
            duration_hours=duration_hours,
            start_date=start_date,
            start_time=start_time,
            end_date=end_date,
            end_time=end_time,
        )
    )
    use_real_time = plot_adcp_data and adcp_start_date is not None and adcp_start_time is not None
    plot_start = utc_start if utc_start is not None else (
        pd.Timestamp(f"{adcp_start_date} {adcp_start_time}") if use_real_time else None
    )
    data = fvcom_data_extraction.load_fvcom_column(
        fvcom_unit,
        REDUCED_DATA_DIR,
        fvcom_start_mjd_resolved,
        duration_hours_resolved,
        start_date=plot_start.strftime("%Y-%m-%d") if plot_start is not None else None,
        start_time=plot_start.strftime("%H:%M:%S") if plot_start is not None else None,
    )

    print(f"Loaded {data['file_prefix']} data from {REDUCED_DATA_DIR}")

    adcp_data = None
    if plot_adcp_data:
        if not use_real_time:
            print("plot_adcp_data is True but adcp_start_date/adcp_start_time are not set; skipping ADCP overlay.")
        else:
            adcp_start_timestamp = pd.Timestamp(f"{adcp_start_date} {adcp_start_time}")
            elapsed_seconds = data["scalars"]["elapsed_s"].max()
            stop_timestamp = adcp_start_timestamp + pd.to_timedelta(elapsed_seconds, unit="s")
            print(f"Loading ADCP data for {fvcom_unit}...")
            print(f"   ADCP data time range: {adcp_start_date} {adcp_start_time} to {stop_timestamp:%Y-%m-%d %H:%M:%S}")
            adcp_data = adcp_data_extraction.load_adcp_data(
                fvcom_unit,
                start_date=adcp_start_date,
                stop_date=stop_timestamp.strftime("%Y-%m-%d"),
                start_time=adcp_start_time,
                stop_time=stop_timestamp.strftime("%H:%M:%S"),
            )
            if utc_start is not None:
                adcp_data_extraction.map_adcp_to_timeline(
                    adcp_data, utc_start, adcp_start_timestamp
                )

    # ---- average speed profile vs height above seafloor ----
    # The sigma layers move with the tide, so average per sigma level and use
    # each level's time-averaged height.
    by_level = data["profiles"].groupby("sigma_level")
    avg_profile = by_level["speed"].mean()
    avg_z = by_level["z_m"].mean()
    order = np.argsort(avg_z.values)

    fig, ax = plt.subplots(figsize=(4.5, 8))
    ax.plot(avg_profile.values[order], avg_z.values[order], ":", linewidth=3, color="tab:orange", label=data["plot_label"])
    if adcp_data is not None:
        adcp_avg_profile = adcp_data["profiles"].groupby("z_m")["speed"].mean().sort_index()
        ax.plot(adcp_avg_profile.values, adcp_avg_profile.index, "o", markerfacecolor="none", color="black", label="ADCP")
        ax.legend()
    ax.set_xlabel("avg speed [m/s]", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel("height above seafloor [m]", fontsize=Y_LABEL_FONT_SIZE)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title(data["title_range"], fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_profile_{case_name}.png", dpi=150)

    # ---- time series of speed at each height in time_series_heights ----
    adcp_available_z = adcp_data["profiles"]["z_m"].unique() if adcp_data is not None else None

    if not time_series_heights:
        raise ValueError("time_series_heights must contain at least one height or 'avg'.")

    fig, axes = plt.subplots(
        len(time_series_heights),
        1,
        figsize=(10, 2.75 * len(time_series_heights)),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    colors = {"fvcom": "tab:orange", "adcp": "black"}

    for ax, height in zip(axes, time_series_heights):
        if height == "avg":
            column_avg = data["profiles"].groupby("time")["speed"].mean().sort_index()
            ax.plot(column_avg.index, column_avg.values, ":", label="FVCOM", linewidth=3, color=colors["fvcom"])
            if adcp_data is not None:
                adcp_column_avg = adcp_data["profiles"].groupby("time")["speed"].mean().sort_index()
                ax.plot(adcp_column_avg.index, adcp_column_avg.values, "o", markerfacecolor="none",
                        label="ADCP", color=colors["adcp"])
            subplot_label = "column average"
        else:
            level = int(avg_z.index[np.argmin(np.abs(avg_z.values - height))])
            series = data["profiles"][data["profiles"]["sigma_level"] == level].sort_values("time")
            ax.plot(series["time"], series["speed"], ":", linewidth=3, color=colors["fvcom"], label="FVCOM")
            height_labels = [f"FVCOM {avg_z.loc[level]:.1f} m"]
            if adcp_data is not None:
                adcp_nearest_z = adcp_available_z[np.argmin(np.abs(adcp_available_z - height))]
                adcp_series = adcp_data["profiles"][adcp_data["profiles"]["z_m"] == adcp_nearest_z].sort_values("time")
                ax.plot(adcp_series["time"], adcp_series["speed"], "o", markerfacecolor="none",
                        label="ADCP", color=colors["adcp"])
                height_labels.append(f"ADCP {adcp_nearest_z:.1f} m")
            subplot_label = "; ".join(height_labels)

        ax.set_ylim(bottom=0)
        ax.set_ylabel("speed [m/s]", fontsize=Y_LABEL_FONT_SIZE)
        ax.set_title(subplot_label, fontsize=TITLE_FONT_SIZE)
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONT_SIZE)

    axes[0].legend(loc="upper right")
    axes[-1].set_xlabel(
        "time" if adcp_data is not None else "time since window start [s]",
        fontsize=X_LABEL_FONT_SIZE,
    )
    fig.suptitle(data["title_range"], fontsize=TITLE_FONT_SIZE)
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_timeseries_{case_name}.png", dpi=150)

    print(f"Wrote {data['file_prefix']}_profile_{case_name}.png and "
          f"{data['file_prefix']}_timeseries_{case_name}.png to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()