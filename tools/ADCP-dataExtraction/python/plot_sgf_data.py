"""
Plot Kynema-SGF Rosario line-sampling data (Sea Spider or Stablemoor) for a
user-selectable SGF unit, mirroring plot_adcp_data.py's average speed profile
plot (z=0 at seafloor, after sgf_data_extraction.py shifts the SGF z=0 water
surface datum to seafloor datum).

Edit sgf_unit, SGF_DATA_DIR, start_date, and start_time below to select the
SGF unit, the post_processing directory to load, and (optionally) the
real-world UTC datetime corresponding to simulation time=0.

Optionally set end_date/end_time to truncate the SGF data to a fixed
real-world end datetime (e.g. to compare multiple cases over the same
window, regardless of how long each case actually ran). If left as None,
the full extent of the loaded SGF data is used.

If plot_adcp_data is True (requires start_date/start_time to be set so the
SGF sim time can be mapped to real-world time), the measured ADCP data for
the same unit is loaded for the real-world time interval spanned by the SGF
data (start_date/start_time through end_date/end_time if set, else through a
stop date/time derived from the SGF data's max simulation time) and overlaid
as 'o' markers on both the profile and time series plots.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import adcp_data_extraction
import sgf_data_extraction

OUTPUT_DIR = Path(__file__).parent / "figures"

# ---- SGF ADCP unit ----
# "ss" for Sea Spider, "stbm" for Stablemoor
sgf_unit = "stbm"

# ---- data location ----
# Directory containing the line_sampling##### folders 
case_name = "max_lev3"
SGF_DATA_DIR = Path(
   f"/scratch/mkuhn/estuary_flows/milestone/{case_name}/post_processing"
)

# ---- times ----
# start_date = None  # e.g. "2024-10-22", required for plotting ADCP data
# start_time = None  # e.g. "22:00:00", required for plotting ADCP data
# end_date = None  # optional e.g. "2024-10-23"
# end_time = None  # optional e.g. "10:00:00"
start_date = "2024-10-22"
start_time = "22:00:00"
end_date = "2024-10-23"
end_time = "01:18:00"
plot_adcp_data = True  # Set to False to skip plotting ADCP data

# Specify heights above seafloor for time series plots.
# Include "avg" to plot the column-averaged speed
time_series_heights = [15, 73, "avg"]


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    data = sgf_data_extraction.load_sgf_data(
        sgf_unit, SGF_DATA_DIR, start_date=start_date, start_time=start_time
    )

    print(f"Loaded {data['file_prefix']} data from {SGF_DATA_DIR}")

    if end_date is not None and end_time is not None:
        if start_date is None or start_time is None:
            print("end_date/end_time are set but start_date/start_time are not; ignoring end_date/end_time.")
        else:
            end_timestamp = pd.Timestamp(f"{end_date} {end_time}")
            data["scalars"] = data["scalars"][data["scalars"]["time"] <= end_timestamp].reset_index(drop=True)
            data["profiles"] = data["profiles"][data["profiles"]["time"] <= end_timestamp].reset_index(drop=True)
            if data["scalars"].empty or data["profiles"].empty:
                raise ValueError("No SGF data remains after applying end_date/end_time.")
            data["title_range"] = (
                f"{data['scalars']['time'].min():%Y-%m-%d %H:%M:%S} to "
                f"{data['scalars']['time'].max():%Y-%m-%d %H:%M:%S}"
            )

    adcp_data = None
    if plot_adcp_data:
        if start_date is None or start_time is None:
            print("plot_adcp_data is True but start_date/start_time are not set; skipping ADCP overlay.")
        else:
            stop_timestamp = data["scalars"]["time"].max()
            print(f"Loading ADCP data for {sgf_unit}...")
            print(f"   ADCP data time range: {start_date} {start_time} to {stop_timestamp.strftime('%Y-%m-%d %H:%M:%S')}")
            adcp_data = adcp_data_extraction.load_adcp_data(
                sgf_unit,
                start_date=start_date,
                stop_date=stop_timestamp.strftime("%Y-%m-%d"),
                start_time=start_time,
                stop_time=stop_timestamp.strftime("%H:%M:%S"),
            )

    # ---- average speed profile vs height above seafloor ----
    avg_profile = data["profiles"].groupby("z_m")["speed"].mean().sort_index()

    fig, ax = plt.subplots(figsize=(4.5, 8))
    ax.plot(avg_profile.values, avg_profile.index, "b-", linewidth=3, label=data["plot_label"])
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
    figname = f"{data['file_prefix']}_profile_{case_name}.png"
    fig.savefig(OUTPUT_DIR / figname, dpi=150)

    # ---- time series of speed at each height in time_series_heights ----
    available_z = data["profiles"]["z_m"].unique()
    adcp_available_z = adcp_data["profiles"]["z_m"].unique() if adcp_data is not None else None

    fig, ax = plt.subplots(figsize=(10, 4))
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, h in enumerate(time_series_heights):
        color = color_cycle[i % len(color_cycle)]
        if h == "avg":
            column_avg = data["profiles"].groupby("time")["speed"].mean().sort_index()
            ax.plot(column_avg.index, column_avg.values, label="column average", linewidth=3, color=color)
            if adcp_data is not None:
                adcp_column_avg = adcp_data["profiles"].groupby("time")["speed"].mean().sort_index()
                ax.plot(adcp_column_avg.index, adcp_column_avg.values, "o", markerfacecolor="none", label="ADCP column average", color=color)
        else:
            nearest_z = available_z[np.argmin(np.abs(available_z - h))]
            series = data["profiles"][data["profiles"]["z_m"] == nearest_z].sort_values("time")
            ax.plot(series["time"], series["speed"], label=f"z = {nearest_z:.1f} m above seafloor", linewidth=3, color=color)
            if adcp_data is not None:
                adcp_nearest_z = adcp_available_z[np.argmin(np.abs(adcp_available_z - h))]
                adcp_series = adcp_data["profiles"][adcp_data["profiles"]["z_m"] == adcp_nearest_z].sort_values("time")
                ax.plot(adcp_series["time"], adcp_series["speed"], "o", markerfacecolor="none", label=f"ADCP z = {adcp_nearest_z:.1f} m", color=color)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("speed [m/s]")
    ax.set_title(f"{data['plot_label']} speed time series\n{data['title_range']}")
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    figname = f"{data['file_prefix']}_timeseries_{case_name}.png"
    fig.savefig(OUTPUT_DIR / figname, dpi=150)

    print(f"Wrote {data['file_prefix']}_profile.png and {data['file_prefix']}_timeseries.png to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()
