"""
Plot Kynema-SGF Rosario line-sampling data (Sea Spider or Stablemoor) for a
user-selectable SGF unit, mirroring plot_adcp_data.py's average speed profile
plot (z=0 at seafloor, after sgf_data_extraction.py shifts the SGF z=0 water
surface datum to seafloor datum).

Edit sgf_unit, SGF_DATA_DIR, start_date, and start_time below to select the
SGF unit, the post_processing directory to load, and (optionally) the
real-world UTC datetime corresponding to simulation time=0.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import sgf_data_extraction

OUTPUT_DIR = Path(__file__).parent / "figures"

# ---- SGF ADCP unit ----
# "ss" for Sea Spider, "stbm" for Stablemoor
sgf_unit = "ss"

# ---- data location ----
# Directory containing the line_sampling##### folders 
SGF_DATA_DIR = Path(
    "/projects/hfm/dmontgomery/1st_run_static_surf_refinement_fsdamp_1/post_processing"
)

# ---- user-editable simulation start time (real-world UTC) ----
# Set both to convert the "time" column from raw simulation seconds to real
# datetimes, comparable to the ADCP timestamps. Leave as None to keep "time"
# as raw simulation seconds.
start_date = None  # e.g. "2024-10-22"
start_time = None  # e.g. "22:00:00"

# Specify heights above seafloor for time series plots.
# Include "avg" to plot the column-averaged speed
time_series_heights = [30, 60, "avg"]


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    data = sgf_data_extraction.load_sgf_data(
        sgf_unit, SGF_DATA_DIR, start_date=start_date, start_time=start_time
    )

    # ---- average speed profile vs height above seafloor ----
    avg_profile = data["profiles"].groupby("z_m")["speed"].mean().sort_index()

    fig, ax = plt.subplots(figsize=(3, 8))
    ax.plot(avg_profile.values, avg_profile.index, "b-", linewidth=3)
    ax.set_xlabel("avg speed [m/s]")
    ax.set_ylabel("height above seafloor [m]")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title(f"{data['plot_label']} avg speed profile\n{data['title_range']}")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_profile.png", dpi=150)

    # ---- time series of speed at each height in time_series_heights ----
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
    ax.set_ylim(bottom=0)
    ax.set_ylabel("speed [m/s]")
    ax.set_title(f"{data['plot_label']} speed time series\n{data['title_range']}")
    ax.legend()
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"{data['file_prefix']}_timeseries.png", dpi=150)

    print(f"Wrote {data['file_prefix']}_profile.png and {data['file_prefix']}_timeseries.png to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()
