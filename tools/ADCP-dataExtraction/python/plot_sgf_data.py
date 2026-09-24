"""
Plot Kynema-SGF Rosario line-sampling data (Sea Spider or Stablemoor) for a
user-selectable SGF unit, mirroring plot_adcp_data.py's average speed profile
plot (z=0 at seafloor, after sgf_data_extraction.py shifts the SGF z=0 water
surface datum to seafloor datum).

Edit adcp_unit, SGF_DATA_DIR, start_date, and start_time below to select the
SGF unit, the simulation's post_processing directory, and the shared native
UTC datetime corresponding to SGF/FVCOM simulation time=0. Configure
adcp_start_date/adcp_start_time separately; ADCP values are shifted by elapsed
time onto the shared timeline for combined plots.

Optionally set end_date/end_time to truncate the SGF data to a fixed
real-world end datetime (e.g. to compare multiple cases over the same
window, regardless of how long each case actually ran). If left as None,
the full extent of the loaded SGF data is used.

If plot_adcp_data is True, the measured ADCP data for the same elapsed window
is loaded using adcp_start_date/adcp_start_time and overlaid as 'o' markers on
both the profile and time series plots.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import adcp_data_extraction
import fvcom_data_extraction
import sgf_data_extraction

OUTPUT_DIR = Path(__file__).parent / "figures"

plot_sgf_data = True  # Set to False to skip plotting SGF data
plot_adcp_data = True  # Set to False to skip plotting ADCP data
plot_fvcom_data = True  # Set to True to overlay FVCOM data

# ---- ADCP unit ----
# "ss" for Sea Spider, "stbm" for Stablemoor
adcp_unit = "ss"

# ---- SGF data location ----
# Directory containing the line_sampling##### folders 
case_name = "max_lev4"
SGF_DATA_DIR = Path(
   f"/scratch/mkuhn/estuary_flows/milestone/{case_name}/post_processing"
)

# ---- Shared SGF/FVCOM UTC time ----
# All user-entered times in this script are UTC. The raw FVCOM MJD is derived
# from start_date/start_time with the 8-hour FVCOM-to-UTC conversion.
start_date = "2015-04-07"
start_time = "01:00:00"
end_date = None
end_time = None
#end_date = "2015-04-07"
#end_time = "06:52:00"

# ---- ADCP time mapped to the shared SGF/FVCOM start ----
# ADCP timestamps are shifted onto the shared timeline for combined plots.
adcp_start_date = "2024-10-04"
adcp_start_time = "12:00:00"

# ---- FVCOM overlay ----
# Directory holding the reduced FVCOM arrays (PS_*.npy). The FVCOM window is
# mapped to start_date/start_time above so it aligns by elapsed time with SGF.
FVCOM_REDUCED_DATA_DIR = Path(
    "/projects/hfm/churchfield/rosario/FVCOM/WA_puget_sound/v1.0.0/00_raw/"
    "Puget_Sound_corrected/20150331_20150501/reduced"
)


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

    if not any((plot_sgf_data, plot_fvcom_data, plot_adcp_data)):
        raise ValueError("Enable at least one of plot_sgf_data, plot_fvcom_data, or plot_adcp_data.")

    # SGF data is still loaded when disabled from the plot because it defines
    # the comparison window when end_date/end_time are not configured.
    data = sgf_data_extraction.load_sgf_data(
        adcp_unit, SGF_DATA_DIR, start_date=start_date, start_time=start_time
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

    comparison_stop_timestamp = (
        pd.Timestamp(f"{end_date} {end_time}")
        if end_date is not None and end_time is not None and start_date is not None and start_time is not None
        else data["scalars"]["time"].max()
    )

    fvcom_data = None
    if plot_fvcom_data:
        if start_date is None or start_time is None:
            print("plot_fvcom_data is True but start_date/start_time are not set; skipping FVCOM overlay.")
        else:
            fvcom_start_mjd, fvcom_duration_hours, _ = (
                fvcom_data_extraction.resolve_fvcom_window(
                    start_date=start_date,
                    start_time=start_time,
                    end_date=comparison_stop_timestamp.strftime("%Y-%m-%d"),
                    end_time=comparison_stop_timestamp.strftime("%H:%M:%S"),
                )
            )
            print(f"Loading {fvcom_duration_hours:.2f} hours of FVCOM data...")
            fvcom_data = fvcom_data_extraction.load_fvcom_column(
                adcp_unit,
                FVCOM_REDUCED_DATA_DIR,
                fvcom_start_mjd,
                fvcom_duration_hours,
                start_date=start_date,
                start_time=start_time,
            )
            fvcom_data["scalars"] = fvcom_data["scalars"].loc[
                fvcom_data["scalars"]["time"] <= comparison_stop_timestamp
            ].reset_index(drop=True)
            fvcom_data["profiles"] = fvcom_data["profiles"].loc[
                fvcom_data["profiles"]["time"] <= comparison_stop_timestamp
            ].reset_index(drop=True)
            if fvcom_data["scalars"].empty or fvcom_data["profiles"].empty:
                raise ValueError("No FVCOM data remains before the comparison end time.")
            fvcom_data["time"] = fvcom_data["scalars"]["time"].to_numpy()
            fvcom_data["title_range"] = (
                f"{fvcom_data['scalars']['time'].min():%Y-%m-%d %H:%M:%S} to "
                f"{fvcom_data['scalars']['time'].max():%Y-%m-%d %H:%M:%S}"
            )
            print(f"Loaded {fvcom_data['file_prefix']} data from {FVCOM_REDUCED_DATA_DIR}")

    adcp_data = None
    if plot_adcp_data:
        if start_date is None or start_time is None or adcp_start_date is None or adcp_start_time is None:
            print("plot_adcp_data is True but a shared or ADCP start timestamp is not set; skipping ADCP overlay.")
        else:
            shared_start_timestamp = pd.Timestamp(f"{start_date} {start_time}")
            adcp_start_timestamp = pd.Timestamp(f"{adcp_start_date} {adcp_start_time}")
            adcp_stop_timestamp = adcp_start_timestamp + (
                comparison_stop_timestamp - shared_start_timestamp
            )
            print(f"Loading ADCP data for {adcp_unit}...")
            print(f"   ADCP data time range: {adcp_start_timestamp:%Y-%m-%d %H:%M:%S} to {adcp_stop_timestamp:%Y-%m-%d %H:%M:%S}")
            adcp_data = adcp_data_extraction.load_adcp_data(
                adcp_unit,
                start_date=adcp_start_timestamp.strftime("%Y-%m-%d"),
                stop_date=adcp_stop_timestamp.strftime("%Y-%m-%d"),
                start_time=adcp_start_timestamp.strftime("%H:%M:%S"),
                stop_time=adcp_stop_timestamp.strftime("%H:%M:%S"),
            )
            adcp_data_extraction.map_adcp_to_timeline(
                adcp_data, shared_start_timestamp, adcp_start_timestamp
            )

    plot_labels = []
    plot_title_range = data["title_range"]
    if plot_sgf_data:
        plot_labels.append(data["plot_label"])
    if fvcom_data is not None:
        plot_labels.append(fvcom_data["plot_label"])
        if not plot_sgf_data:
            plot_title_range = fvcom_data["title_range"]
    if adcp_data is not None:
        plot_labels.append("ADCP")
        if not plot_sgf_data and fvcom_data is None:
            plot_title_range = adcp_data["title_range"]
    if not plot_labels:
        raise ValueError("No requested dataset could be plotted with the current time configuration.")
    output_prefix = (
        data["file_prefix"] if plot_sgf_data else
        fvcom_data["file_prefix"] if fvcom_data is not None else
        f"ADCP_{adcp_data['unit']}"
    )

    # ---- average speed profile vs height above seafloor ----
    fig, ax = plt.subplots(figsize=(5.25, 8))
    if plot_sgf_data:
        avg_profile = data["profiles"].groupby("z_m")["speed"].mean().sort_index()
        ax.plot(avg_profile.values, avg_profile.index, "-", linewidth=3, color="tab:blue", label=data["plot_label"])
    if fvcom_data is not None:
        fvcom_by_level = fvcom_data["profiles"].groupby("sigma_level")
        fvcom_avg_profile = fvcom_by_level["speed"].mean()
        fvcom_avg_z = fvcom_by_level["z_m"].mean()
        fvcom_order = np.argsort(fvcom_avg_z.values)
        ax.plot(
            fvcom_avg_profile.values[fvcom_order],
            fvcom_avg_z.values[fvcom_order],
            linestyle=":",
            linewidth=3,
            color="tab:orange",
            label=fvcom_data["plot_label"],
        )
    if adcp_data is not None:
        adcp_avg_profile = adcp_data["profiles"].groupby("z_m")["speed"].mean().sort_index()
        ax.plot(adcp_avg_profile.values, adcp_avg_profile.index, "o", markerfacecolor="none", color="black", label="ADCP")
    if fvcom_data is not None or adcp_data is not None:
        ax.legend()
    ax.set_xlabel("avg speed [m/s]", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel("height above seafloor [m]", fontsize=Y_LABEL_FONT_SIZE)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title(plot_title_range, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout()
    figname = f"{output_prefix}_profile_{case_name}.png"
    fig.savefig(OUTPUT_DIR / figname, dpi=150)

    # ---- time series of speed at each height in time_series_heights ----
    available_z = data["profiles"]["z_m"].unique() if plot_sgf_data else None
    adcp_available_z = adcp_data["profiles"]["z_m"].unique() if adcp_data is not None else None
    if fvcom_data is not None:
        fvcom_by_level = fvcom_data["profiles"].groupby("sigma_level")
        fvcom_avg_z = fvcom_by_level["z_m"].mean()

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
    colors = {"sgf": "tab:blue", "fvcom": "tab:orange", "adcp": "black"}

    for ax, h in zip(axes, time_series_heights):
        if h == "avg":
            if plot_sgf_data:
                column_avg = data["profiles"].groupby("time")["speed"].mean().sort_index()
                ax.plot(column_avg.index, column_avg.values, label="SGF", linewidth=3, color=colors["sgf"])
            if fvcom_data is not None:
                fvcom_column_avg = fvcom_data["profiles"].groupby("time")["speed"].mean().sort_index()
                ax.plot(fvcom_column_avg.index, fvcom_column_avg.values, linestyle=":", linewidth=3,
                        label="FVCOM", color=colors["fvcom"])
            if adcp_data is not None:
                adcp_column_avg = adcp_data["profiles"].groupby("time")["speed"].mean().sort_index()
                ax.plot(adcp_column_avg.index, adcp_column_avg.values, "o", markerfacecolor="none",
                        label="ADCP", color=colors["adcp"])
            subplot_label = "column average"
        else:
            height_labels = []
            if plot_sgf_data:
                nearest_z = available_z[np.argmin(np.abs(available_z - h))]
                series = data["profiles"][data["profiles"]["z_m"] == nearest_z].sort_values("time")
                ax.plot(series["time"], series["speed"], label="SGF", linewidth=3, color=colors["sgf"])
                height_labels.append(f"SGF {nearest_z:.1f} m")
            if fvcom_data is not None:
                fvcom_level = int(fvcom_avg_z.index[np.argmin(np.abs(fvcom_avg_z.values - h))])
                fvcom_series = fvcom_data["profiles"][
                    fvcom_data["profiles"]["sigma_level"] == fvcom_level
                ].sort_values("time")
                ax.plot(fvcom_series["time"], fvcom_series["speed"], linestyle=":", linewidth=3,
                        label="FVCOM", color=colors["fvcom"])
                height_labels.append(f"FVCOM {fvcom_avg_z.loc[fvcom_level]:.1f} m")
            if adcp_data is not None:
                adcp_nearest_z = adcp_available_z[np.argmin(np.abs(adcp_available_z - h))]
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
    axes[-1].set_xlabel("time", fontsize=X_LABEL_FONT_SIZE)
    fig.suptitle(plot_title_range, fontsize=TITLE_FONT_SIZE)
    fig.autofmt_xdate()
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    figname = f"{output_prefix}_timeseries_{case_name}.png"
    fig.savefig(OUTPUT_DIR / figname, dpi=150)

    print(
        f"Wrote {output_prefix}_profile_{case_name}.png and "
        f"{output_prefix}_timeseries_{case_name}.png to {OUTPUT_DIR}"
    )
    plt.show()


if __name__ == "__main__":
    main()
