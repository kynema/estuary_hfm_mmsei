"""Compare multiple Kynema-SGF Rosario cases at one ADCP deployment.

Set ``case_names`` to the SGF cases to compare. Each case must contain a
``post_processing`` directory under ``SGF_CASES_DIR``. All cases share the
same UTC start timestamp, optional UTC end timestamp, and optional FVCOM/ADCP
overlays. When no end timestamp is supplied, the comparison ends at the
shortest loaded SGF case so every SGF curve covers the complete plotted range.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import adcp_data_extraction
import fvcom_data_extraction
import sgf_data_extraction

OUTPUT_DIR = Path(__file__).parent / "figures"

# ---- ADCP unit ----
# "ss" for Sea Spider, "stbm" for Stablemoor
adcp_unit = "ss"

# ---- SGF cases ----
# Each case is loaded from SGF_CASES_DIR / case_name / "post_processing".
SGF_CASES_DIR = Path("/scratch/mkuhn/estuary_flows/milestone")
case_names = ["max_lev1", "max_lev2", "max_lev3", "max_lev4"]

# ---- Shared SGF/FVCOM UTC time ----
# All user-entered times in this script are UTC.
start_date = "2015-04-07"
start_time = "01:00:00"
# Set both values to select a fixed common comparison end. Leave both as None
# to use the final timestamp of the shortest loaded SGF case.
end_date = None
end_time = None

# ---- Optional ADCP overlay ----
plot_adcp_data = True
# ADCP timestamps are mapped by elapsed time onto the shared SGF/FVCOM range.
adcp_start_date = "2024-10-04"
adcp_start_time = "12:00:00"

# ---- Optional FVCOM overlay ----
plot_fvcom_data = True
FVCOM_REDUCED_DATA_DIR = Path(
    "/projects/hfm/churchfield/rosario/FVCOM/WA_puget_sound/v1.0.0/00_raw/"
    "Puget_Sound_corrected/20150331_20150501/reduced"
)

# Specify heights above seafloor for time-series plots. Include "avg" to plot
# the column-averaged speed.
time_series_heights = [15, 73, "avg"]

# ---- Plot font sizes ----
X_LABEL_FONT_SIZE = 16
Y_LABEL_FONT_SIZE = 16
TICK_LABEL_FONT_SIZE = 14
TITLE_FONT_SIZE = 14


def _truncate_to_end(data, end_timestamp, case_name):
    """Restrict a loaded SGF dataset to the shared comparison end."""
    data["scalars"] = data["scalars"].loc[
        data["scalars"]["time"] <= end_timestamp
    ].reset_index(drop=True)
    data["profiles"] = data["profiles"].loc[
        data["profiles"]["time"] <= end_timestamp
    ].reset_index(drop=True)
    if data["scalars"].empty or data["profiles"].empty:
        raise ValueError(f"No SGF data remains for case '{case_name}' in the comparison window.")
    data["time"] = data["scalars"]["time"].to_numpy()


def _set_title_range(data):
    """Update a data dictionary's display range after time filtering."""
    data["title_range"] = (
        f"{data['scalars']['time'].min():%Y-%m-%d %H:%M:%S} to "
        f"{data['scalars']['time'].max():%Y-%m-%d %H:%M:%S}"
    )


def _elapsed_hours(times, start_timestamp):
    """Convert datetime values to elapsed hours from the shared start."""
    return (pd.DatetimeIndex(times) - start_timestamp).total_seconds() / 3600.0


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not case_names:
        raise ValueError("case_names must contain at least one SGF case name.")
    if (end_date is None) != (end_time is None):
        raise ValueError("Set both end_date and end_time, or leave both as None.")

    shared_start = pd.Timestamp(f"{start_date} {start_time}")
    sgf_cases = []
    for case_name in case_names:
        case_dir = SGF_CASES_DIR / case_name / "post_processing"
        data = sgf_data_extraction.load_sgf_data(
            adcp_unit, case_dir, start_date=start_date, start_time=start_time
        )
        print(f"Loaded SGF case '{case_name}' from {case_dir}")
        sgf_cases.append((case_name, data))

    if end_date is not None:
        comparison_stop = pd.Timestamp(f"{end_date} {end_time}")
        if comparison_stop <= shared_start:
            raise ValueError("end_date/end_time must be after start_date/start_time.")
    else:
        comparison_stop = min(data["scalars"]["time"].max() for _, data in sgf_cases)

    for case_name, data in sgf_cases:
        _truncate_to_end(data, comparison_stop, case_name)
        _set_title_range(data)
    comparison_stop = min(data["scalars"]["time"].max() for _, data in sgf_cases)
    title_range = f"{shared_start:%Y-%m-%d %H:%M:%S} to {comparison_stop:%Y-%m-%d %H:%M:%S}"

    fvcom_data = None
    if plot_fvcom_data:
        fvcom_start_mjd, fvcom_duration_hours, _ = fvcom_data_extraction.resolve_fvcom_window(
            start_date=start_date,
            start_time=start_time,
            end_date=comparison_stop.strftime("%Y-%m-%d"),
            end_time=comparison_stop.strftime("%H:%M:%S"),
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
        _truncate_to_end(fvcom_data, comparison_stop, "FVCOM")
        _set_title_range(fvcom_data)

    adcp_data = None
    if plot_adcp_data:
        if adcp_start_date is None or adcp_start_time is None:
            print("plot_adcp_data is True but adcp_start_date/adcp_start_time are not set; skipping ADCP overlay.")
        else:
            adcp_start = pd.Timestamp(f"{adcp_start_date} {adcp_start_time}")
            adcp_stop = adcp_start + (comparison_stop - shared_start)
            print(f"Loading ADCP data for {adcp_unit}...")
            adcp_data = adcp_data_extraction.load_adcp_data(
                adcp_unit,
                start_date=adcp_start.strftime("%Y-%m-%d"),
                stop_date=adcp_stop.strftime("%Y-%m-%d"),
                start_time=adcp_start.strftime("%H:%M:%S"),
                stop_time=adcp_stop.strftime("%H:%M:%S"),
            )
            adcp_data_extraction.map_adcp_to_timeline(adcp_data, shared_start, adcp_start)

    if not time_series_heights:
        raise ValueError("time_series_heights must contain at least one height or 'avg'.")

    sgf_colors = plt.get_cmap("tab10").colors
    case_tag = "_".join(case_names)

    # ---- average speed profile vs. height above seafloor ----
    fig, ax = plt.subplots(figsize=(5.25, 8))
    for index, (case_name, data) in enumerate(sgf_cases):
        avg_profile = data["profiles"].groupby("z_m")["speed"].mean().sort_index()
        ax.plot(
            avg_profile.values,
            avg_profile.index,
            linewidth=3,
            color=sgf_colors[index % len(sgf_colors)],
            label=case_name,
        )
    if fvcom_data is not None:
        by_level = fvcom_data["profiles"].groupby("sigma_level")
        avg_profile = by_level["speed"].mean()
        avg_z = by_level["z_m"].mean()
        order = np.argsort(avg_z.values)
        ax.plot(
            avg_profile.values[order], avg_z.values[order], ":", linewidth=3,
            color="black", label="FVCOM",
        )
    if adcp_data is not None:
        avg_profile = adcp_data["profiles"].groupby("z_m")["speed"].mean().sort_index()
        ax.plot(
            avg_profile.values, avg_profile.index, "o", markerfacecolor="none",
            color="black", label="ADCP",
        )
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.01))
    ax.set_xlabel("avg speed [m/s]", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel("height above seafloor [m]", fontsize=Y_LABEL_FONT_SIZE)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_title(title_range, fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(OUTPUT_DIR / f"multiple_sgf_{case_tag}_profile.png", dpi=150)

    # ---- time series at specified heights ----
    fig, axes = plt.subplots(
        len(time_series_heights),
        1,
        figsize=(10, 2.75 * len(time_series_heights)),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]
    sgf_available_z = [data["profiles"]["z_m"].unique() for _, data in sgf_cases]
    adcp_available_z = adcp_data["profiles"]["z_m"].unique() if adcp_data is not None else None
    if fvcom_data is not None:
        fvcom_avg_z = fvcom_data["profiles"].groupby("sigma_level")["z_m"].mean()

    for ax, height in zip(axes, time_series_heights):
        height_labels = []
        subplot_label = "Column Average"
        if height != "avg":
            first_case_z = sgf_available_z[0]
            representative_z = first_case_z[np.argmin(np.abs(first_case_z - height))]
            height_labels.append(f"SGF {representative_z:.1f} m")
        for index, ((case_name, data), available_z) in enumerate(zip(sgf_cases, sgf_available_z)):
            if height == "avg":
                series = data["profiles"].groupby("time")["speed"].mean().sort_index()
            else:
                nearest_z = available_z[np.argmin(np.abs(available_z - height))]
                series = data["profiles"].loc[
                    data["profiles"]["z_m"] == nearest_z
                ].sort_values("time")
                series = series.set_index("time")["speed"]
            ax.plot(
                _elapsed_hours(series.index, shared_start), series.values, linewidth=3,
                color=sgf_colors[index % len(sgf_colors)], label=case_name,
            )

        if fvcom_data is not None:
            if height == "avg":
                series = fvcom_data["profiles"].groupby("time")["speed"].mean().sort_index()
            else:
                level = int(fvcom_avg_z.index[np.argmin(np.abs(fvcom_avg_z.values - height))])
                series = fvcom_data["profiles"].loc[
                    fvcom_data["profiles"]["sigma_level"] == level
                ].sort_values("time").set_index("time")["speed"]
                height_labels.append(f"FVCOM {fvcom_avg_z.loc[level]:.1f} m")
            ax.plot(
                _elapsed_hours(series.index, shared_start), series.values,
                ":", linewidth=3, color="black", label="FVCOM",
            )

        if adcp_data is not None:
            if height == "avg":
                series = adcp_data["profiles"].groupby("time")["speed"].mean().sort_index()
            else:
                nearest_z = adcp_available_z[np.argmin(np.abs(adcp_available_z - height))]
                series = adcp_data["profiles"].loc[
                    adcp_data["profiles"]["z_m"] == nearest_z
                ].sort_values("time").set_index("time")["speed"]
                height_labels.append(f"ADCP {nearest_z:.1f} m")
            ax.plot(
                _elapsed_hours(series.index, shared_start), series.values,
                "o", markerfacecolor="none", color="black", label="ADCP",
            )

        if height != "avg":
            subplot_label = "; ".join(height_labels)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("speed [m/s]", fontsize=Y_LABEL_FONT_SIZE)
        ax.set_title(subplot_label, fontsize=TITLE_FONT_SIZE)
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONT_SIZE)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.01))
    axes[-1].set_xlabel("hours", fontsize=X_LABEL_FONT_SIZE)
    fig.suptitle(title_range, fontsize=TITLE_FONT_SIZE)
    fig.tight_layout(rect=(0, 0.1, 1, 0.94))
    fig.savefig(OUTPUT_DIR / f"multiple_sgf_{case_tag}_timeseries.png", dpi=150)

    print(
        f"Wrote multiple_sgf_{case_tag}_profile.png and "
        f"multiple_sgf_{case_tag}_timeseries.png to {OUTPUT_DIR}"
    )
    plt.show()


if __name__ == "__main__":
    main()
