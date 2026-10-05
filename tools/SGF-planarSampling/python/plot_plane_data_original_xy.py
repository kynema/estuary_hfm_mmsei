"""
Contour-plot Kynema-SGF PlaneSampler output (AMReX particles arranged in a
plane), the plane-sampling counterpart to plot_sgf_data.py.

Kynema-SGF writes one folder per output time under the post_processing
directory, named <sampling group><#####> (e.g. plane_sampling_slow00120 for
the "plane_sampling_slow" group in milestone/rosario.inp). Each folder holds:
  - sampling_info.yaml: simulation time (s) and the list of samplers
    (index/label/type)
  - particles/: AMReX particle data read via AmrexParticleFile from
    ../../FVCOM-dataExtraction/amrex_particle.py

A PlaneSampler may contain several parallel planes (the "offsets" input), so
the plane normal is detected from the sampled coordinates and one offset is
selected for plotting.

Edit the user settings below to select the sampling group, sampler label,
output time(s), field, and offset.
"""

import glob
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# AmrexParticleFile lives with the FVCOM extraction tools.
FVCOM_TOOLS_DIR = Path(__file__).parent.parent.parent / "FVCOM-dataExtraction"
if FVCOM_TOOLS_DIR.exists():
    sys.path.insert(0, str(FVCOM_TOOLS_DIR))

from amrex_particle import AmrexParticleFile

OUTPUT_DIR = Path(__file__).parent / "figures"

# ---- data location ----
case_name = "max_lev1"
PLANE_DATA_DIR = Path(
    f"/scratch/mkuhn/estuary_flows/milestone/{case_name}/post_processing"
)

# Sampling group (the incflo.post_processing entry) and the label of the
# PlaneSampler within that group.
sampling_group = "plane_sampling_slow"
sampler_label = "ps"
# sampling_group = "plane_sampling_end"
# sampler_label = "pe"

# ---- output times to plot ----
# "last", "all", or a list of indices into the sorted output folders
output_selection = "all"

# ---- field to contour ----
# A particle column name (e.g. "velocityx", "vof") or "speed" for the
# horizontal speed hypot(velocityx, velocityy).
field = "tke"

# ---- which plane of a multi-offset PlaneSampler ----
# Index into the sorted unique positions along the plane normal.
offset_index = 0

# ---- contour appearance ----
n_levels = 40
color_map = "viridis"
color_limits = (0.0, 0.5) if field == "perturb_speed" else ((0.0, 0.3) if field == "tke" else (0.0, 2.0))  # e.g. (0.0, 2.0) to fix the color scale across times

# ---- optional real-world time ----
# Shared native SGF/FVCOM UTC datetime corresponding to simulation time = 0,
# used only for titles. This is the simulation's own native calendar (e.g.
# matching the FVCOM initialization used to start the Kynema-SGF run), not
# the real-world ADCP measurement calendar -- see
# adcp_data_extraction.map_adcp_to_timeline() and plot_sgf_data.py's
# start_date/adcp_start_date for how the two calendars are related.
start_date = "2015-04-07"
start_time = "01:00:00"

# Real-world ADCP measurement UTC datetime that start_date/start_time (sim
# time = 0) corresponds to, used to also show the matching ADCP calendar date.
adcp_start_date = "2024-10-04"
adcp_start_time = "12:00:00"

# ---- deployment-location velocity arrows ----
plot_deployment_arrows = True
arrow_length_scale = 500.  # matplotlib quiver "scale"; None lets matplotlib auto-scale
arrow_max_distance = None  # skip a location if no sampled point is within this distance (m); None disables the check

AXIS_NAMES = {"xco": "x", "yco": "y", "zco": "z"}

STBM_loc = [2798.3098930663546, 2510.166347393766]
SS_loc = [2931.081119452836, 2521.688631718047]


def discover_output_folders(root_dir, group):
    """Return the sampling output folders for a group, sorted by the numeric
    suffix in their name (which is not zero-padded to a fixed width)."""
    folders = [Path(f) for f in glob.glob(str(Path(root_dir) / f"{group}*")) if Path(f).is_dir()]
    folders = [f for f in folders if f.name.removeprefix(group).isdigit()]
    folders = sorted(folders, key=lambda f: int(f.name.removeprefix(group)))
    if not folders:
        raise FileNotFoundError(f"No {group}##### folders found in {root_dir}")
    return folders


def load_plane(folder, label):
    """Load one output folder and return (DataFrame for label, sim time [s])."""
    with open(folder / "sampling_info.yaml", "r") as fh:
        info = yaml.safe_load(fh)

    sim_time_s = float(info["time"])
    set_id_to_label = {s["index"]: s["label"] for s in info["samplers"]}
    if label not in set_id_to_label.values():
        raise ValueError(f"Sampler '{label}' not in {folder}; available: {sorted(set_id_to_label.values())}")

    df = AmrexParticleFile(str(folder / "particles"))().copy()
    df["label"] = df["set_id"].map(set_id_to_label)
    df = df[df["label"] == label].copy()
    if df.empty:
        raise ValueError(f"No particles for sampler '{label}' in {folder}.")
    return df, sim_time_s


def detect_plane_axes(df):
    """Return (normal, first, second) coordinate column names for an
    axis-aligned sampling plane, identifying the normal as the coordinate with
    the fewest distinct values (one per offset)."""
    counts = {col: len(np.unique(df[col].to_numpy())) for col in ("xco", "yco", "zco")}
    normal = min(counts, key=counts.get)
    in_plane = [col for col in ("xco", "yco", "zco") if col != normal]
    if counts[normal] >= min(counts[col] for col in in_plane):
        raise ValueError(
            "Could not identify an axis-aligned plane normal from the sampled "
            f"coordinates (distinct values per axis: {counts})."
        )
    return normal, in_plane[0], in_plane[1]


def select_offset(df, normal, index):
    """Restrict the data to one of the parallel planes of a PlaneSampler."""
    positions = np.unique(df[normal].to_numpy())
    if not -len(positions) <= index < len(positions):
        raise IndexError(f"offset_index {index} out of range; {len(positions)} plane(s) available.")
    position = positions[index]
    return df[df[normal] == position].copy(), position


def field_values(df, name):
    if name == "speed":
        return np.hypot(df["velocityx"], df["velocityy"])
    if name == "average_speed":
        return np.hypot(df["velocity_mean_avgx"], df["velocity_mean_avgy"])
    if name == "perturb_speed":
        return np.hypot(df["velocityx"]-df["velocity_mean_avgx"], df["velocityy"]-df["velocity_mean_avgy"])
    if name == "tke":
        return (0.5 * (np.power(df["velocityx"]-df["velocity_mean_avgx"],2) + np.power(df["velocityy"]-df["velocity_mean_avgy"],2) + np.power(df["velocityz"]-df["velocity_mean_avgz"],2)))
    if name not in df.columns:
        raise ValueError(f"Field '{name}' not found; available: {sorted(df.columns)}")
    return df[name]


def nearest_deployment_points(df, first, second):
    """For each deployment location, find the sampled point (row of df)
    nearest in the plane's in-plane coordinates. Returns a dict of
    label -> (row, distance).
    """
    locations = {"STBM": STBM_loc, "SS": SS_loc}
    points = df[[first, second]].to_numpy()
    nearest = {}
    for label, location in locations.items():
        distances = np.hypot(points[:, 0] - location[0], points[:, 1] - location[1])
        i = int(np.argmin(distances))
        if arrow_max_distance is not None and distances[i] > arrow_max_distance:
            print(f"Skipping {label} arrow: nearest sampled point is {distances[i]:.1f} m away (> arrow_max_distance).")
            continue
        nearest[label] = (df.iloc[i], distances[i])
    return nearest


def plot_plane(df, sim_time_s, folder_name):
    normal, first, second = detect_plane_axes(df)
    df, normal_position = select_offset(df, normal, offset_index)
    df["field"] = field_values(df, field)

    # The sampler grid is structured, so pivot it rather than triangulating.
    grid = df.pivot_table(index=second, columns=first, values="field")
    x = grid.columns.to_numpy()
    y = grid.index.to_numpy()
    values = np.ma.masked_invalid(grid.to_numpy())

    vmin, vmax = color_limits if color_limits is not None else (None, None)
    aspect = "equal" if normal == "zco" else "auto"

    fig, ax = plt.subplots(figsize=(9, 7))
    levels = np.linspace(vmin, vmax, n_levels+1)
    contour = ax.contourf(x, y, values, levels=levels, cmap=color_map, vmin=vmin, vmax=vmax)
    fig.colorbar(contour, ax=ax, label=field)
    ax.plot(STBM_loc[0], STBM_loc[1], "ro", label="STBM")
    ax.plot(SS_loc[0], SS_loc[1], "bo", label="SS")
    # ax.vlines(STBM_loc[0], ymin=y.min(), ymax=y.max(), colors="r", linestyles="--", label="STBM, y+" + str(normal_position-STBM_loc[1]) + "m")
    # ax.vlines(SS_loc[0], ymin=y.min(), ymax=y.max(), colors="b", linestyles="--", label="SS, y+" + str(normal_position-SS_loc[1]) + "m")

    if plot_deployment_arrows:
        # velocityx/y/z are in AMR-Wind's fixed x/y/z frame, so project onto
        # whichever two components correspond to the plane's in-plane axes.
        vel_component = {"xco": "velocityx", "yco": "velocityy", "zco": "velocityz"}
        avg_vel_component = {"xco": "velocity_mean_avgx", "yco": "velocity_mean_avgy", "zco": "velocity_mean_avgz"}
        i = 0
        for label, (row, distance) in nearest_deployment_points(df, first, second).items():
            uavg = row[avg_vel_component[first]]
            vavg = row[avg_vel_component[second]]
            u = row[vel_component[first]]
            v = row[vel_component[second]]
            if (field == "perturb_speed"):
                u -= uavg
                v -= vavg
            if (i == 0):
                l_str = "perturb velocity" if (field == "perturb_speed") else "velocity"
                avgl_str = "avg velocity"
            else:
                l_str = ""
                avgl_str = ""
            i += 1
            # scale=1/arrow_length_scale with scale_units="xy" makes plotted arrow
            # length (in data units) equal speed * arrow_length_scale.
            ax.quiver(
                row[first], row[second], uavg, vavg,
                color="w", scale=1.0 / arrow_length_scale, scale_units="xy",
                angles="xy", zorder=5, label=avgl_str
            )
            ax.quiver(
                row[first], row[second], u, v,
                color="k", scale=1.0 / arrow_length_scale, scale_units="xy",
                angles="xy", zorder=5, label=l_str
            )

    ax.legend()
    ax.set_xlabel(f"{AXIS_NAMES[first]} [m]")
    ax.set_ylabel(f"{AXIS_NAMES[second]} [m]")
    ax.set_aspect(aspect)

    if start_date is not None and start_time is not None:
        stamp = pd.Timestamp(f"{start_date} {start_time}") + pd.to_timedelta(sim_time_s, unit="s")
        time_label = f"(t = {sim_time_s:.1f} s)\nFVCOM: {stamp:%Y-%m-%d %H:%M:%S}"
        if adcp_start_date is not None and adcp_start_time is not None:
            adcp_stamp = pd.Timestamp(f"{adcp_start_date} {adcp_start_time}") + pd.to_timedelta(sim_time_s, unit="s")
            time_label += f"\nADCP: {adcp_stamp:%Y-%m-%d %H:%M:%S}"
    else:
        time_label = f"t = {sim_time_s:.1f} s"
    ax.set_title(
        f"{field} at {AXIS_NAMES[normal]} = {normal_position:.1f} m {time_label}"
    )
    fig.tight_layout()

    figname = f"{folder_name}_off{offset_index}_{case_name}_{field}.png"
    fig.savefig(OUTPUT_DIR / figname, dpi=150)
    plt.close(fig)
    return figname


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    folders = discover_output_folders(PLANE_DATA_DIR, sampling_group)
    if output_selection == "last":
        selected = folders[-1:]
    elif output_selection == "all":
        selected = folders
    else:
        selected = [folders[i] for i in output_selection]

    written = []
    for folder in selected:
        df, sim_time_s = load_plane(folder, sampler_label)
        print(f"Loaded {len(df)} particles from {folder.name} (t = {sim_time_s:.1f} s)")
        written.append(plot_plane(df, sim_time_s, folder.name))

    print(f"Wrote {', '.join(written)} to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()
