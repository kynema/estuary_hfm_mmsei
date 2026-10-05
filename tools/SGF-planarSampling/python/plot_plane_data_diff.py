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
case_path1 = "max_lev4_turbine_runs/time0"
PLANE_DATA_DIR1 = Path(
    f"/scratch/mkuhn/estuary_flows/milestone/{case_path1}/post_processing"
)
case_path2 = "max_lev4_focused_reruns/time0"
PLANE_DATA_DIR2 = Path(
    f"/scratch/mkuhn/estuary_flows/milestone/{case_path2}/post_processing"
)

case_name = "max_lev4_turbine_diff_time0"

# Sampling group (the incflo.post_processing entry) and the label of the
# PlaneSampler within that group.
sampling_group = "plane_sampling"
sampler_label = "xy"
# sampling_group = "plane_sampling_end"
# sampler_label = "pe"

# ---- output times to plot ----
# "last", "all", or a list of indices into the sorted output folders
output_selection = "all"
stride = 5

# ---- field to contour ----
# A particle column name (e.g. "velocityx", "vof") or "speed" for the
# horizontal speed hypot(velocityx, velocityy).
field = "speed"

# ---- which plane of a multi-offset PlaneSampler ----
# Index into the sorted unique positions along the plane normal.
offset_index = 1

# ---- contour appearance ----
n_levels = 40
color_map = "viridis"
color_limits = (-1.5, 0.5) if field == "speed" or field == "average_speed" else None  # e.g. (0.0, 2.0) to fix the color scale across times

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


def field_values(df1, df2, name):
    mask = np.where(df1["terrain_blank"] < 0.5, 0.0, 1e8)
    if name == "speed":
        return (np.hypot(df1["velocityx"], df1["velocityy"]) - \
               np.hypot(df2["velocityx"], df2["velocityy"])) + mask
    if name == "average_speed":
        return (np.hypot(df1["velocity_mean_avgx"], df1["velocity_mean_avgy"]) - \
               np.hypot(df2["velocity_mean_avgx"], df2["velocity_mean_avgy"])) + mask
    if name not in df1.columns:
        raise ValueError(f"Field '{name}' not found; available: {sorted(df1.columns)}")
    return (df1[name] - df2[name]) + mask


def plot_plane_diff(df1, df2,sim_time_s, folder_name):
    normal, first, second = detect_plane_axes(df1)
    df1, normal_position = select_offset(df1, normal, offset_index)
    df2, _ = select_offset(df2, normal, offset_index)
    df1["field"] = field_values(df1, df2, field)

    # The sampler grid is structured, so pivot it rather than triangulating.
    grid = df1.pivot_table(index=second, columns=first, values="field")
    x = grid.columns.to_numpy()
    y = grid.index.to_numpy()
    values = np.ma.masked_invalid(grid.to_numpy())

    vmin, vmax = color_limits if color_limits is not None else (None, None)
    aspect = "equal" if normal == "zco" else "auto"

    fig, ax = plt.subplots(figsize=(9, 7))
    if color_limits is not None:
        levels = np.linspace(vmin, vmax, n_levels+1)
    contour = ax.contourf(x, y, values, levels=levels if color_limits is not None else n_levels, \
                          cmap=color_map, vmin=vmin, vmax=vmax)
    fig.colorbar(contour, ax=ax, label=field)

    #ax.legend()
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
        f"{field} deficit at {AXIS_NAMES[normal]} = {normal_position:.1f} m {time_label}"
    )
    fig.tight_layout()

    figname = f"{folder_name}_off{offset_index}_{case_name}_{field}.png"
    fig.savefig(OUTPUT_DIR / figname, dpi=150)
    plt.close(fig)
    return figname


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    folders1 = discover_output_folders(PLANE_DATA_DIR1, sampling_group)
    folders2 = discover_output_folders(PLANE_DATA_DIR2, sampling_group)
    if output_selection == "last":
        selected1 = folders1[-1:]
        selected2 = folders2[-1:]
    elif output_selection == "all":
        selected1 = folders1
        selected2 = folders2
    else:
        selected1 = [folders1[i] for i in output_selection]
        selected2 = [folders2[i] for i in output_selection]
    size = min(len(selected1), len(selected2))
    selected1 = selected1[:size]
    selected2 = selected2[:size]

    written = []
    iter = 0
    for folder1, folder2 in zip(selected1, selected2):
        if iter % stride != 0:
            iter += 1
            continue
        df1, sim_time_s1 = load_plane(folder1, sampler_label)
        df2, sim_time_s2 = load_plane(folder2, sampler_label)
        print(f"Loaded {len(df1)} particles from {folder1.name} (t = {sim_time_s1:.1f} s)")
        print(f"Loaded {len(df2)} particles from {folder2.name} (t = {sim_time_s2:.1f} s)")
        written.append(plot_plane_diff(df1, df2, sim_time_s1, folder1.name))
        iter += 1

    print(f"Wrote {', '.join(written)} to {OUTPUT_DIR}")
    plt.show()


if __name__ == "__main__":
    main()
