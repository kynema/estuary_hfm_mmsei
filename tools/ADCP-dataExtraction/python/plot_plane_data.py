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
case_name = "max_lev3_dmont"
PLANE_DATA_DIR = Path(
    f"/scratch/mkuhn/estuary_flows/milestone/{case_name}/post_processing"
)

# Sampling group (the incflo.post_processing entry) and the label of the
# PlaneSampler within that group.
sampling_group = "plane_sampling_slow"
sampler_label = "ps"

# ---- output times to plot ----
# "last", "all", or a list of indices into the sorted output folders
output_selection = "last"

# ---- field to contour ----
# A particle column name (e.g. "velocityx", "vof") or "speed" for the
# horizontal speed hypot(velocityx, velocityy).
field = "speed"

# ---- which plane of a multi-offset PlaneSampler ----
# Index into the sorted unique positions along the plane normal.
offset_index = 0

# ---- contour appearance ----
n_levels = 40
color_map = "viridis"
color_limits = None  # e.g. (0.0, 2.0) to fix the color scale across times

# ---- optional real-world time ----
# UTC datetime corresponding to simulation time = 0, used only for titles.
start_date = "2024-10-02"
start_time = "22:00:00"

AXIS_NAMES = {"xco": "x", "yco": "y", "zco": "z"}


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
    if name not in df.columns:
        raise ValueError(f"Field '{name}' not found; available: {sorted(df.columns)}")
    return df[name]


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
    contour = ax.contourf(x, y, values, levels=n_levels, cmap=color_map, vmin=vmin, vmax=vmax)
    fig.colorbar(contour, ax=ax, label=field)
    ax.set_xlabel(f"{AXIS_NAMES[first]} [m]")
    ax.set_ylabel(f"{AXIS_NAMES[second]} [m]")
    ax.set_aspect(aspect)

    if start_date is not None and start_time is not None:
        stamp = pd.Timestamp(f"{start_date} {start_time}") + pd.to_timedelta(sim_time_s, unit="s")
        time_label = f"{stamp:%Y-%m-%d %H:%M:%S} (t = {sim_time_s:.1f} s)"
    else:
        time_label = f"t = {sim_time_s:.1f} s"
    ax.set_title(
        f"{sampling_group}.{sampler_label} {field} at {AXIS_NAMES[normal]} = {normal_position:.1f} m\n{time_label}"
    )
    fig.tight_layout()

    figname = f"{sampling_group}_{sampler_label}_{field}_{folder_name}_{case_name}.png"
    fig.savefig(OUTPUT_DIR / figname, dpi=150)
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
