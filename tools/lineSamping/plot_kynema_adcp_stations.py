#!/usr/bin/env python3

import glob
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.io import loadmat

#--------------------------------------------#
# SETTINGS
FVCOM_TOOLS_DIR = Path(
    "./estuary_hfm_mmsei/tools/FVCOM-dataExtraction"
)
ADCP_DIR = Path("./adcp_data")

ADCP_FILES = {
    "Stablemoor": "STBM_Sig500s_Rosario.mat",
    "Sea Spider": "SeaSpider_Sig250_Rosario.mat",
}

LINE_GROUP = "line_sampling"
LINE_LABELS = {
    "Stablemoor": "STBM",
    "Sea Spider": "SS",
}

# Optional legend labels for Kynema cases.
# Keys are the case-directory names passed through sys.argv.
# Leave empty or omit a case to use its directory name.
CASE_LABELS = {
    # "case_run_coarse_4nodes": "Coarse",
    # "case_run_fine_4nodes": "Fine",
}

SIM_START = pd.Timestamp("2024-10-02 22:00:00")

HAB_M = {
    "Stablemoor": 10.0,
    "Sea Spider": 10.0,
}

# Used only if terrain_blank was not written by LineSampler.
BED_Z = {
    "Stablemoor": None,
    "Sea Spider": None,
}

VOF_WET = 0.5

USE_CACHE = True
REBUILD_CACHE = False
CACHE_DIRNAME = ".kynema_cache"
CACHE_FILENAME = "line_station_timeseries.csv"

OUT_DIR = Path("./kynema_adcp_comparison")
OUT_CSV = OUT_DIR / "kynema_station_timeseries.csv"

SAVE_FIG = True
SHOW_FIG = False
DPI = 150

ZOOM_FIRST_DAY = True
ZOOM_HOURS = 12.0

# Fixed y-axis limits for the Kynema water-depth panel.
# Set to None to use automatic scaling.
WATER_DEPTH_YLIM = (70.0, 100.0)

#--------------------------------------------#

if FVCOM_TOOLS_DIR.exists():
    sys.path.insert(0, str(FVCOM_TOOLS_DIR))

from amrex_particle import AmrexParticleFile

#--------------------------------------------#

def discover_folders(case_dir):
    root = Path(case_dir) / "post_processing"
    folders = [
        Path(f) for f in glob.glob(str(root / f"{LINE_GROUP}*"))
        if Path(f).is_dir()
    ]
    folders = [
        f for f in folders
        if f.name.removeprefix(LINE_GROUP).isdigit()
    ]
    folders.sort(key=lambda f: int(f.name.removeprefix(LINE_GROUP)))
    if not folders:
        raise FileNotFoundError(
            f"No {LINE_GROUP}##### folders found in {root}"
        )
    return folders

#--------------------------------------------#

def load_line_folder(folder):
    with open(folder / "sampling_info.yaml", "r") as fh:
        info = yaml.safe_load(fh)

    labels = {s["index"]: s["label"] for s in info["samplers"]}
    df = AmrexParticleFile(str(folder / "particles"))().copy()
    df["label"] = df["set_id"].map(labels)

    return df, float(info["time"])

#--------------------------------------------#

def surface_z(z, vof):
    z = np.asarray(z, float)
    vof = np.asarray(vof, float)

    order = np.argsort(z)
    z = z[order]
    vof = vof[order]

    crossings = []
    for i in range(len(z) - 1):
        f0 = vof[i] - 0.5
        f1 = vof[i + 1] - 0.5

        if not np.isfinite(f0) or not np.isfinite(f1):
            continue
        if f0 == 0:
            crossings.append(z[i])
        elif f0 * f1 < 0:
            crossings.append(
                z[i] + (0.5 - vof[i]) *
                (z[i + 1] - z[i]) /
                (vof[i + 1] - vof[i])
            )

    return float(max(crossings)) if crossings else np.nan

#--------------------------------------------#

def infer_bed_z(df, station):
    if "terrain_blank" not in df.columns:
        if BED_Z[station] is None:
            raise RuntimeError(
                f"{station}: terrain_blank is not in LineSampler output. "
                f"Set BED_Z['{station}'] near the top of the script."
            )
        return float(BED_Z[station])

    q = df.sort_values("zco")
    z = q["zco"].to_numpy(float)
    blank = q["terrain_blank"].to_numpy(float)

    terrain = np.isfinite(blank) & (blank >= 0.5)
    if not np.any(terrain):
        raise RuntimeError(f"{station}: no terrain_blank == 1 points found.")

    z_terrain = np.nanmax(z[terrain])
    fluid_above = z[(z > z_terrain) & ~terrain]

    if len(fluid_above):
        return 0.5 * (z_terrain + np.nanmin(fluid_above))

    return float(z_terrain)

#--------------------------------------------#

def interpolate_profile(z, q, target):
    z = np.asarray(z, float)
    q = np.asarray(q, float)
    good = np.isfinite(z) & np.isfinite(q)

    if np.sum(good) < 2:
        return np.nan

    z = z[good]
    q = q[good]
    order = np.argsort(z)
    z = z[order]
    q = q[order]

    if target < z[0] or target > z[-1]:
        return np.nan

    return float(np.interp(target, z, q))

#--------------------------------------------#

def process_station(df, station):
    label = LINE_LABELS[station]
    q = df[df["label"] == label].copy()

    if q.empty:
        raise RuntimeError(f"Sampler '{label}' not found for {station}.")

    q = q.sort_values("zco")
    bed = infer_bed_z(q, station)
    surf = surface_z(q["zco"], q["vof"])

    z = q["zco"].to_numpy(float)
    vof = q["vof"].to_numpy(float)
    u = q["velocityx"].to_numpy(float)
    v = q["velocityy"].to_numpy(float)

    wet = (
        np.isfinite(z)
        & np.isfinite(vof)
        & (z > bed)
        & (z < surf)
        & (vof >= VOF_WET)
    )

    if "terrain_blank" in q.columns:
        wet &= q["terrain_blank"].to_numpy(float) < 0.5

    if np.sum(wet) < 2:
        ubar = np.nan
        vbar = np.nan
    else:
        ubar = float(np.nanmean(u[wet]))
        vbar = float(np.nanmean(v[wet]))

    hub_z = bed + HAB_M[station]
    hub_u = interpolate_profile(z[wet], u[wet], hub_z)
    hub_v = interpolate_profile(z[wet], v[wet], hub_z)

    return {
        "station": station,
        "bed_z": bed,
        "surface_z": surf,
        "water_depth": surf - bed,
        "hab_m": HAB_M[station],
        "hub_z": hub_z,
        "hub_u": hub_u,
        "hub_v": hub_v,
        "hub_speed": np.hypot(hub_u, hub_v),
        "depthavg_u": ubar,
        "depthavg_v": vbar,
        "depthavg_speed": np.hypot(ubar, vbar),
    }

#--------------------------------------------#

def process_case(case_dir):
    case_dir = Path(case_dir).resolve()
    case_name = case_dir.name
    folders = discover_folders(case_dir)

    cache_dir = case_dir / CACHE_DIRNAME
    cache_file = cache_dir / CACHE_FILENAME

    cached = pd.DataFrame()

    if USE_CACHE and cache_file.exists() and not REBUILD_CACHE:
        cached = pd.read_csv(cache_file, parse_dates=["datetime"])

        current = {folder.name for folder in folders}
        cached = cached[cached["sampling_folder"].isin(current)].copy()

        # Always use the current case-directory name, even if the folder
        # was renamed after the cache was originally created.
        cached["case"] = case_name

    done = set(cached["sampling_folder"]) if not cached.empty else set()
    todo = [folder for folder in folders if folder.name not in done]

    print(f"{case_name}: {len(folders)} line-sampling folders")
    print(f"  Cached: {len(done)}")
    print(f"  To read: {len(todo)}")

    rows = []

    for i, folder in enumerate(todo, 1):
        df, time_s = load_line_folder(folder)
        print(f"[{i:4d}/{len(todo):4d}] {case_name}/{folder.name}")

        for station in LINE_LABELS:
            row = process_station(df, station)
            row.update({
                "case": case_name,
                "sampling_folder": folder.name,
                "time_s": time_s,
                "datetime": SIM_START + pd.to_timedelta(time_s, unit="s"),
            })
            rows.append(row)

    new_data = pd.DataFrame(rows)

    if cached.empty:
        result = new_data
    elif new_data.empty:
        result = cached
    else:
        result = pd.concat([cached, new_data], ignore_index=True)

    if result.empty:
        return result

    result["datetime"] = pd.to_datetime(result["datetime"])
    result["case"] = case_name
    result = result.drop_duplicates(
        subset=["sampling_folder", "station"],
        keep="last",
    )
    result = result.sort_values(["station", "time_s"])

    if USE_CACHE:
        cache_dir.mkdir(parents=True, exist_ok=True)
        result.to_csv(cache_file, index=False)
        print(f"  Cache: {cache_file}")

    return result

#--------------------------------------------#

def flatten(x, prefix="", out=None):
    if out is None:
        out = {}

    if isinstance(x, dict):
        for k, v in x.items():
            if not str(k).startswith("__"):
                key = f"{prefix}.{k}" if prefix else str(k)
                flatten(v, key, out)
    else:
        try:
            a = np.asarray(x)
            if a.size:
                out[prefix] = a
        except Exception:
            pass

    return out

#--------------------------------------------#

def find_key(flat, names):
    for name in names:
        hits = [
            k for k in flat
            if k.lower().split(".")[-1] == name.lower()
        ]
        if len(hits) == 1:
            return hits[0]
    return None

#--------------------------------------------#

def matlab_time(x):
    x = np.asarray(x).squeeze()

    if x.dtype.kind in "OUS":
        return pd.DatetimeIndex(pd.to_datetime(x))

    x = x.astype(float)
    med = np.nanmedian(x)

    if 700000 < med < 900000:
        return pd.DatetimeIndex(
            pd.to_datetime(x - 719529.0, unit="D", origin="unix")
        )

    if 40000 < med < 100000:
        return pd.DatetimeIndex(
            pd.Timestamp("1858-11-17") + pd.to_timedelta(x, unit="D")
        )

    if 1e9 < med < 3e9:
        return pd.DatetimeIndex(pd.to_datetime(x, unit="s"))

    raise RuntimeError("Could not identify ADCP time format.")

#--------------------------------------------#

def as_time_bins(a, nt):
    a = np.asarray(a, float).squeeze()

    if a.ndim == 1:
        if len(a) != nt:
            raise RuntimeError("ADCP array length does not match time.")
        return a[:, None]

    axes = [i for i, n in enumerate(a.shape) if n == nt]
    if not axes:
        raise RuntimeError(
            f"No dimension of ADCP array {a.shape} matches time length {nt}."
        )

    return np.moveaxis(a, axes[0], 0).reshape(nt, -1)

#--------------------------------------------#

def adcp_vertical_coordinate(z, nt, nbin):
    z = np.asarray(z, float).squeeze()

    if z.ndim == 1:
        if len(z) != nbin:
            raise RuntimeError("ADCP HAB length does not match velocity bins.")
        return z

    if z.size == nbin:
        return z.reshape(-1)

    axes = [i for i, n in enumerate(z.shape) if n == nt]
    if axes:
        zz = np.moveaxis(z, axes[0], 0).reshape(nt, -1)
        if zz.shape[1] == nbin:
            return zz

    raise RuntimeError("Could not align ADCP HAB coordinates.")

#--------------------------------------------#

def interp_adcp_profile(a, coord, target):
    a = np.asarray(a, float)
    out = np.full(a.shape[0], np.nan)

    for i in range(a.shape[0]):
        z = coord if np.ndim(coord) == 1 else coord[i]
        out[i] = interpolate_profile(z, a[i], target)

    return out

#--------------------------------------------#

def load_adcp(station):
    path = ADCP_DIR / ADCP_FILES[station]
    if not path.exists():
        raise FileNotFoundError(f"ADCP file not found: {path}")

    d = flatten(loadmat(path, simplify_cells=True))

    kt = find_key(d, ["t", "time", "datetime", "datenum", "mtime"])
    ku = find_key(d, ["u", "east", "east_vel", "east_velocity", "vel_east"])
    kv = find_key(d, ["v", "north", "north_vel", "north_velocity", "vel_north"])
    kz = find_key(d, ["z", "hab", "height_above_bed"])

    if None in (kt, ku, kv, kz):
        raise RuntimeError(
            f"{station}: could not identify ADCP time/u/v/HAB variables."
        )

    time = matlab_time(d[kt])
    u = as_time_bins(d[ku], len(time))
    v = as_time_bins(d[kv], len(time))
    coord = adcp_vertical_coordinate(d[kz], len(time), u.shape[1])

    valid = (
        np.mean(np.isfinite(u), axis=0) >= 0.90
    ) & (
        np.mean(np.isfinite(v), axis=0) >= 0.90
    )

    if not np.any(valid):
        valid = np.any(np.isfinite(u) & np.isfinite(v), axis=0)

    ubar = np.nanmean(u[:, valid], axis=1)
    vbar = np.nanmean(v[:, valid], axis=1)

    hub_u = interp_adcp_profile(u, coord, HAB_M[station])
    hub_v = interp_adcp_profile(v, coord, HAB_M[station])

    return pd.DataFrame({
        "datetime": time,
        "depthavg_speed": np.hypot(ubar, vbar),
        "hub_speed": np.hypot(hub_u, hub_v),
    })

#--------------------------------------------#

def make_station_plot(model, obs, station, zoom_hours=None):
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    for case, q in model[model["station"] == station].groupby("case", sort=False):
        case_label = CASE_LABELS.get(case, case)
        axs[0].plot(q["datetime"], q["depthavg_speed"], label=f"{case_label} Kynema")
        axs[1].plot(q["datetime"], q["hub_speed"], label=f"{case_label} Kynema")
        axs[2].plot(q["datetime"], q["water_depth"], label=f"{case_label} Kynema")

    axs[0].plot(obs["datetime"], obs["depthavg_speed"], "--", label="ADCP")
    axs[1].plot(obs["datetime"], obs["hub_speed"], "--", label="ADCP")

    axs[0].set_ylabel("Depth-avg speed [m/s]")
    axs[1].set_ylabel(f"Speed at {HAB_M[station]:g} m HAB [m/s]")
    axs[2].set_ylabel("Water depth [m]")
    axs[2].set_xlabel("Time [UTC]")

    if WATER_DEPTH_YLIM is not None:
        axs[2].set_ylim(WATER_DEPTH_YLIM)

    if zoom_hours is not None:
        zoom_start = SIM_START
        zoom_end = SIM_START + pd.Timedelta(hours=zoom_hours)
        axs[2].set_xlim(zoom_start, zoom_end)
        title = f"{station} - first {zoom_hours:g} hours"
        suffix = f"_first_{zoom_hours:g}h"
    else:
        title = station
        suffix = ""

    for ax in axs:
        ax.grid(True, alpha=0.25)
        ax.legend()

    fig.suptitle(title)
    fig.tight_layout()

    outfile = OUT_DIR / (
        f"{station.lower().replace(' ', '_')}"
        f"_kynema_vs_adcp{suffix}.png"
    )

    if SAVE_FIG:
        fig.savefig(outfile, dpi=DPI, bbox_inches="tight")
    if SHOW_FIG:
        plt.show()

    plt.close(fig)

#--------------------------------------------#

def plot_station(model, station):
    obs = load_adcp(station)

    if ZOOM_FIRST_DAY:
        make_station_plot(
            model,
            obs,
            station,
            zoom_hours=ZOOM_HOURS,
        )
    else:
        make_station_plot(
            model,
            obs,
            station,
        )

#--------------------------------------------#

def main():
    if len(sys.argv) < 2:
        print(f"Usage: python {Path(sys.argv[0]).name} case1 [case2 ...]")
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    frames = [process_case(case) for case in sys.argv[1:]]
    model = pd.concat(frames, ignore_index=True)
    model = model.sort_values(["case", "station", "time_s"])
    model.to_csv(OUT_CSV, index=False)

    for station in LINE_LABELS:
        plot_station(model, station)

    print(f"Wrote {OUT_CSV}")
    print(f"Wrote figures to {OUT_DIR}")

#--------------------------------------------#

if __name__ == "__main__":
    main()

