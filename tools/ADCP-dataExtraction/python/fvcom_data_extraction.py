"""Load FVCOM velocity data at the Rosario ADCP deployment locations.

The loader returns the same common data-dictionary shape used by the ADCP
and SGF extraction modules.  FVCOM sigma-layer heights are converted to
meters above the local seafloor, allowing direct comparison in plots.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# inflowPrepMMC (lat/lon -> UTM) lives with the FVCOM extraction tools.
FVCOM_TOOLS_DIR = Path(__file__).parent.parent.parent / "FVCOM-dataExtraction"
if FVCOM_TOOLS_DIR.exists():
    sys.path.insert(0, str(FVCOM_TOOLS_DIR))

import inflowPrepMMC as mmc


# ---- deployment locations (lat, lon), matching the FVCOM notebook ----
DEPLOYMENT_LATLON = {
    "stbm": (48.56260, -122.76690),
    "ss": (48.56270, -122.76510),
}
UTM_ZONE = "10 U"
MJD_EPOCH = pd.Timestamp("1858-11-17")
FVCOM_TO_UTC_HOURS = 8.0


def datetime_to_mjd(date, time=None):
    """Convert a date/time value to Modified Julian Day.

    ``date`` may be a date string, datetime-like value, or ``Timestamp``.
    When provided, ``time`` is combined with ``date`` before conversion.
    """
    timestamp = pd.Timestamp(f"{date} {time}") if time is not None else pd.Timestamp(date)
    return (timestamp - MJD_EPOCH).total_seconds() / 86400.0


def utc_datetime_to_fvcom_mjd(date, time=None):
    """Convert a UTC date/time to the raw MJD convention used by FVCOM.

    FVCOM raw MJD timestamps are ``FVCOM_TO_UTC_HOURS`` behind UTC, matching
    the offset used by the netCDF comparison workflow.
    """
    timestamp = pd.Timestamp(f"{date} {time}") if time is not None else pd.Timestamp(date)
    return datetime_to_mjd(timestamp - pd.Timedelta(hours=FVCOM_TO_UTC_HOURS))


def resolve_fvcom_window(
    *,
    start_mjd=None,
    end_mjd=None,
    duration_hours=None,
    start_date=None,
    start_time=None,
    end_date=None,
    end_time=None,
):
    """Resolve an exclusive native-MJD or UTC FVCOM time window.

    Return ``(raw_start_mjd, duration_hours, utc_start)``.  UTC windows
    require complete start and end date/time values and are converted to the
    raw FVCOM MJD convention. MJD windows require ``start_mjd`` and exactly
    one of ``end_mjd`` or ``duration_hours``. ``utc_start`` is ``None`` for
    MJD windows.
    """
    has_mjd = start_mjd is not None or end_mjd is not None
    utc_values = (start_date, start_time, end_date, end_time)
    has_utc = any(value is not None for value in utc_values)

    if has_mjd and has_utc:
        raise ValueError(
            "Specify the FVCOM window with either MJD values or UTC "
            "start_date/start_time/end_date/end_time, not both."
        )

    if has_utc:
        if any(value is None for value in utc_values):
            raise ValueError(
                "UTC FVCOM windows require start_date, start_time, end_date, "
                "and end_time."
            )
        utc_start = pd.Timestamp(f"{start_date} {start_time}")
        utc_end = pd.Timestamp(f"{end_date} {end_time}")
        duration = (utc_end - utc_start).total_seconds() / 3600.0
        if duration <= 0:
            raise ValueError("The UTC FVCOM end timestamp must be after the start timestamp.")
        return utc_datetime_to_fvcom_mjd(utc_start), duration, utc_start

    if start_mjd is None:
        raise ValueError(
            "Set start_mjd for an MJD window, or provide a complete UTC window."
        )
    if end_mjd is not None and duration_hours is not None:
        raise ValueError("Set either end_mjd or duration_hours for an MJD window, not both.")
    if end_mjd is not None:
        duration = (end_mjd - start_mjd) * 24.0
    elif duration_hours is not None:
        duration = duration_hours
    else:
        raise ValueError("An MJD FVCOM window requires end_mjd or duration_hours.")
    if duration <= 0:
        raise ValueError("The MJD FVCOM end timestamp must be after the start timestamp.")
    return start_mjd, duration, None


def load_fvcom_column(unit, reduced_dir, start_mjd, duration_hours,
                      start_date=None, start_time=None):
    """Extract the FVCOM column nearest a deployment location.

    When ``start_date`` and ``start_time`` are provided, they define the
    real-world timestamp corresponding to the beginning of the selected
    FVCOM window. Otherwise, the returned ``time`` column contains elapsed
    seconds and ``title_range`` reports native FVCOM UTC timestamps.
    """
    unit = unit.strip().lower()
    if unit not in DEPLOYMENT_LATLON:
        raise ValueError(f"Unknown fvcom_unit '{unit}'; expected 'ss' or 'stbm'.")

    reduced_dir = Path(reduced_dir)
    time_mjd = np.load(reduced_dir / "PS_time.npy")
    x = np.load(reduced_dir / "PS_x.npy")
    y = np.load(reduced_dir / "PS_y.npy")
    h_bathymetry = -np.load(reduced_dir / "PS_h.npy")
    nv = np.load(reduced_dir / "PS_nv.npy")
    siglay = np.load(reduced_dir / "PS_siglay.npy")
    u = np.load(reduced_dir / "PS_u.npy")
    v = np.load(reduced_dir / "PS_v.npy")
    zeta_node = np.load(reduced_dir / "PS_zeta.npy")

    # Water-surface elevation at cell centers, averaged from the nodes.
    zeta = np.empty((len(time_mjd), len(x)))
    for t in range(zeta.shape[0]):
        zeta[t, :] = np.mean(zeta_node[t, nv], axis=0)

    lat, lon = DEPLOYMENT_LATLON[unit]
    coords = mmc.inflowPrepMMC()
    utm_x, utm_y, utm_zone = coords.LatLonToUTM(lat, lon)
    if utm_zone != UTM_ZONE:
        print(f"Warning: UTM zone of {unit.upper()} ({utm_zone}) does not match {UTM_ZONE}.")

    column = int(np.argmin(np.hypot(x - utm_x, y - utm_y)))
    distance = np.hypot(x[column] - utm_x, y[column] - utm_y)
    print(f"Nearest FVCOM column: {column}  distance to deployment location: {distance:.1f} m")

    end_mjd = start_mjd + duration_hours / 24.0
    i0 = int(np.argmin(np.abs(time_mjd - start_mjd)))
    i1 = int(np.argmin(np.abs(time_mjd - end_mjd))) + 1
    window = slice(i0, i1 + 1)

    sigma = siglay[:, column]
    seafloor = h_bathymetry[column]
    surface = zeta[window, column]
    z = seafloor - (surface[:, None] - seafloor) * sigma[None, :]
    z_above_seafloor = z - seafloor

    u_col = u[window, :, column]
    v_col = v[window, :, column]
    n_time, n_sigma = u_col.shape

    elapsed_s = (time_mjd[window] - time_mjd[i0]) * 24.0 * 3600.0
    if start_date is not None and start_time is not None:
        times = pd.Timestamp(f"{start_date} {start_time}") + pd.to_timedelta(elapsed_s, unit="s")
    else:
        times = pd.Series(elapsed_s)

    profiles = pd.DataFrame(
        {
            "time": np.repeat(np.asarray(times), n_sigma),
            "elapsed_s": np.repeat(elapsed_s, n_sigma),
            "sigma_level": np.tile(np.arange(n_sigma), n_time),
            "z_m": z_above_seafloor.ravel(),
            "u": u_col.ravel(),
            "v": v_col.ravel(),
        }
    )
    profiles["speed"] = np.hypot(profiles["u"], profiles["v"])

    scalars = pd.DataFrame({"time": np.asarray(times), "elapsed_s": elapsed_s})
    scalars["mjd"] = time_mjd[window]

    if start_date is not None and start_time is not None:
        title_range = (
            f"{scalars['time'].min():%Y-%m-%d %H:%M:%S} to "
            f"{scalars['time'].max():%Y-%m-%d %H:%M:%S}"
        )
    else:
        start_utc = MJD_EPOCH + pd.to_timedelta(time_mjd[i0], unit="D")
        stop_utc = MJD_EPOCH + pd.to_timedelta(time_mjd[i1], unit="D")
        title_range = f"{start_utc:%Y-%m-%d %H:%M:%S} to {stop_utc:%Y-%m-%d %H:%M:%S} (FVCOM UTC)"

    label = unit.upper()
    return {
        "unit": label,
        "file_prefix": f"FVCOM_{label}",
        "plot_label": f"FVCOM {label}",
        "time": np.asarray(times),
        "scalars": scalars,
        "profiles": profiles,
        "waves": None,
        "title_range": title_range,
    }