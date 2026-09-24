# ADCP, SGF, and FVCOM Data Extraction — Python Post-Processing

This directory contains a collection of Python scripts for post-processing and plotting ADCP, Kynema-SGF, and FVCOM data. The scripts are designed to be run in the `kynema-env` conda environment (see `estuary_hfm_mmsei/environment.yml`).

These scripts output figures that are written to `figures/`.

## FVCOM MJD and UTC Time Windows

The reduced FVCOM data `PS_time.npy` values use raw Modified Julian Days (MJD),
measured in days from 1858-11-17 00:00:00 UTC. The FVCOM comparison workflow
uses an additional 8-hour FVCOM-to-UTC offset: a raw FVCOM timestamp is 8
hours behind the UTC time used for SGF/FVCOM/ADCP comparison plots. This is
handled by `fvcom_data_extraction.utc_datetime_to_fvcom_mjd()`.

## ADCP data extraction 

There are two scripts associated with ADCP data extraction and plotting:
- `adcp_data_extraction.py` — reusable data-loading module, not intended to be run directly. Loads and filters either of the two Rosario ADCP datasets and returns a common dict shape (`unit`, `file_prefix`, `plot_label`, `time`, `scalars`, `profiles`, `waves`, `title_range`) so downstream scripts don't need to special-case the underlying `.mat` file layout.
- `plot_adcp_data.py` — plotting script that calls into `adcp_data_extraction.py`. Produces average speed profile vs. height above seafloor, time series of speed at specified heights, and (for STBM only) surface wave time series and wave height vs. peak period scatter. User-editable inputs are near the top of the script.

#### `adcp_data_extraction.py`

- `load_ss_data(...)` — loads `../Matlab/data/SeaSpider_Sig250_Rosario.mat`
  (Sea Spider seafloor tripod, uplooking Sig250). `profiles` speed is the
  **horizontal-only** magnitude, $\sqrt{u^2+v^2}$, matching
  `../Matlab/RosarioSeaSpiderAnalysis.m`'s `speed = (u.^2+v.^2).^.5`.
- `load_stbm_data(...)` — loads `../Matlab/data/STBM_Sig500s_Rosario.mat`
  (STBM mid-water mooring, up- and down-looking Sig500s), using the
  fixed-height-above-seabed gridded data (`fixedz_profiles`). `profiles`
  speed is the **3D** magnitude, $\sqrt{u^2+v^2+w^2}$, matching
  `../Matlab/RosarioStablemoorAnalysis.m`'s `speed = (u.^2+v.^2+w.^2).^.5`
  (a conservative correction for mooring tilt/trim; per the Matlab script's
  comment this does not change the AEP result vs. horizontal-only speed).
  Also returns a `waves` DataFrame (or `None`) from the `Waves` struct array
  in the `.mat` file, which is populated from an externally preprocessed
  `STBM_Sig500up_preprocessed.mat` file (see the "surface waves" section of
  the Matlab script — this preprocessing step requires
  `sigProcess_all.m`/`UVZwaves.m` from the NortekCodes/SWIFTcodes repos,
  which are not part of this workspace).
- `load_adcp_data(adcp_unit, ...)` — dispatches to `load_ss_data` or
  `load_stbm_data` based on a case-insensitive `adcp_unit` string (`"ss"` or
  `"stbm"`).
- `hub_height_above_seafloor(data, hubdepth, waterdepth)` — computes hub
  height above the seafloor, using the measured `waterdepth` scalar for SS or
  the given fixed `waterdepth` for STBM.
- `reflection_zone()` — returns the fixed `(start, end)` heights above
  seafloor (83–94 m) of the surface-reflection zone shaded and excluded in
  the profile plots, matching the `area([...],[94 94],83)` shading in the
  Matlab scripts.

#### `plot_adcp_data.py`

Plotting script that calls into `adcp_data_extraction.py`. Produces:

1. **Average speed profile** vs. height above seafloor (z=0 at seafloor),
   matching the first panel of `SeaSpider_Sig250_Rosario_profiles.png` /
   `STBM_Sig500s_Rosario_profiles.png`. The fixed surface-reflection zone
   from `reflection_zone()` is shaded red and excluded from the plotted
   profile line.
2. **Time series of speed** at each height in `time_series_heights` (default:
   hub depth, near-bottom height, and column average).
3. For **STBM only**, if wave data is present: **surface wave time series**
   (significant wave height, peak wave period), matching the `fh21` plot in
   `../Matlab/RosarioStablemoorAnalysis.m`, plus a **wave height vs. peak
   period** scatter, matching the `fh22` plot. If `Waves` is not present in
   the `.mat` file (or has no data in the requested range), these plots are
   skipped with a printed message.

### User-editable inputs

Near the top of `plot_adcp_data.py`:

- `adcp_unit` — `"ss"`/`"SS"` for Sea Spider or `"stbm"`/`"STBM"` for
  Stablemoor.
- `start_date`, `stop_date`, `start_time` — narrow the analysis window. Leave
  as `None` to use the full data record — in that case, plot titles show the
  actual min/max timestamps found in the data rather than the literal
  strings "start"/"end".
- Turbine/site params: `R` (rotor radius), `hubdepth`, `waterdepth`,
  `near_bottom_height`, and `time_series_heights`.

## SGF data extraction
#### `sgf_data_extraction.py`

Reusable data-loading module, not intended to be run directly. Loads
Kynema-SGF simulation line-sampling output (particle data written to
`line_sampling#####` folders under a simulation's `post_processing`
directory, one folder per output time) and returns the same common dict
shape used by `adcp_data_extraction.py` (`unit`, `file_prefix`,
`plot_label`, `time`, `scalars`, `profiles`, `waves`, `title_range`), so
`plot_sgf_data.py` can reuse the same plotting logic as `plot_adcp_data.py`.

- Each `line_sampling#####` folder contains a `sampling_info.yaml` (simulation
  time in seconds, plus the list of `LineSampler`s and their labels, e.g.
  `SS`/`STBM`) and a `particles/` subfolder with AMReX native particle binary
  data, read via `AmrexParticleFile` from
  `../../FVCOM-dataExtraction/amrex_particle.py`.
- Kynema-SGF defines `z=0` at the water surface (sampled `zco` is negative
  below the surface), whereas the ADCP data defines `z=0` at the seafloor. To
  align the two datasets, `z_m` (height above seafloor) is computed as
  `zco + waterdepth`, where `waterdepth` is the local still-water depth:
  STBM uses the fixed `adcp_data_extraction.STBM_WATERDEPTH` constant (the
  same surveyed depth used by `plot_adcp_data.py`); SS uses the mean measured
  `waterdepth` scalar from `adcp_data_extraction.load_ss_data()`.
- `velocityx`/`velocityy`/`velocityz` are renamed to `u`/`v`/`w`, and `speed`
  is the horizontal-only magnitude, $\sqrt{u^2+v^2}$.
- Particles with `vof` below `vof_threshold` (default `0.5`, i.e. mostly air,
  above the free surface) are excluded from `profiles`, since velocity there
  is not physically meaningful for comparison to the ADCP data. Pass
  `vof_threshold=None` to disable this filtering.
- `load_ss_sgf_data(root_dir, ...)` / `load_stbm_sgf_data(root_dir, ...)` —
  load the SS or STBM `LineSampler` data across all output times found under
  `root_dir` (a simulation's `post_processing` directory).
- `load_sgf_data(sgf_unit, root_dir, ...)` — dispatches to the above based on
  a case-insensitive `sgf_unit` string (`"ss"` or `"stbm"`).
- Optional `start_date`/`start_time` arguments (e.g. `"2024-10-22"`/
  `"22:00:00"`) give the real-world UTC datetime corresponding to simulation
  time=0, converting the `time` column from raw simulation seconds to real
  datetimes comparable to the ADCP timestamps. If omitted, `time` stays as
  raw simulation seconds.

#### `plot_sgf_data.py`

Plotting script that calls into `sgf_data_extraction.py`, mirroring
`plot_adcp_data.py`. Produces:

1. **Average speed profile** vs. height above seafloor (z=0 at seafloor).
2. **Time series of speed** at each height in `time_series_heights` (default:
   30 m, 60 m, and column average).

### User-editable inputs

Near the top of `plot_sgf_data.py`:

- `sgf_unit` — `"ss"`/`"SS"` for Sea Spider or `"stbm"`/`"STBM"` for
  Stablemoor.
- `SGF_DATA_DIR` — the simulation's `post_processing` directory containing
  the `line_sampling#####` folders. This location is expected to move in the
  future; update by hand.
- `start_date`, `start_time` — the real-world UTC datetime corresponding to
  simulation time=0. Leave as `None` to keep the `time` column as raw
  simulation seconds.
- `time_series_heights` — heights above seafloor (m) to plot in the time
  series, plus optionally `"avg"` for the column-averaged speed.

#### `plot_multiple_sgf_data.py`

Compare the cases in `case_names` using a shared UTC window, with optional
FVCOM and ADCP overlays. Cases are read from
`SGF_CASES_DIR / case_name / "post_processing"`.

## FVCOM data extraction

- `fvcom_data_extraction.py` loads reduced FVCOM columns and provides the
  MJD/UTC conversion helpers.
- `plot_fvcom_data.py` plots one FVCOM column, optionally with ADCP data.
