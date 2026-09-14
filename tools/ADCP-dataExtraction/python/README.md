# ADCP Data Extraction — Python Post-Processing

This directory contains Python ports of a subset of the Matlab analysis scripts in `../Matlab/`, for producing plots from the "clean" summary `.mat` files that those Matlab scripts generate. Run with the `kynema-env` conda environment.

These scripts output figures that are written to `figures/`.

## adcp_data_extraction.py

Reusable data-loading module, not intended to be run directly. Loads and
filters either of the two Rosario ADCP datasets and returns a common dict
shape (`unit`, `file_prefix`, `plot_label`, `time`, `scalars`, `profiles`,
`waves`, `title_range`) so downstream scripts don't need to special-case the
underlying `.mat` file layout.

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

## plot_adcp_data.py

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

## User-editable inputs

Near the top of `plot_adcp_data.py`:

- `adcp_unit` — `"ss"`/`"SS"` for Sea Spider or `"stbm"`/`"STBM"` for
  Stablemoor.
- `start_date`, `stop_date`, `start_time` — narrow the analysis window. Leave
  as `None` to use the full data record — in that case, plot titles show the
  actual min/max timestamps found in the data rather than the literal
  strings "start"/"end".
- Turbine/site params: `R` (rotor radius), `hubdepth`, `waterdepth`,
  `near_bottom_height`, and `time_series_heights`.

