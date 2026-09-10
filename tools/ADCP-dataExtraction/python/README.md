# ADCP Data Extraction — Python Post-Processing

This directory contains Python ports of a subset of the Matlab analysis scripts in `../Matlab/`, for producing plots from the "clean" summary `.mat` files that those Matlab scripts generate. Run with the `kynema-env` conda environment.

These scripts output figures that are written to `figures/`.

## RosarioSeaSpiderAnalysis.py

Python port of `../Matlab/RosarioSeaSpiderAnalysis.m`. Loads
`../Matlab/data/SeaSpider_Sig250_Rosario.mat` (Sea Spider seafloor tripod,
uplooking Sig250) and produces:

1. **Average speed profile** vs. height above seafloor (z=0 at seafloor),
   matching the first panel of `SeaSpider_Sig250_Rosario_profiles.png`. Speed
   here is the **horizontal-only** magnitude, $\sqrt{u^2+v^2}$, matching the
   Matlab script's `speed = (u.^2+v.^2).^.5`. A fixed surface-reflection zone
   (z = 83 to 94 m above seafloor, matching the `area([...],[94 94],83)`
   shading in the Matlab script) is shaded red and excluded from the plotted
   profile line.
2. **Time series of speed** at each height in `time_series_heights` (default:
   near-bottom height, configurable near the top of the script). If a height
   equals `hubdepth`, the exact `hubspeed` series from the `.mat` file is
   plotted and labeled "hub depth" instead of interpolating from the profile
   grid.

## RosarioStablemoorAnalysis.py

Python port of `../Matlab/RosarioStablemoorAnalysis.m`. Loads
`../Matlab/data/STBM_Sig500s_Rosario.mat` (STBM mid-water mooring,
up- and down-looking Sig500s) and produces:

1. **Average speed profile** vs. height above seafloor, using the
   fixed-height-above-seabed gridded data (`fixedz_profiles`), matching the
   first panel of `STBM_Sig500s_Rosario_profiles.png`. Speed here is the
   **3D** magnitude, $\sqrt{u^2+v^2+w^2}$, matching the Matlab script's
   `speed = (u.^2+v.^2+w.^2).^.5` (a conservative correction for mooring
   tilt/trim; per the Matlab script's comment this does not change the AEP
   result vs. horizontal-only speed). Same fixed 83–94 m surface-reflection
   shading as the Sea Spider script.
2. **Time series of speed** at each height in `time_series_heights` (default:
   hub depth and near-bottom height), with the same hub-depth handling as
   above.
3. **Surface wave time series** (significant wave height, peak wave period),
   matching the `fh21` plot in the Matlab script. Uses the `Waves` struct
   array in the `.mat` file, which is populated from an externally
   preprocessed `STBM_Sig500up_preprocessed.mat` file (see the "surface
   waves" section of the Matlab script — this preprocessing step requires
   `sigProcess_all.m`/`UVZwaves.m` from the NortekCodes/SWIFTcodes repos,
   which are not part of this workspace). If `Waves` is not present in the
   `.mat` file, this plot is skipped with a printed message.
4. **Wave height vs. peak period** scatter, matching the `fh22` plot.

## Date range selection

Both scripts have user-editable `start_date`, `stop_date`, and `start_time`
variables near the top. Leave them as `None` to use the full data record —
in that case, plot titles show the actual min/max timestamps found in the
data rather than the literal strings "start"/"end".
