# ADCP Data Extraction Tools
This directory contains a set of tools for extracting and processing ADCP data for analysis. The tools are designed to work with the data files provided in the [Thomson et al. 2025 report](https://mhkdr.openei.org/files/609/Rosario_Data_Report_final.pdf). 

The ADCP's used in the study are a Sea Spider seafloor tripod (uplooking Sig250) and a STBM mid-water mooring (up- and down-looking Sig500s), whose locations are shown in the table below:

| Device | Latitude | Longitude |
|--------|----------|-----------|
| Sea Spider | 48.56270 | -122.76510 |
| Stablemoor (STBM) | 48.56260 | -122.76690 |

## Matlab Scripts
The Matlab scripts included in the `Matlab` subdirectory are used to load the raw data files from the `avgd` and `burst` directories within [`Nortek Matlab Raw.zip`](https://mhkdr.openei.org/submissions/609), and produce the .mat files STBM Sig 500.mat and Sea Spider Sig 250.mat in the `data` directory. The scripts will also produce a number of figures that are included in the report.

## Python Scripts
The Python scripts included in the `python` subdirectory are used to load the .mat files `STBM Sig 500.mat` and `Sea Spider Sig 250.mat` from `Matlab/data`, and extract desired data for comparison to the kynema-sgf simulation data at over a specified range of dates. Additionally, there are scripts to extract and plot the SGF simulation data from the `post_processing` directory of a kynema-sgf simulation run. The Python scripts are designed to be run in the `kynema-env` conda environment.

- `adcp_data_extraction.py` / `plot_adcp_data.py`: load and plot the Rosario ADCP datasets (SS and STBM).
- `sgf_data_extraction.py` / `plot_sgf_data.py`: load and plot Kynema-SGF simulation line-sampling output (SS and STBM LineSamplers) for comparison to the ADCP data. Kynema-SGF defines z=0 at the water surface, so `sgf_data_extraction.py` shifts the sampled z-coordinates by the local still-water depth (STBM uses the fixed surveyed depth, SS uses the mean measured ADCP waterdepth) so that both datasets share a common z=0 at the seafloor. Edit the `sgf_unit`, `SGF_DATA_DIR`, `start_date`, and `start_time` constants at the top of `plot_sgf_data.py` to select the unit, the simulation's `post_processing` directory, and (optionally) the real-world UTC datetime corresponding to simulation time=0.