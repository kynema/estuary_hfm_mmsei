# ADCP Data Extraction Tools
This directory contains a set of tools for extracting and processing ADCP data for analysis. The tools are designed to work with the data files provided in the [Thomson et al. 2025 report](https://mhkdr.openei.org/files/609/Rosario_Data_Report_final.pdf). 

## Matlab Scripts
The Matlab scripts included in the `Matlab` subdirectory are used to load the raw data files from the `avgd` and `burst` directories within [`Nortek Matlab Raw.zip`](https://mhkdr.openei.org/submissions/609), and produce the .mat files STBM Sig 500.mat and Sea Spider Sig 250.mat in the `data` directory. The scripts will also produce a number of figures that are included in the report.

## Python Scripts
The Python scripts included in the `Python` subdirectory are used to load the .mat files STBM Sig 500.mat and Sea Spider Sig 250.mat from the `Matlab/data`, and extract desired data for comparison to the kynema-sgf simulation data at over a specified range of dates.