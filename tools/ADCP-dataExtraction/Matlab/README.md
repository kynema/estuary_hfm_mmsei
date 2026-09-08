# Matlab scripts and data from Thomson et al. 2025

This directory includes the Matlab scripts and a placeholder for the data used in the [Thomson et al. 2025 report](https://mhkdr.openei.org/files/609/Rosario_Data_Report_final.pdf). The Matlab scripts are designed to extract and process ADCP data for analysis. The scripts that are necessary for our study are included in this repository, but the data files themselves are not included due to size constraints. The data files can be downloaded from the [OpenEI MHK Data Repository](https://mhkdr.openei.org/submissions/609). Specifically, you'll want to download the following data files into a subdirectory `Matlab/data`:

- STBM Sig 500.mat
- Sea Spider Sig 250.mat

The included Matlab scripts have been edited so the paths work correctly with the data files in the `data` directory. These scripts will load the raw data files from the `avgd` and `burst` directories within [`Nortek Matlab Raw.zip`](https://mhkdr.openei.org/submissions/609), then produce the .mat files STBM Sig 500.mat and Sea Spider Sig 250.mat in the `data` directory. The scripts will also produce a number of figures that are included in the report.

