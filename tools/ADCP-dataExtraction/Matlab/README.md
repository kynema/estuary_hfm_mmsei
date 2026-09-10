# Matlab scripts and data from Thomson et al. 2025

This directory includes the Matlab scripts and a placeholder for the data used in the [Thomson et al. 2025 report](https://mhkdr.openei.org/files/609/Rosario_Data_Report_final.pdf). The Matlab scripts are designed to extract and process raw ADCP data, outputting several figures and cleaned up data files:

- STBM Sig 500.mat
- Sea Spider Sig 250.mat

These scripts are _not_ necessary for our analysis, but they are included in this repository just in case we need them. Additionally, the raw data files are not included due to size constraints. The raw data files (`Nortek Matlab Raw.zip`) can be downloaded from the [OpenEI MHK Data Repository](https://mhkdr.openei.org/submissions/609). The included Matlab scripts have been edited so the paths work correctly with the data files in the `data` directory. These scripts will load the raw data files from the `avgd` and `burst` directories within [`Nortek Matlab Raw.zip`](https://mhkdr.openei.org/submissions/609), then produce the .mat files STBM Sig 500.mat and Sea Spider Sig 250.mat in the `data` directory. The scripts will also produce a number of figures that are included in the report.

