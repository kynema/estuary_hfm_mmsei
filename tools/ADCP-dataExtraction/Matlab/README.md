# Matlab scripts and data from Thomson et al. 2025

This directory includes the Matlab scripts cleaned data used in the [Thomson et al. 2025 report](https://mhkdr.openei.org/files/609/Rosario_Data_Report_final.pdf). The Matlab scripts are designed to extract and process raw ADCP data, outputting several figures and cleaned up data files:

- STBM Sig 500.mat
- Sea Spider Sig 250.mat

These scripts are _not_ necessary for our analysis, but they are included in this repository just in case we need them. Additionally, the raw data files are not included due to size constraints. The raw data files (`Nortek Matlab Raw.zip`) can be downloaded from the [OpenEI MHK Data Repository](https://mhkdr.openei.org/submissions/609). The included Matlab scripts have been edited so the paths work correctly with the data files in the `data` directory. These scripts will load the raw data files from the `data/avgd` and `data/burst` directories within [`Nortek Matlab Raw.zip`](https://mhkdr.openei.org/submissions/609), then produce the .mat files STBM Sig 500.mat and Sea Spider Sig 250.mat in the `data` directory. The scripts will also produce a number of figures that are included in the report. Note that the OpenEI data repository is missing `STBM_Sig500up_preprocessed.mat`, which is required for the surface wave analysis in `RosarioStablemoorAnalysis.m`. This file was obtained from Jim Thomson and is included in this repository for completeness. 

