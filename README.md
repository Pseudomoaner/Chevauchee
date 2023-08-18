# Chevauchée
### Modelling and image analysis approaches for unpicking motility-enhanced CDI

This repository contains all code necessary to reproduce the results presented in the manuscript 'Cell motility greatly empowers bacterial contact weapons', by Sean C. Booth, Oliver J. Meacock and Kevin R. Foster.

'Donut' and 'plotSpread' are packages used for plotting purposes, available on Matlab central:
https://ch.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot
https://ch.mathworks.com/matlabcentral/fileexchange/56833-donut

The other directories are organised as follows:

## Image analysis
Contains Matlab and R scripts and imageJ macros for processing imaging data from motile surficial colonies. Contains two sub-directories:

### intermixing
Contains scripts to estimate the intermixing and packing fraction of surficial colonies over time. The **Batch_analyseImageContactData.m** script performs analyses using functions in this directory to extract intermixing and packing fraction statistics. **createSeedImg.m** is used to create a visualisation of the initial state of a surficial colony, just after inoculation.

### PIV
Contains scripts to pre-process timelapse videos, run PIVlab particle image velocimetry analysis to measure the cell movement velocity, then extract and analyze the resultant data. **Timelapse_bg_subtract.ijm** is an imageJ script that subtracts the background, improving the image quality. **PIVlab_commandline_batch_process** runs the PIVlab analysis and contains the exact analysis parameters that were used. **PIVlabDataExtract_V3.m** then iterates through the resultant files to extract the data as .csv files. PIV_analyze_V2.R then loads these files in R and computes representative summary statistics (median velocities).

## ODEs
Contains Matlab scripts for running the coarse-grained continuum model. The **Plotting** sub-directory contains scripts used to generate the data visualisations in the manuscript and supplementary movies, while the **DiffusiveModel** sub-directory contains the core modelling code. In the latter, the **diffusiveModel.m** script can be used to generate a single simulation with user-specified settings. The **diffusiveVariableVModel.m** and **diffusiveVariableVModelParamSweep.m** scripts run this model based on experimentally-derived velocity and inoculation density measurements.

## QuasiS3Rruns
Contains Matlab scripts for running parameter sweeps of the QuasiS3R self-propelled rod model (https://github.com/Pseudomoaner/QuasiS3R). The sub-directory **Plotting** contains functions for reproducing all IBM-associated plots in the manuscript.

Directories named 'Archives' contain old versions of code not presented in the manuscript.
