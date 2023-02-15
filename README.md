# Master Thesis
Thesis: How does missing data affect measurement error estimation in composite data sets when using MILC: A simulation study
Performed by: Iris van Santen
Supervised by:
Unstitution: 

# Introduction
This repository provides all information and necessary files to replicate the simulation study.

# Software requirements
The whole thesis is coded with R. 

# Instructions
The simulation study consists of three parts:
1. Generate data sets found in folder
2. Apply MILC found in folder
3. Statistical Analysis found in folder

A fourth folder contains scripts to generate the tables and figures given in the thesis report and a fifth folder contains files to reproduce the manuscript of the thesis report.
The repository contains the following files:

Folders/Files | Description
--- | ---
Data_generation | Folder with scripts for data generation
/MAIN_data_generation | Script to generate the data for the simulation study
/FUN_dataset_simulation | Script with function for dataset simulation
MILC_application | Folder with script for MILC application
/MAIN_MILC | Script to apply MILC on the simulate data sets
/FUN_LC | Script with function for LC model estimation
Statistical_analysis | Folder with scripts for the statistical analysis
/MAIN_StatAn | Scripts to perform statistical analysis
/FUN_Performance_Measures | Script with function for the performance measures
/Output_generation | Folder with scripts for the output given in the thesis report

Below are the steps given to reproduce this simulation study:
## Step 1: Data generation 
Folder: Data_generation
1. The working directory should be set in the right folder.

## Step 2: MILC application
Folder: MILC application
1. Move the datafile generated in the previous step in this folder.  

## Step 3: Statistical Analysis
Folder: Statistical Analysis
1. Move the datafiles generated in the previous step in this folder.

## Step 4: Output generation 
Folder: Output_generation
1. Move the datafiles fenerated in the previous step in this folder.

# Contact
For any questions. 
