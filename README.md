# crf-radiative-forcing-pipeline
The data reduction pipeline for the paper Cloud Radiative Response to Galactic Cosmic Ray Variability by Svensmark and Shaviv

---

## Input Data

The pipeline uses **monthly MODIS satellite products** containing atmospheric and cloud parameters relevant for radiative transfer calculations.

Typical variables include:

- Cloud optical depth
- Cloud fraction
- Effective radius
- Surface albedo
- Atmospheric state parameters

Due to file size constraints, the MODIS datasets are **not included** in this repository.

They can be obtained from:

NASA LAADS DAAC  
https://ladsweb.modaps.eosdis.nasa.gov/

---

## Program Components

### 1. Parameter Extraction

This program reads MODIS monthly datasets and extracts the variables required by the radiative transfer model.

It performs:

- Data selection and filtering
- Unit conversion
- Reformatting into model input format

**Output:** Model-ready parameter files.

---

### 2. Radiative Transfer Calculations

Radiative fluxes are computed using the **Fu–Liou radiative transfer model**.

This step:

- Ingests extracted atmospheric parameters
- Allows modification of greenhouse gas concentrations
- Computes upward and downward radiative fluxes

**Output includes:**

- TOA shortwave and longwave fluxes
- Net radiative forcing

The model can be downloaded from:
http://insert.link

---

### 3. Output Processing

Scripts compute:

- Net TOA forcing fields
- Temporal averages
- Diagnostic radiative quantities

---

### 4. Plotting

Plotting scripts reproduce the figures presented in the paper, including:

- Global forcing maps
- Time-series analysis
- Regional forcing comparisons

---

## Requirements

The code requires:

- IDL 
- Fu–Liou radiative transfer model installation
  
---

## Reproducibility Instructions

To reproduce the results:

1. Download the required MODIS monthly datasets.
2. Run the parameter extraction program.
3. Execute the radiative transfer model.
4. Process output files to compute net TOA forcing.
5. Run plotting scripts to generate figures.

---

## Data Availability

MODIS datasets are publicly available from NASA archives:

https://ladsweb.modaps.eosdis.nasa.gov/

---

## Citation

If you use this repository, please cite the original paper:
ZZZ

