# crf-radiative-forcing-pipeline
The data reduction pipeline for the paper "Cloud Radiative Response to Galactic Cosmic Ray Variability" by Svensmark and Shaviv.

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

NASA MODIS Data:  
https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MOD08_M3/

A list of the files downloaded can be found in "List-of-MODIS-files.txt".

The analysis of the MODIS data with the Fu Liou model is then compared with the CERES data, which can be downloaded from:

CERES DATA:
https://ceres-tool.larc.nasa.gov/ord-tool/jsp/SYN1degEd41Selection.jsp

A list of the files downloaded can be found in "List-of-CERES-files.txt".

---

## Additional Code Requirement

The Fu-Liou model can be downloaded from:
https://web.archive.org/web/20100527145310/http://snowdog.larc.nasa.gov/cgi-bin/rose/flp200503/flp200503.cgi

One of the files simple.f90 was modified, and is named simpleMODIS.f90 and can be found in this respository. 

---

## The Pipeline

The first step is to download the data and model above. Then modify the Fu-Liou code with simpleMODIS.f90.

The next step is to run the IDL program: Extract_Modis_Parameters.pro which produces several files that the Fu Liou code then reads. 

The Fu Liou code generates several result files which contain the geographic and monthly dependence of the forcing, normalized to several galactic cosmic ray reductions.

The additional IDL programs plot the geographic distribution as a function of month. 

---

## Citation

If you use this repository, please cite the original paper:
to-be-inserted-after-publication

