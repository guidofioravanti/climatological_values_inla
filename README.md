## Interpolating climate variables by using INLA and the SPDE approach (Last Update: February 28 2024)

This repository contains code and resources that can be used to reproduce the results presented in the paper: "*Interpolating climate variables by using INLA and the SPDE approach*".

### R code

The code filename is: **model_review_february2024.R**. The code is fully commented in order to make it as clear as possible.
The code has been updated in February 2024. 

- The original code relied on the R "raster" and "sp" libraries.
- Tne new code uses the R "terra" and "sf" libraries.
- The new code uses the "fmesher" package for the mesh definition.
- The covariates are stored as raster "tif" files (instead of SpatialPixelsDataFrame objects).

The script illustrates how:
- to define the model by using the R "inlabru" package
- to predict using the "predict" function from "inlabru"
- to create tif output files and draw some maps.

In some parts, the R script is not very elegant but (we hope) it should be useful.

### Input data

#### Climatological values from in-situ stations

The datasets provided here (dati_tmax_6.RDS, dati_tmax_7.RDS, dati_tmax_8.RDS) contain the maximum temperature (tmax) monthly series for the summer season (June, July and August). The stations' metadata are available in text format (fixed_anagrafica.tmax.csv).  

#### Covariates

The covariates are encoded as ".tif" raster files:

- **std_costa.tif**: linear distance from the coast
- **std_dem.tif**: digital elevation model
- **std_latitudine.tif**: latitude raster
  
Note: the covariates are standardized ("std"). 

#### Other files

- **template.tif**: this raster file serves as a template which defines the study domain extent and the spatial resolution of the output maps
- **config_july2023.yml**: config file (some parameters for this script)
- **analisi-covariate.Rmd**: an RMarkdown document to document the main model outputs (this .Rmd script uses the "brinla" package)

#### Last but not least

The R script has been tested using:
- R version 4.3.2 (2023-10-31) -- "Eye Holes"
- INLA (24.2.9)
- inlabru (2.10.1)
- fmesher (0.1.5)
- terra (1.7.71)
- sf (1.0.15)

### Authors :writing_hand:

* Guido Fioravanti, European Commission's Joint Research Center <guido.fioravanti@ec.europa.eu>
* Sara Martino, Norwegian University of Science and Technology, Norway
* Michela Cameletti, Universita degli Studi di Bergamo, Italy
* Andrea Toreti, European Commission's Joint Research Center


