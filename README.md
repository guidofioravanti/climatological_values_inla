## Interpolating climate variables by using INLA and the SPDE approach

This repository contains code and resources that can be used to reproduce the results presented in the paper: "*Interpolating climate variables by using INLA and the SPDE approach*".

### R code

The code filename is: **model_july2023.R**. The code is fully commented in order to make it as clear as possible.

### Input data

#### Climatological values from in-situ stations

The datasets provided here (dati_tmax_6.RDS, dati_tmax_7.RDS, dati_tmax_8.RDS) contain the maximum temperature (tmax) monthly series for the summer season (June, July and August). The stations' metadata are available in text format (fixed_anagrafica.tmax.csv).  

#### Covariates

The covariates are encoded as ".RDS" files:

- **std_costa.RDS**: linear distance from the coast
- **std_dem.RDS**: digital elevation model
- **std_latitudine.RDS**: latitude raster
  
Note: the covariates are standardized. 

#### Other files

- **griglia.tif**: this raster file serves as a template which defines the study domain extent and the spatial resolution of the output maps
- **config_july2023.yml**: config file
- **analisi-covariate.Rmd**: an RMarkdown document to document the main model outputs

### Authors :writing_hand:

* Guido Fioravanti, European Commission's Joint Research Center <guido.fioravanti@ec.europa.eu>
* Sara Martino, Norwegian University of Science and Technology, Norway
* Michela Cameletti, Universita degli Studi di Bergamo, Italy
* Andrea Toreti, European Commission's Joint Research Center


