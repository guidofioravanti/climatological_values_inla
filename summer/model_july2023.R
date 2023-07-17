rm(list=objects())
library("INLA")
library("tidyverse")
library("sf")
library("sp")
library("raster")
library("inlabru")
library("config")

#####################################################
# sf objects for the Italian regions are in my personal package: regioniItalia
# to install it run the following command: 
# devtools::install_github("guidofioravanti/regioniItalia")
#####################################################
library("regioniItalia") 

####################################################
# The working directory must be the name of the season:
# winter (December, January and February)
# spring (March, April, May)
# summer (June, July, August)
#autumn  (September, October, November)
#####################################################
basename(getwd())->Season
# read the configuration file "config.yml". The parametres read from the config file depend on
# the name of the working directory (summer, in this example)
purrr::partial(.f=config::get,config=Season,file="./config_july2023.yml")->myget

#PARAM can be: "tmax" (maximum temperature) or "tmin" (minimum temperature)
PARAM<-myget(value="param")
#This is not relevant for running the code: MODELLO<-myget(value="model")
number_of_years<-as.integer(myget(value="number_of_years"))
annoF<-as.integer(myget(value="last_year"))
annoI<-annoF-number_of_years+1 #first year
MONTHS<-as.integer(myget(value="months"))

#if we want to run the code using only the "Aeronautica Militare" network
ONLY_AERONAUTICA<-myget(value="only_aeronautica")

# We remove the small islands from the dataset: Pantelleria, Ponza, Ustica, M. calamita 
ISOLETTE<-c("aeronautica_6897","aeronautica_6841","aeronautica_6874","aeronautica_6790")

#load the pardiso license
inla.setOption(pardiso.license = "~/pardiso/licenza.txt")
bru_options_set(bru_verbose = 1)

CRS("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")->mycrs

#the code should be updated by using the "terra" package!
#griglia (grid) is a raster template on only 0s
raster::raster("griglia.tif")->griglia
#italia_senza_isole is an sf object from my personal library: regioniItalia
#italia_senza_isole means "Italy without islands"
st_transform(italia_senza_isole,crs = mycrs)->maschera

#latitude
readRDS("std_latitudine.RDS")->lat
#dem
readRDS("std_dem.RDS")->dem
#linear distance from the coast 
readRDS("std_costa.RDS")->costa

#run the model for each month in the season
purrr::walk(MONTHS,.f=function(MESE){ #MESE is MONTH
  
  #input file name
  glue::glue("dati_{PARAM}_{MESE}.RDS")->nomeFileDati
  #output file name
  glue::glue("inla_{PARAM}_{MESE}.RDS")->nomeFileOut
  
  readRDS(nomeFileDati) %>%
    filter(!SiteID %in% ISOLETTE) %>% #remove small islands
    filter(yy>=annoI & yy<=annoF) %>% #filter years
    mutate(yy=yy-annoI+1) %>% #transform yy (years) into values from 1 to number_of_years
    mutate(banda=yy)->df #create a new variable "banda"

  if(ONLY_AERONAUTICA){
    
    df %>%
      filter(grepl("^aeronautica_.+",SiteID))->df
    
  }#end if ONLY_AERONAUTICA
  
  #these codes identify bad quality time series which we remove from the input dataset
  codici_eliminare<-c("sardegna_CA069B522","sardegna_NU085B531","sardegna_CA069B522","sardegna_CA048B551","piemonte_5651","liguria_5868","liguria_5818")
  
  df %>%
    filter(! SiteID %in% codici_eliminare)->df
  
  saveRDS(df,glue::glue("subDati_{PARAM}_{MESE}.RDS"))

  #read stations metadata
  read_delim(glue::glue("fixed_anagrafica.{PARAM}.csv"),delim=";",col_names = TRUE) %>%
    filter(id %in% unique(df$SiteID)) %>%
    rename(x=Longitude,y=Latitude)->ana
  
  st_as_sf(ana,coords=c("x","y"),crs=4326)->sfAna
  st_transform(sfAna,crs=32632)->utm_sfAna
  st_transform(utm_sfAna,crs=mycrs)->utm_sfAna
  
  #we could create the mesh using the stations from the Italian Peninsula  and the stations from Sardinia region separated (two nonconvex.hull objects)
  st_coordinates(utm_sfAna %>% filter(!grepl("sardegna",regione2,ignore.case = TRUE)))->puntiStivale
  st_coordinates(utm_sfAna %>% filter(grepl("sardegna",regione2,ignore.case = TRUE)))->puntiSardegna
  #in this code we use all the points (Sardinia+Italian Peninsula) together
  st_coordinates(utm_sfAna)->puntiItalia
  
  #inla.nonconvex.hull(points=puntiStivale)->stivale
  #inla.nonconvex.hull(points=puntiSardegna)->isola
  #INLA::inla.mesh.2d(boundary=list(list(stivale,isola)),max.edge=c(75,200),cutoff = 25,offset=c(5),min.angle=25,crs = mycrs)->mesh_modello
  #---> originale INLA::inla.mesh.2d(loc=puntiItalia,max.edge=c(75,300),cutoff = 50,offset=c(5),min.angle=25,crs = mycrs)->mesh_modello
  INLA::inla.mesh.2d(loc=puntiItalia,max.edge=c(50,300),cutoff = 25,offset=c(150),min.angle=25,crs = mycrs)->mesh_modello
  saveRDS(mesh_modello,glue::glue("mesh_{MESE}_{PARAM}.RDS"))

  inla.spde2.pcmatern(mesh=mesh_modello,constr = TRUE,prior.range = c(800,0.7),prior.sigma=c(1,0.6))->spde
 
  left_join(df,ana,by=c("SiteID"="id")) %>%
    filter(!is.na(x) & !is.na(y))->df
  
  coordinates(df)=~x+y
  proj4string(df)<-CRS("+init=epsg:4326")
  spTransform(df,CRS("+init=epsg:32632"))->utm_df
  spTransform(utm_df,mycrs)->utm_df
  
  utm_df$Intercept<-1
  
  if(!file.exists(glue::glue("bru_{MESE}_{PARAM}.RDS"))){
    bru(mensile~-1+
            Intercept+
            yy+
            anno(yy,model="iid")+
            std_dem(main=dem,model="linear")+
            std_latitudine(main=lat,model="linear")+
            std_costa(main=costa,model="linear")+
            #stazione(SiteID,model="iid")+
            latent(coordinates,model=spde,group=banda,ngroup=number_of_years,control.group = list(model="iid")),
         data=utm_df,
         family = "gaussian",
         options = list(verbose=F,control.compute=list(cpo=F,dic=F,waic=F)))->inla.out
    
    saveRDS(inla.out,glue::glue("bru_{MESE}_{PARAM}.RDS"))
    
  }else{
    readRDS(glue::glue("bru_{MESE}_{PARAM}.RDS"))->inla.out
  }  
  
  rmarkdown::render(input="analisi-covariate.Rmd",output_file = glue::glue("analisi-covariate_{stagione}_{MESE}.html"),params = list(result=inla.out,stagione=stagione,mese=MESE,param=PARAM))
  
  pixels(mesh_modello,nx=350,ny=350,mask=as_Spatial(maschera))->pixels_italia  
  cprod(pixels_italia,data.frame(banda=seq_len(number_of_years),yy=seq_len(number_of_years)))->pixels_italia_tutti
  #predictions
  predizioni <-predict(inla.out,pixels_italia_tutti,~data.frame(banda=banda,yy=yy,myformula=Intercept+yy+anno+std_dem+std_latitudine+latent))  
  saveRDS(predizioni,glue::glue("predizioni_{MESE}_{PARAM}.RDS"))
  
})
