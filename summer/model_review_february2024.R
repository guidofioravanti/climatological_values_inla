#Last Review: February 28 2024
#The code relies on "sf" and "terra" packages, removed the dependencies on "sp" and "raster"!
#At the end of the script an example of how create a tif with the mean prediction of the posterior distribution
#and how to draw the maps at low resolution (obviously you need more memory than my laptot)
rm(list=objects())
library("ggplot2")
library("ggspatial")
library("patchwork") 
library("scico") #scientific colour palettes from Fabio Crameri
library("readr")
library("dplyr")
library("purrr")
library("glue")
library("terra")
library("sf")
library("inlabru")
library("INLA")
library("fmesher") #the new INLA way for creating meshes
library("brinla") #for the html report: must be downloaded from github (not from my github)
library("rmarkdown") #for the html report
library("config") #parameters for this script

#####################################################
# sf objects for the Italian regions are in my personal package: regioniItalia
# to install the "regioniItalia" package: 
# devtools::install_github("guidofioravanti/regioniItalia")
#####################################################
library("regioniItalia") 

#load pardiso (follow the instructions from https://panua.ch/pardiso/ to get a free license)
inla.pardiso()

#in my case, I have to set "num.threads=1" to avoid crash in INLA
inla.setOption(inla.mode="compact",num.threads =1) 


####################################################
# The working directory name must be the one of the 4 seasons (lower case):
# winter (December, January and February)
# spring (March, April, May)
# summer (June, July, August)
#autumn  (September, October, November)
#####################################################
basename(getwd())->Season

#According to the Season, the script decides which months we want to elaborate

# Read the configuration file "config.yml". The parameters read from the config file depend on
# the name of the working directory (summer, in this example)... 
purrr::partial(.f=config::get,config=Season,file="./config_july2023.yml")->myget

#PARAM can be: "tmax" (maximum temperature) or "tmin" (minimum temperature)
PARAM<-myget(value="param")

#Other parameters: 
# - number_of_years
# - yearE (End year, last year in the data series)
yearE<-as.integer(myget(value="last_year"))
number_of_years<-as.integer(myget(value="number_of_years"))
yearS<-yearE- number_of_years +1 #Start year, yearS 


#months in the season
MONTHS<-as.integer(myget(value="months"))

#If we want to run the code using only the "Aeronautica Militare" network (this is a national network covering the whole Italian Peninsula, including the islands...but not very dense)
ONLY_AERONAUTICA<-myget(value="only_aeronautica")

#Some meteorological stations are located on small Italian islands
#We decided remove the following small islands from the dataset: Pantelleria, Ponza, Ustica, M. calamita 
ISOLETTE<-c("aeronautica_6897","aeronautica_6841","aeronautica_6874","aeronautica_6790")


#epsg:32632 in KM
CRS("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")->mycrs


#template is a raster of 0s. We use "template" to rasterize the predictions. 
#In this example, we use a low resolution grid (my laptop is old and has not enough memory to manage high resolution grids)
rast("template.tif")->template
terra::res(template)<-10 #low resolution for the final maps

#italia_senza_isole ("Italy without islands") is an sf object from my personal library: regioniItalia (see above)
#italia_senza_isole represents the Italian territory including only the two major islands: Sardinia and Sicily
st_transform(italia_senza_isole,crs = mycrs)->mymask

#Ok..I cannot explain why I created the following sf objects (latitude, dem costa), there is no sense! These should be SpatRaster objects...
#When I wrote the script inlabru supported sp and rast objects, while sf and SpatRaster objects gave me different issues..maybe this explains this choice?
#In any case, be smarter than me! 

#in the following objects "std" indicates that the values have been standardized

#latitude
rast("std_latitudine.tif")->lat

#dem (altitude)
rast("std_dem.tif")->dem

#linear distance from the coast (costa is the Italian for coast) 
rast("std_costa.tif")->costa


#run the model for each month in the season.... purrr::walk, avoid "for" loops
purrr::walk(MONTHS,.f=function(MONTH){ #MONTH is MONTH
  
  #input file name
  glue::glue("dati_{PARAM}_{MONTH}.RDS")->nomeFileDati
  
  #output file name
  glue::glue("inla_{PARAM}_{MONTH}.RDS")->nomeFileOut
  
  readRDS(nomeFileDati) %>%
    filter(!SiteID %in% ISOLETTE) %>% #remove small islands
    filter(yy>=yearS & yy<=yearE) %>% #filter years: thedataset can contain data out of the period of interest
    mutate(yy=yy-yearS+1) %>% #transform yy (years) into values from 1 to number_of_years.. 
    mutate(banda=yy)->df #create a new variable "banda"
  
  
  #if TRUE, run the analysis using only the data from Aeronatutica Militare network (there is no sense in doing that, just in the exploratory analysis)
  if(ONLY_AERONAUTICA){
    
    df %>%
      filter(grepl("^aeronautica_.+",SiteID))->df
    
  }#end if ONLY_AERONAUTICA
  
  #these codes identify bad quality time series which we remove from the input dataset
  to_remove<-c("sardegna_CA069B522","sardegna_NU085B531","sardegna_CA069B522","sardegna_CA048B551","piemonte_5651","liguria_5868","liguria_5818")
  
  df %>%
    filter(! SiteID %in% to_remove)->df
  
  saveRDS(df,glue::glue("subDati_{PARAM}_{MONTH}.RDS"))

  #read the metadata of the stations of interest (the file contains information about stations not available here)
  read_delim(glue::glue("fixed_anagrafica.{PARAM}.csv"),delim=";",col_names = TRUE) %>%
    filter(id %in% unique(df$SiteID)) %>%
    rename(x=Longitude,y=Latitude)->ana #ana is the short version for "anagrafica"
  
  #create an sf objects: 
  st_as_sf(ana,coords=c("x","y"),crs=4326)->sfAna
  st_transform(sfAna,crs=32632)->utm_sfAna #there is no need for this intermediate step
  st_transform(utm_sfAna,crs=mycrs)->utm_sfAna
  
  #in this code we use all the points (Sardinia+Italian Peninsula) together
  st_coordinates(utm_sfAna)->puntiItalia
  
  #We could create the mesh using the stations from the Italian Peninsula  and the stations from Sardinia region separated (two nonconvex.hull objects)
  #st_coordinates(utm_sfAna %>% filter(!grepl("sardegna",regione2,ignore.case = TRUE)))->puntiStivale #Stivale is boot
  #st_coordinates(utm_sfAna %>% filter(grepl("sardegna",regione2,ignore.case = TRUE)))->puntiSardegna  #Sardegna is Sardinia
  #inla.nonconvex.hull(points=puntiStivale)->stivale
  #inla.nonconvex.hull(points=puntiSardegna)->isola
  #INLA::inla.mesh.2d(boundary=list(list(stivale,isola)),max.edge=c(75,200),cutoff = 25,offset=c(5),min.angle=25,crs = mycrs)->mesh_modello
  
  #The above-commented code is the approach we followed for the creation of the mesh for the PM10 data in Italy (see the paper). There we
  #wanted to create two dense meshes in Italy and in Sardinia and a larger mesh over the sea. Here we consider the temperature a continuous field,
  #so we create a unique mesh covering Italy and Sardinia (for Sicily no problem, it is very close to the boot)

  fmesher::fm_mesh_2d_inla(loc.domain=st_as_sfc(utm_sfAna),max.edge=c(50,300),cutoff = 25,offset=c(150),min.angle=25,crs = mycrs)->mesh_modello
  
  #visualize the mesh and the points:
  st_transform(italia_senza_isole,crs=mycrs)->italia_senza_isole_utm
  
  png("visualize_me_this.png",width=1000,height=1000)
  plot(mesh_modello)
  plot(st_geometry(italia_senza_isole_utm),add=TRUE)
  plot(st_geometry(utm_sfAna),add=TRUE,pch=21,bg="red")
  dev.off()
  saveRDS(mesh_modello,glue::glue("mesh_{MONTH}_{PARAM}.RDS"))
  
  
  #PC-prior on the Matern field
  inla.spde2.pcmatern(mesh=mesh_modello,constr = TRUE,prior.range = c(800,0.7),prior.sigma=c(1,0.6))->matern_temperature
 
  
  left_join(df,ana,by=c("SiteID"="id")) %>%
    filter(!is.na(x) & !is.na(y))->df
  
  st_transform(st_as_sf(df,coords=c("x","y"),crs=4326),crs=mycrs)->utm_df
  
  #likelihood
  #mensile is "monthly" in Italian..are the temperature values we want to interpolate
  mylike<-inlabru::like(family = "gaussian",formula=mensile~Intercept+yy+anno+std_dem+std_latitudine+std_costa+spde_field,data=utm_df,include=c("Intercept","yy","anno","std_dem","std_latitudine","std_costa","spde_field"))
  
  #components
  mycmp<- ~ -1+
             Intercept(1)+
             yy+ #an integer from 1 to the number of years (this is for the temporal linear trend)
             anno(yy,model="iid")+ #anno (year) is a label for yy: this is the random effect for years to capture extra variability 
             std_dem(main=dem,model="linear")+ #std_dem is a label
             std_latitudine(main=lat,model="linear")+ #std_latitudine is a label
             std_costa(main=costa,model="linear")+ #std_costa is a label
             spde_field(geometry,model=matern_temperature,group=banda,ngroup=number_of_years,control.group = list(model="iid")) #spde_field is a label
  
             
  bru(mycmp,mylike,options=list(verbose=TRUE))->inla.out
  #saveRDS(inla.out,glue::glue("bru_{MONTH}_{PARAM}.RDS"))
    
  #generate an html report
  rmarkdown::render(input="analisi-covariate.Rmd",output_file = glue::glue("analisi-covariate_{Season}_{MONTH}.html"),params = list(result=inla.out,stagione=Season,mese=MONTH,param=PARAM))

  #predictions: it is important having access to the environment where the model has been defined...otherwise..no way!  
  #for the paper we used dims=c(350,350), here I use 200 x 200 ... otherwise my laptop crashes
  fm_pixels(mesh_modello,dims=c(200,200),mask=st_geometry(mymask))->pixels_italy
  fm_cprod(pixels_italy,data.frame(banda=seq_len(number_of_years),yy=seq_len(number_of_years)))->cpix
  
  #predictions
  pred<-predict(inla.out,cpix,~Intercept+yy+anno+std_dem+std_latitudine+std_costa+spde_field)
  #saveRDS(pred,glue::glue("pred_{MONTH}_{PARAM}.RDS"))
  
  #Ok we have the predictions! WOW!!!
  #create a tif for the mean of the posterior distribution 
  purrr::map(1:number_of_years,.f=function(ix){
    
    pred[pred$banda==ix,c("mean")]->myprediction 
    vect(myprediction)->vprediction
    
    #rasterize the SpatVector
    rasterize(vprediction,template,field="mean",fun="mean")  
    
  }) %>% rast()->mystack
  
  #write a multylayer tif  
  writeRaster(mystack,glue::glue("pred_{MONTH}_{PARAM}.tif"),gdal="COMPRESS=LZW",overwrite=TRUE)
  
  #draw the maps using ggplot,ggspatial, scico and patchwork
  purrr::map(52:60,.f=function(ix){
    
    mystack[[ix]]->mygrid

    ggplot()+
      ggspatial::layer_spatial(mygrid)+
      ggspatial::layer_spatial(italia_senza_isole_utm,fill="transparent")+
      scale_x_continuous(limits=c(0,1500),expand=c(0,0))+
      scale_y_continuous(limits=c(4000,5500),expand=c(0,0))+      
      scale_fill_scico(limits=c(5,40),breaks=seq(5,40,by=5),palette = "lajolla",name="",
                       na.value="transparent",
                       direction=-1,
                       guide=guide_colorbar(barheight = ggplot2::unit(6,"cm")))+
      theme_bw()+
      theme(axis.text = element_blank(),axis.ticks = element_blank(),
            panel.grid = element_blank())
    
    
  })->mapList
  
  reduce(mapList[1:8],.f=`+`)+plot_layout(guides="collect",ncol = 2)->finalPlot
  
  #plot the last 8 maps 
  png(glue::glue("map_mean_{Season}_{MONTH}.png"),width=1024,height=1024)
  print(finalPlot)
  dev.off()
  
})
