####Required libraries####
library(bcdata)
library(sf)
library(sp)
library(raster)
library(tidyverse)
library(rgdal)
library(rgrass7)

####Custom libraries####
#Package for calculating upstream statistics (devtools::install_github("HunterGleason/GRASSBasinStats"))
library(GRASSBasinStats)
#Package for obtaining bcdata layers in batch (devtools::install_github("HunterGleason/GetBCDataBatch"))
library(GetBCDataBatch)
#Package for obtaining CDED DEM, may be incorporated into bcdata package devtools::install_github(????))
library(GrabCDED)

####Define Global Vars####
#BC Albs
bcalb_epsg <- 3005
bcalb_prj4<-'+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs'


####Load in user specified watershed vector layer, or use BC watershed code####
user_ws<-readline("Is a watershed vector layer being provided? Enter T of F. If F, a BC watershed code will be promted for.: ")

#Prompt for grass data base path 
grass_Db<-readline("Please provide path of GRASS GIS gisDbase, e.g., '/home/hunter/grassdata': ")
#Promt user for name of the grass location
grass_loc<-readline("Please provide desired name of GRASS GIS Location, e.g.,Wetland_PARS: ")


print("Reading in watershed layer...")
#If user has basin shapefile to use 
if(as.logical(user_ws))
{
  #Get path to basin shapefile
  ws_pth<-readline("Please provide path to watershed vector layer to use: ")
  
  #Read basin as sf 
  ws<-st_as_sf(sf::st_read(ws_pth))
  
  #Reproject if required
  if(sf::st_crs(ws)$epsg!=bcalb_epsg)
  {
    ws<-st_transform(ws,bcalb_epsg)
  }
  
#If use is using a BC watershed code 
}else{
  #Prompt for code 
  ws_code<-readline("Please provide BC watershed code to use: ")
  
  #Use bcdata package to download basin as sf 
  ws <- bcdc_query_geodata("51f20b1a-ab75-42de-809d-bf415a0f9c62", crs = bcalb_epsg) %>%
    filter(WATERSHED_GROUP_CODE == ws_code) %>%
    collect() %>% st_as_sf()
}


####Get BC CDED tiles intersecting watershed, creat DEM####
print("Obtaining digital elevation model, this may take a while...")

#Uses CDED package to get DEM intersecting the study basin 
dem_ws<-cded_get(ws)

#Reproject if required 
if(crs(dem_ws)@projargs!=bcalb_prj4)
{
  dem_ws<-raster::projectRaster(dem_ws,crs=bcalb_prj4)
}

####Load in wetland classification results, mask to watershed vector layer####
#Prompt user for path to wetland classification raster (3-class) 
wwt_pth<-readline("Please provide path to wetland model classification raster, extent assumed to contain watershed AOI: ")

#Read classification using raster package 
wwt<-raster(wwt_pth)

#Reproject if required. 
if(crs(wwt)@projargs!=bcalb_prj4)
{
  wwt<-raster::projectRaster(wwt,crs=bcalb_prj4)
}

#Resample DEM to wetland classification resolution
dem_ws<-raster::resample(dem_ws,wwt)

#Crop and mask wetland classification to watershed
wwt_ws<-wwt %>%
  raster::crop(ws) %>%
  raster::mask(ws)


####Convert wwt_ws to a data frame with x-y coords, filter out upland pixels####

xyz <- as.data.frame(rasterToPoints(wwt_ws))

colnames(xyz)<-c("x","y","Value")

wwt_df<-xyz %>%
  dplyr::mutate(UID=row_number()) %>%
  dplyr::mutate(out_name = paste('basin_',UID,sep="")) %>%
  dplyr::filter(Value==2 | Value==3)

####Get disturbance / land use layers####

#Streams
bcdc_describe_feature("92344413-8035-4c08-b996-65a9b3f62fca")

#Describe current disturbance layers of interest
#Cutblocks
bcdc_describe_feature("b1b647a6-f271-42e0-9cd0-89ec24bce9f7")
#Roads
bcdc_describe_feature("bb060417-b6e6-4548-b837-f9060d94743e")
#Fires
bcdc_describe_feature("22c7cb44-1463-48f7-8e47-88857f207702")

#Define layer parameter vectors for input to GetBCDataBatch package 
sn<-c("92344413-8035-4c08-b996-65a9b3f62fca",NA,"WATERSHED_GROUP_CODE == 'PARS'")
cb.young<-c("b1b647a6-f271-42e0-9cd0-89ec24bce9f7","ws","HARVEST_YEAR > 2010")
cb.adult<-c("b1b647a6-f271-42e0-9cd0-89ec24bce9f7","ws","HARVEST_YEAR <= 2010 && HARVEST_YEAR > 2000")
cb.mature<-c("b1b647a6-f271-42e0-9cd0-89ec24bce9f7","ws","HARVEST_YEAR <= 2000")
rd<-c("bb060417-b6e6-4548-b837-f9060d94743e","ws",NA)
fr.young<-c("22c7cb44-1463-48f7-8e47-88857f207702","ws","FIRE_YEAR > 2010")
fr.adult<-c("22c7cb44-1463-48f7-8e47-88857f207702","ws","FIRE_YEAR <= 2010 && FIRE_YEAR > 2000")
fr.mature<-c("22c7cb44-1463-48f7-8e47-88857f207702","ws","FIRE_YEAR <= 2000")

#Create layer data frame and rename columns 
lyr_tabl<-as.data.frame(rbind(sn,cb.young,cb.adult,cb.mature,rd,fr.young,fr.adult,fr.mature))
colnames(lyr_tabl)<-c('code','geom_ext','filter_exp')

#DL layers using GetBCDataBatch
dist_lst<-GetBCDataBatch::get_bcdata_batch(lyr_tabl)

####Rasterize disturbance layers for calculating basing statistics####
sn<-fasterize::fasterize(dist_lst[[1]] %>% sf::st_buffer(dist=(sum(res(dem_ws)/2))), dem_ws, background = 0) %>% raster::crop(ws)
cb.young<-fasterize::fasterize(dist_lst[[2]],dem_ws, background = 0) %>% raster::crop(ws)
cb.adult<-fasterize::fasterize(dist_lst[[3]],dem_ws, background = 0) %>% raster::crop(ws)
cb.mature<-fasterize::fasterize(dist_lst[[4]],dem_ws, background = 0) %>% raster::crop(ws)
rd<-fasterize::fasterize(dist_lst[[5]] %>% sf::st_buffer(dist=(sum(res(dem_ws)/2))), dem_ws, background = 0) %>% raster::crop(ws)
fr.young<-fasterize::fasterize(dist_lst[[6]],dem_ws, background = 0) %>% raster::crop(ws)
fr.adult<-fasterize::fasterize(dist_lst[[7]],dem_ws, background = 0) %>% raster::crop(ws)
fr.mature<-fasterize::fasterize(dist_lst[[8]],dem_ws, background = 0) %>% raster::crop(ws)

####Set of GRASS Environment####
initGRASS(gisBase = "/usr/lib/grass78",
          home=tempdir(),
          gisDbase =grass_Db,
          location = grass_loc,
          mapset = 'PERMANENT',
          override = T,
          remove_GISRC=T)

#Set g.region to 'dem_ws' parameters 
use_sp()
#write dem_ws to GRASS env. 
writeRAST(as(dem_ws,"SpatialGridDataFrame"),'dem',ignore.stderr = T,overwrite = T)
#Set grass region parameters to DEM layer 
execGRASS('g.region',parameters = list(raster='dem'))


#Write DEM, wetland classification and FWA stream raster to GRASS environment
writeRAST(as(wwt_ws,"SpatialGridDataFrame"),'wwt',ignore.stderr = T, overwrite = T)
writeRAST(as(sn,"SpatialGridDataFrame"),'fwa_stream_rast',ignore.stderr = T, overwrite = T)

#Write binary disturbance raters to GRASS environment 
writeRAST(as(rd,"SpatialGridDataFrame"),'roads',ignore.stderr = T,overwrite = T)
writeRAST(as(cb.young,"SpatialGridDataFrame"),'cb_young',ignore.stderr = T,overwrite = T)
writeRAST(as(cb.adult,"SpatialGridDataFrame"),'cb_adult',ignore.stderr = T,overwrite = T)
writeRAST(as(cb.mature,"SpatialGridDataFrame"),'cb_mature',ignore.stderr = T,overwrite = T)
writeRAST(as(fr.young,"SpatialGridDataFrame"),'fr_young',ignore.stderr = T,overwrite = T)
writeRAST(as(fr.adult,"SpatialGridDataFrame"),'fr_adult',ignore.stderr = T,overwrite = T)
writeRAST(as(fr.mature,"SpatialGridDataFrame"),'fr_mature',ignore.stderr = T,overwrite = T)

####Run r.watershed in GRASS to create DEM derivatives and stream network####
execGRASS("r.watershed",parameters = list(elevation='dem',threshold=10,accumulation='acc',tci='topo_idx',spi='strm_pow',drainage='dir',stream='stream_r',length_slope='slop_lngth',slope_steepness='steep'),flags = c('overwrite','quiet'))

####Example of summarizing basin statistics for sample (n=1000) of wetland pixels for the binary cutblock(adult) layer, use 26 threads####

#sample wwt_df pixels (n=1000)
wwt_df_samp<-wwt_df[sample(c(1:nrow(wwt_df)),1000),]

#Example of GRASSBasinStats use for binary cutblock(adult) layer 
test<-GRASSBasinStats::get_basin_stats(basin_df=wwt_df_samp,procs=26,stat_rast='cb_adult')


####Where from here, convert pixels to polygons, clumping ????####


