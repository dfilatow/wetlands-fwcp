library(bcdata)##https://github.com/bcgov/bcdata
library(sf)
library(sp)
library(raster)
library(tidyverse)
library(GrabCDED)
library(rgdal)
library(rgrass7)
library(openSTARS)

####Define Global Vars####
bcalb_epsg <- 3005
bcalb_prj4<-'+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs'


####Load in user specified watershed vector layer####
user_ws<-readline("Is a watershed vector layer being provided? Enter T of F. If F, a BC watershed code will be promted for.: ")

print("Reading in watershed layer...")
if(as.logical(user_ws))
{
  ws_pth<-readline("Please provide path to watershed vector layer to use: ")
  
  ws<-st_as_sf(sf::st_read(ws_pth))
  
  if(sf::st_crs(ws)$epsg!=bcalb_epsg)
  {
    ws<-st_transform(ws,bcalb_epsg)
  }
  
  
}else{
  ws_code<-readline("Please provide BC watershed code to use: ")
  
  ws <- bcdc_query_geodata("51f20b1a-ab75-42de-809d-bf415a0f9c62", crs = bcalb_epsg) %>%
    filter(WATERSHED_GROUP_CODE == ws_code) %>%
    collect() %>% st_as_sf()
}


####Get BC CDED tiles intersecting watershed, creat DEM####
print("Obtaining digital elevation model, this may take a while...")
dem_ws<-cded_get(ws)

if(crs(dem_ws)@projargs!=bcalb_prj4)
{
  dem_ws<-raster::projectRaster(dem_ws,crs=bcalb_prj4)
}

####Load in wetland classification results, mask to watershed vector layer####
wwt_pth<-readline("Please provide path to wetland model classification raster, extent assumed to contain watershed AOI: ")
wwt<-raster(wwt_pth)
if(crs(wwt)@projargs!=bcalb_prj4)
{
  wwt<-raster::projectRaster(wwt,crs=bcalb_prj4)
}

#Crop and mask wetland classification to watershed
wwt_ws<-wwt %>%
  raster::crop(ws) %>%
  raster::mask(ws)

#Remove upland pixels
wwt_ws[wwt_ws==1]<-NA


####Convert wwt_ws to a data frame with x-y coords####
wwt_df <- raster::as.data.frame(wwt_ws,xy=T,na.rm=T) %>% 
  dplyr::rename(OID=X20190926.103025_map_recl_OID,Value=X20190926.103025_map_recl_Value,Count=X20190926.103025_map_recl_Count) %>% 
  dplyr::mutate(UID=row_number()) %>%
  dplyr::mutate(out_name = paste('basin_',UID,sep=""))




####Get disturbance / land use layers####

#Create sf lines: stream network, sn for watershed code ws.code
#https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-stream-network
#warning: this takes 1-2 minutes depending on your choice of watershed
bcdc_describe_feature("92344413-8035-4c08-b996-65a9b3f62fca")
sn <- bcdc_query_geodata("92344413-8035-4c08-b996-65a9b3f62fca", crs = bcalb) %>%
  filter( WATERSHED_GROUP_CODE == ws.code) %>%
  collect() %>% st_as_sf()





initGRASS(gisBase = "/usr/lib/grass78",
          home=tempdir(),
          gisDbase ='/home/hunter/grassdata',
          location = 'WetlandTestNewer',
          override = T,
          remove_GISRC=T)

use_sp()

dem_ws_pth<-paste(tempdir(),'/dem_ws.tif',sep="")
writeRaster(dem_ws,dem_ws_pth)

openSTARS::setup_grass_environment(dem_ws_pth)

# 
# wwt_pnt<-rasterToPoints(wwt_ws,spatial = T)
# coordinates(wwt_df) <- ~x+y
# st_write(st_as_sf(wwt_pnt),"/home/hunter/Downloads/wetlnd_pnts.gpkg")
# cl<-makeCluster(25)
# registerDoParallel(cl)
# 
# hm<-foreach(c(i=1:nrow(wwt_df)),.packages=c('rgrass7'),.export=c('wwt_df')) %dopar% + calc_basin(wwt_df[i,1],wwt_df[i,2],wwt_df[i,7])



writeRAST(as(dem_ws,"SpatialGridDataFrame"),'dem',ignore.stderr = T)


execGRASS("r.watershed",parameters = list(elevation='dem',threshold=10,accumulation='acc',tci='topo_idx',spi='strm_pow',drainage='dir',stream='stream_r',length_slope='slop_lngth',slope_steepness='steep'),flags = c('overwrite','quiet'))


delin_basin<-function(basin_df,procs,out_dir)
{
  cmd_tbl<-c()
  
  cmd_tbl[1]<-"#!/bin/bash"
  
  cmds<-as.data.frame(paste('r.stream.basins direction=dir@PERMANENT coordinates=',basin_df$x,",",basin_df$y,' basins=basin_',basin_df$UID,' --overwrite --quiet &',sep=""))
  
  cmd_tbl<-as.data.frame(rbind(cmd_tbl,cmds))
  
  cmd_tbl$row_idx<-c(1:nrow(cmd_tbl))
  
  colnames(cmd_tbl)<-c('cmd','ridx')
  
  wait_v<-(c(1:floor(nrow(basin_df)/procs))*procs)+0.5
  
  wait_v<-as.data.frame(cbind(rep("wait",length(wait_v)),wait_v))
  
  colnames(wait_v)<-c('cmd','ridx')
  
  cmd_tbl<-rbind(cmd_tbl,wait_v)
  
  cmd_tbl$ridx<-as.numeric(cmd_tbl$ridx)
  
  cmd_tbl<-as.data.frame(cmd_tbl)[order(cmd_tbl$ridx),]
  
  cmd_tbl<-cmd_tbl$cmd
  
  script_pth=paste(out_dir,'/grass_delin_basin.sh',sep="")
  
  write(cmd_tbl,ncolumns = 1,script_pth)
  
  system(paste('chmod u+x ',script_pth,sep=""))
  
  cmd_str<-paste('sh ',script_pth,sep="")
  
  print("Delineating basins, this may take a while ...")
  system(cmd_str)
}

delin_basin(wwt_df,22,"/home/hunter/Downloads/Poo")
