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

grass_Db<-readline("Please provide path of GRASS GIS gisDbase, e.g., '/home/hunter/grassdata': ")
grass_loc<-readline("Please provide desired name of GRASS GIS Location, e.g.,Wetland_PARS: ")


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

#Resample DEM to wetland classification resolution
dem_ws<-raster::resample(dem_ws,wwt)

#Crop and mask wetland classification to watershed
wwt_ws<-wwt %>%
  raster::crop(ws) %>%
  raster::mask(ws)


####Convert wwt_ws to a data frame with x-y coords####

xyz <- as.data.frame(rasterToPoints(wwt_ws))

colnames(xyz)<-c("x","y","Value")

wwt_df<-xyz %>%
  dplyr::mutate(UID=row_number()) %>%
  dplyr::mutate(out_name = paste('basin_',UID,sep="")) %>%
  dplyr::filter(Value==2 | Value==3)

####Get disturbance / land use layers####

#Create sf lines: stream network, sn for watershed code ws.code
#https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-stream-network
#warning: this takes 1-2 minutes depending on your choice of watershed
bcdc_describe_feature("92344413-8035-4c08-b996-65a9b3f62fca")
sn <- bcdc_query_geodata("92344413-8035-4c08-b996-65a9b3f62fca") %>%
  filter( WATERSHED_GROUP_CODE == ws_code) %>%
  collect() %>% st_as_sf()



####Set of GRASS Environment####
initGRASS(gisBase = "/usr/lib/grass78",
          home=tempdir(),
          gisDbase =grass_Db,
          location = grass_loc,
          override = T,
          remove_GISRC=T)

use_sp()

dem_ws_pth<-paste(tempdir(),'/dem_ws.tif',sep="")
writeRaster(dem_ws,dem_ws_pth)

#Set g.region to that of dem_ws using a convenient openSTARS function 
openSTARS::setup_grass_environment(dem_ws_pth)

#Write DEM and wetland classification to wetland environment
writeRAST(as(dem_ws,"SpatialGridDataFrame"),'dem',ignore.stderr = T,overwrite = T)
writeRAST(as(wwt_ws,"SpatialGridDataFrame"),'wwt',ignore.stderr = T, overwrite = T)

#Run r.watershed in GRASS to create DEM derivitives and stream network 
execGRASS("r.watershed",parameters = list(elevation='dem',threshold=10,accumulation='acc',tci='topo_idx',spi='strm_pow',drainage='dir',stream='stream_r',length_slope='slop_lngth',slope_steepness='steep'),flags = c('overwrite','quiet'))

#Set GRASS mask to classification raster 
execGRASS("r.mask",parameters = list(raster='wwt'))


execGRASS("r.mapcalc",flags = c("overwrite"),
          parameters = list(
            expression="wwt_acc = acc"
          ))

execGRASS("r.mask", flags = c("r"))


wwt_acc<-readRAST('wwt_acc',ignore.stderr = T)

wwt_acc<-raster(wwt_acc) %>%
  raster::crop(ws) %>%
  raster::mask(ws)

wwt_acc<-as.data.frame(wwt_acc,xy=T,na.rm=T)

wwt_acc<-wwt_acc[complete.cases(wwt_acc),]

wwt_acc<-wwt_acc[sample(c(1:nrow(wwt_acc)),nrow(wwt_acc)),]

wwt_acc<-wwt_acc %>%
  dplyr::mutate(UID=row_number()) %>%
  dplyr::mutate(out_name = paste('basin_',UID,sep=""))

delin_basin<-function(basin_df,procs)
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
  
  cmd_tbl<-c(cmd_tbl,"exit 0")
  
  script_pth=paste(tempdir(),'/grass_delin_basin.sh',sep="")
  
  write(cmd_tbl,ncolumns = 1,script_pth)
  
  system(paste('chmod u+x ',script_pth,sep=""))
  
  cmd_str<-paste('sh ',script_pth,sep="")
  
  print("Delineating basins, this may take a while ...")
  system(cmd_str)
}

#wwt_df_samp<-wwt_df[sample(c(1:nrow(wwt_df)),100),]

#delin_basin(wwt_df_samp,26)

calc_basin_stats<-function(basin_df,stat_rast,procs,out_dir)
{
  cmd_tbl<-c()
  
  cmd_tbl[1]<-"#!/bin/bash"
  
  cmds<-as.data.frame(paste('r.univar -g --overwrite map=',stat_rast,' zones=basin_',basin_df$UID,' output=',out_dir,'/basin_stats_',basin_df$UID,'.txt separator=newline &',sep=""))
  
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
  
  cmd_tbl<-c(cmd_tbl,"exit 0")
  
  script_pth=paste(tempdir(),'/grass_basin_stats.sh',sep="")
  
  write(cmd_tbl,ncolumns = 1,script_pth)
  
  system(paste('chmod u+x ',script_pth,sep=""))
  
  cmd_str<-paste('sh ',script_pth,sep="")
  
  print("Calculating basin stats, this may take a while ...")
  system(cmd_str)
}

#calc_basin_stats(wwt_df_samp,'slop_lngth',26,"/home/hunter/Downloads/Poo")

stats_to_df<-function(basin_df,stat_dir,stat)
{
  basin_vec<-c()
  stat_vec<-c()
  
  for(r in c(1:nrow(basin_df)))
  {
    stats<-readr::read_lines(paste(stat_dir,'/basin_stats_',basin_df$UID[r],'.txt',sep=""))
    
    basin_vec[r]<-basin_df$UID[r]
    
    if(stat=="n")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[2],"=")[[1]][2])
    }
    if(stat=="min")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[5],"=")[[1]][2])
    }
    if(stat=="max")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[6],"=")[[1]][2])
    }
    if(stat=="mean")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[8],"=")[[1]][2])
    }
    if(stat=="MAE")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[9],"=")[[1]][2])
    }
    if(stat=="stddev")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[10],"=")[[1]][2])
    }
    if(stat=="var")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[11],"=")[[1]][2])
    }
    if(stat=="sum")
    {
      stat_vec[r]<-as.numeric(strsplit(stats[13],"=")[[1]][2])
    }
    
  }
  
  rslt<-as.data.frame(cbind(basin_vec,stat_vec))
  
  colnames(rslt)<-c("UID",stat)
  
  return(rslt)
  
}


get_basin_stats<-function(basin_df,procs,stat_rast,stat,chunk)
{
  N<-nrow(basin_df)
  
  start_idx<-1
  end_idx<-1
  
  stat_vec<-c()
  
  while(end_idx<=N)
  {
    if((start_idx+chunk)>N)
    {
      end_idx<-N
    }else{
      end_idx<-start_idx+chunk
    }
    
    basin_df_sub<-basin_df[c(start_idx:end_idx),]
    
    delin_basin(basin_df_sub,procs)
    
    stat_dir<-paste(tempdir(),"/temp_basin_stats/",sep="")
    system(paste("mkdir ",stat_dir,sep=""))
    
    calc_basin_stats(basin_df_sub,stat_rast,procs,stat_dir)
    
    tmp<-stats_to_df(basin_df_sub,stat_dir,stat)
    
    stat_vec<-rbind(stat_vec,tmp)
    
    execGRASS('g.remove',parameters = list(type='raster',pattern="basin*"),flags=c('f'))
    
    system(paste("rm -r ",stat_dir,sep=""))
    
    start_idx<-end_idx+1
    
  }
  
}
