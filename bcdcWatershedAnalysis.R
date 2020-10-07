####Required libraries####
library(bcdata)
library(sf)
library(sp)
library(raster)
library(tidyverse)
library(rgdal)
library(rgrass7)
library(stars)
library(svMisc)

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
# if(crs(dem_ws)@projargs!=bcalb_prj4)
# {
#   dem_ws<-raster::projectRaster(dem_ws,crs=bcalb_prj4)
# }

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


wwt_ws[wwt_ws==1] <- NA

wwt_ws_clump <- raster::clump(wwt_ws,directions=8)

wwt_ws_v <- raster::rasterToPolygons(wwt_ws_clump,n=8,na.rm=T,dissolve = T)

wwt_ws_v <- wwt_ws_v %>% sf::st_as_sf()

wwt_ws_v <- wwt_ws_v %>% 
  dplyr::mutate(Area = sf::st_area(.)) %>%
  dplyr::filter(as.numeric(Area)>=5000) %>%
  dplyr::rename(UID=clumps)


####Convert wwt_ws to a data frame with x-y coords, filter out upland pixels####

# xyz <- as.data.frame(rasterToPoints(wwt_ws))
# 
# colnames(xyz)<-c("x","y","Value")
# 
# wwt_df<-xyz %>%
#   dplyr::mutate(UID=row_number()) %>%
#   dplyr::mutate(out_name = paste('basin_',UID,sep="")) %>%
#   dplyr::filter(Value==2 | Value==3)

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

#Write DEM, wetland classification and FWA stream raster to GRASS environment
writeRAST(as(wwt_ws,"SpatialGridDataFrame"),'wwt',ignore.stderr = T, overwrite = T)

#Set grass region parameters to 'wwt' layer 
execGRASS('g.region',parameters = list(raster='wwt'))

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
execGRASS("r.watershed",parameters = list(elevation='dem',threshold=10,accumulation='acc',tci='topo_idx',spi='strm_pow',drainage='dir',stream='stream_r',length_slope='slop_lngth',slope_steepness='steep'),flags = c('overwrite','quiet','a'))

acc_r<-readRAST('acc',ignore.stderr = T)
acc_r<-stars::st_as_stars(acc_r)
st_crs(acc_r)<-bcalb_prj4



wetlnd_pour_pnt<-function(wetland)
{
  
  acc_filt<-as.data.frame(acc_r[wetland])
  
  acc_max<-acc_filt[which.max(acc_filt$acc),]
  
  acc_x<-acc_max$x
  
  acc_y<-acc_max$y
  
  return(c(acc_x,acc_y))
  
}

pour_pnts<-c()

for(feat in c(1:nrow(wwt_ws_v)))
{
  progress(feat, max.value = nrow(wwt_ws_v),progress.bar = TRUE)
  
  UID<-wwt_ws_v[feat,]$UID
  
  x_y<-wetlnd_pour_pnt(wwt_ws_v[feat,])
  
  if(feat==1)
  {
    pour_pnts<-c(x_y[1],x_y[2],UID)
  }else{
    pour_pnts<-rbind(pour_pnts,c(x_y[1],x_y[2],UID))
  }
  
}

pour_pnts<-as.data.frame(pour_pnts)
colnames(pour_pnts)<-c('x','y','UID')

####Example of summarizing basin statistics for wetland pour points for the binary cutblock(adult) layer, use 26 threads####

#Example of GRASSBasinStats use for binary cutblock(adult) layer 
test<-GRASSBasinStats::get_basin_stats(basin_df=pour_pnts,procs=26,stat_rast='cb_adult')

test$pa_cb_adult<-test$SUM/test$N

wwt_ws_v<-wwt_ws_v %>% left_join(test,by='UID')

mapview::mapview(wwt_ws_v,zcol='pa_cb_adult')


####!!!!Where from here, convert pixels to polygons, clumping ????!!!!####
###############################END OF TESTED PART OF SCRIPT###########################################################
######################################################################################################################


zonal(stack(r, r*10), z, 'sum')
##add a section to calculate v1 area for each clump and graph against elevation
##consider using MRVB
##consider making sub-watershed

ww.clump.poly <- rasterToPolygons(v1.clumps,dissolve = TRUE)


###?????????when should this be done in the order????????????

####make polygons of queens case clumped  wetland/water predictions
## upland set to NULL (NA) wetlands are 2 and water is 3
##add a field for the area of each polygon (will be in ha??)
ww.ws.poly <- wwt_ws %>% clump(direction = 8) %>% rasterToPolygons(dissolve = TRUE) %>% 
  smoothr::smooth(method = "ksmooth", smoothness = 2) %>% st_as_sf()
ww.ws.poly$area<- st_area(ww.ws.poly)



wwpolyclumps.ws <- v1.clumps %>% rasterToPolygons(dissolve = TRUE) %>% 
  smoothr::smooth(method = "ksmooth", smoothness = 2) %>% st_as_sf()
wwpolyclumps.ws$area<- st_area(wwpolyclumps.ws)
## TODO add the v1.clumps to the values.brick

################################################
###############explore the data#################
length(unique(sn$LINEAR_FEATURE_ID))#reports the number of unique inteters in the linear feature id field in the sn sf
plot(st_geometry(ws))#plot the geometry of ws
plot(sn["STREAM_ORDER"])#plot the stream network sf coloured by stream order

map1 <- mapview(sn.ms, maxpixels = 30672) # set maxpixels to the number of pixels in your raster
map2 <- mapview(ww.ms, maxpixels = 255486)
map3 <- mapview(sn.ms.r, maxpixels = 255486)
map4 <- map1+map2+map3
map4
###########################################################


sa.ws <- surfaceArea(as(dem.ws, 'SpatialGridDataFrame'))
area.ws <- st_area(ws)
AreaTosurfaceaAreaRatio <-  area.ws/sa.ws

##END TESTED AND EFFICIENT PART OF SCRIPT##############
###################### to this point is working efficiently############
## trying to make wetland clumps into polgyons
values (ww.ms)[values(ww.ms)<2] = NA
ww.ms.poly <- ww.ms %>% clump(direction = 8) %>% rasterToPolygons(dissolve = TRUE) %>% 
  smoothr::smooth(method = "ksmooth", smoothness = 2) %>% st_as_sf()
ww.ms.poly$area<- st_area(ww.ms.poly)

values (ww.ws)[values(ww.ws)<2] = NA
ww.ws.poly <- ww.ws %>% clump(direction = 8) %>% rasterToPolygons(dissolve = TRUE) %>% 
  smoothr::smooth(method = "ksmooth", smoothness = 2) %>% st_as_sf()
ww.ws.poly$area<- st_area(ww.ws.poly)

# work on adding the rasters together. 
## note the rasters do not line up?? why? wwt and bcdem should have the same origin. 
#unless i used the web version of wwt
sn.ms.r[is.na(sn.ms.r[])] <- 0
ww.ms[is.na(ww.ms[])] <- 0
r <- sn.ms.r + ww.ms
NAvalue(r) <- 0 #this doesnt seem to have worked.
mapview(r)

#Crop data to the ms for testing. this is a good method for line and point features
sn.ms <- st_intersection(sn, ms)
dem.ms <- dem.ws %>% raster::crop(ms) %>% raster::mask(ms)
sn.ms.r <- sn.ms %>% st_buffer(dist=20) %>% 
  st_union() %>% st_as_sf() %>% 
  fasterize(dem.ms)######lightning fast. YAYA!!
plot(sn.ms.r)


####################################################################

##############other stuff to work on#################

##make the bcdc read into a funtion with a map sheet and or a ws as input

##IMPROVEMENT: add in steps to make sure ww and dem.ws match origin, projection, pixel size, etc.
## consider stacking

##Read in area of interst aoi.file, transfrom to bcalb projection, union any internal polygons
##buffer with optional buffer distance
aoi <- read_sf(aoi.file,stringsAsFactors = TRUE) %>%
  st_transform(crs = bcalb) %>% 
  st_union() %>%
  st_buffer(dist = buffer) %>% 
  st_as_sf()
#plot(st_geometry(aoi))

dem.aoi <- mask(dem.bc, aoi)
