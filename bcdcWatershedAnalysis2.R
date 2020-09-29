##WWA.R Wetland Watershed Analysis
##author : Deepa Filatow
##initiated: March 12, 2020
##last edited: May 4, 2020

##option to clear memeory
#rm(list = ls())

############################################################################
##Load libraries############################################################
############################################################################

library(bcdata)##https://github.com/bcgov/bcdata
library(sf)
library(tidyverse)
library(mapview)
library(raster)
library(fasterize)
library(rdgal)
library(igraph) #was required by raster::clump with directions=8. got an error and had to install
library(rgeos) #got an error trying to crop and rasterToPolygon
#install.packages("devtools")
#devtools::install_github("mstrimas/smoothr") ##tried to load from CRAN without success
library(smoothr)


############################################################################
### ADD TEST VARIABLES #####################################################
############################################################################

#BC Alb. EPSG Code
bcalb <- 3005

#Path to 3 class wetland classification output
wwt.file <- "Data/Categ/3CategoryPrediction/3CategoryPrediction/20190926-103025_map_recl.tif" #water(3), wetland(2), terrestrial prediction surface(1)

# A test watershed (ws) code...(PARS ws has plenty of data)
ws.code <- "PARS"

#This is the mapsheet code within the PARS test watershed
mapsheet <- "093J070"

#Get the mapsheet from bcdc using above code
ms <- bcdc_query_geodata("a61976ac-d8e8-4862-851e-d105227b6525", crs = bcalb) %>%
  filter( MAP_TILE == mapsheet) %>%
  collect()

#Load the output from the 2020 FWCP-Peace wetland model as a raster
wwt <- raster(wwt.file)

#Define a buffer distance 
BUFF <-  100

############################################################################
### Watershed Group ########################################################
############################################################################


#Create sf polygon: watershed group, ws from bcdc watershed code = ws.code
#https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-watershed-groups
bcdc_describe_feature("51f20b1a-ab75-42de-809d-bf415a0f9c62")

ws <- bcdc_query_geodata("51f20b1a-ab75-42de-809d-bf415a0f9c62", crs = bcalb) %>%
  filter(WATERSHED_GROUP_CODE == ws.code) %>%
  collect() %>% st_as_sf()


##TO DO:Test if conversion to sf is necessary for rasterization below now that raster:: is used

ws.area_m2 <- sum(ws$AREA_HA)*10000 # watershed area for in hectrars later calculations. have to use sum in case of mulipolygons on the coast

#plot(st_geometry(ws))

##Create sf polygon: 20K mapsheet grids that intersect with the target watershed group
##https://catalogue.data.gov.bc.ca/dataset/bcgs-1-20-000-grid
ms.ws <- bcdc_query_geodata("a61976ac-d8e8-4862-851e-d105227b6525", crs = bcalb) %>%
  filter(INTERSECTS(ws)) %>%
  collect()

#plot(st_geometry(ms.ws))


#crop dem.habc to the watershed ws
dem.ws <- GrabCDED::cded_get(ws)#loadem
#NOTE:crop then mask is more efficient and applies the correct extent
#then directly masking from all of bc


############################################################################
### water and wetland values################################################
############################################################################


## crop and mask upland wetland water wwt raster (1 2 3 ) for the watershed (ws) and the mapsheet (ms)

wwt.ws <- wwt %>% raster::crop(ws) %>% raster::mask(ws)#wetland prediction clipped to watershed

wwt.count.df <- wwt.ws %>% freq() %>% as.data.frame() %>% na.omit()#data frame of non NA wetland catagory counts

wwt.ws.dem <- wwt.ws %>% projectRaster(dem.ws, bilinear) #Reproject the raster to the proveded DEM

ww.ws <- wwt.ws.dem

ww.ws<-raster::mask(ww.ws,ww.ws<2,maskvalue=TRUE)


##Create sf lines: stream network, sn for watershed code ws.code
#https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-stream-network
##warning: this takes 1-2 minutes depending on your choice of watershed
bcdc_describe_feature("92344413-8035-4c08-b996-65a9b3f62fca")
sn <- bcdc_query_geodata("92344413-8035-4c08-b996-65a9b3f62fca", crs = bcalb) %>%
  filter(WATERSHED_GROUP_CODE == ws.code) %>%
  collect() %>% st_as_sf()

#plot(sn["STREAM_ORDER"])

sn.length_m <- sum(sn$FEATURE_LENGTH_M)##length of streams in m in the watershed


############################################################################
### disturbance and land use layers#########################################
############################################################################


##Create sf polygon: cutblocks, cb from bcdc that intersects with ws
##https://catalogue.data.gov.bc.ca/dataset/harvested-areas-of-bc-consolidated-cutblocks-
bcdc_describe_feature("b1b647a6-f271-42e0-9cd0-89ec24bce9f7")
cb <- bcdc_query_geodata("b1b647a6-f271-42e0-9cd0-89ec24bce9f7", crs = bcalb) %>%
  filter(INTERSECTS(ws)) %>%
  collect() %>% st_as_sf()## Test if conversion to sf is necessary for rasterization below now that raster:: is used

#plot(st_geometry(cb))

cb.area_m2 <- sum(cb$FEATURE_AREA_SQM)

##Create sf lines: roads, rd from bcdc that intersects with ws
#https://catalogue.data.gov.bc.ca/dataset/digital-road-atlas-dra-master-partially-attributed-roads
##warning: this takes 1-2 minutes depending on your choice of watershed
bcdc_describe_feature("bb060417-b6e6-4548-b837-f9060d94743e")
rd <- bcdc_query_geodata("bb060417-b6e6-4548-b837-f9060d94743e", crs = bcalb) %>%
  filter(INTERSECTS(ws)) %>%
  collect() %>% st_as_sf()

##what about planned roads
##https://catalogue.data.gov.bc.ca/dataset/bcts-planned-roads

##Add in a part to summarize m of road
rd.length_m <- sum(rd$FEATURE_LENGTH_M)

#fasterize disturbance layers
cb.r <- fasterize(cb, dem.ws, background = 0)
rd.r <- rd %>% st_buffer(dist=(sum(res(dem.ws)/2))) %>% fasterize(dem.ws, background = 0)
sn.r <- sn %>% st_buffer(dist=(sum(res(dem.ws)/2))) %>% fasterize(dem.ws, background = 0)

#create a raster of connected wetland and stream rasters = 1 rest null
ww.sn.r <- sn.r + wwt.ws.dem
values(ww.sn.r)[values(ww.sn.r) <= 1] = NA
values(ww.sn.r)[values(ww.sn.r) > 1] = 1

landuse.stack <- stack(cb.r, rd.r)#create a brick of land use rasters lu1
hydrology.stack <- stack(ww.ws, sn.r, ww.sn.r)#create a brick of values v1. in this case wetlands, water and streams
count.landuse <- sum(landuse.stack)# count of landuses 0 for none. Max of n for where all layers overlap
count.hydrology <- sum(hydrology.stack)

# create a hydrology values layer hv where all cells containing values are 1 and all cells withouth values are null
hv <- count.hydrology
values(hv)[values(hv) <= 1] = NA
values(hv)[values(hv) >1 ] = 1
hv.clumps <- raster::clump(hv, directions=8)

##Add connected wetlands and streams v1.clumps to the hydrology raster stack
hydrology.stack <- addLayer(hydrology.stack, hv.clumps)


##try some stuff with v1.clumps

hist(log(freq(hv.clumps))) # where i have used dem = habc elev this is the hystogram of area (ha) of connected wetland clumps

hv.clump.elev <- zonal(z = hv.clumps, x = dem.ws, fun = 'sum', na.rm = TRUE)  

hv.clump.elev





zonal(stack(r, r*10), z, 'sum')
##add a section to calculate v1 area for each clump and graph against elevation
##consider using MRVB
##consider making sub-watershed

ww.clump.poly <- rasterToPolygons(v1.clumps,dissolve = TRUE)


###?????????when should this be done in the order????????????

####make polygons of queens case clumped  wetland/water predictions
## upland set to NULL (NA) wetlands are 2 and water is 3
##add a field for the area of each polygon (will be in ha??)
ww.ws.poly <- ww.ws %>% clump(direction = 8) %>% rasterToPolygons(dissolve = TRUE) %>% 
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
##todo: add steps to alline all rasters to the crs provided. Report an error where crs is missing.
sn.ms.r[is.na(sn.ms.r[])] <- 0
ww.ms[is.na(ww.ms[])] <- 0
r <- sn.ms.r + ww.ms
NAvalue(r) <- 0 #this doesn't seem to have worked.
mapview(r)

#Crop data to the ms for testing. this is a good method for line and point features
sn.ms <- st_intersection(sn, ms)
dem.ms <- dem.ws %>% raster::crop(ms) %>% raster::mask(ms)
sn.ms.r <- sn.ms %>% st_buffer(dist=20) %>% 
  st_union() %>% st_as_sf() %>% 
  fasterize(dem.ms)######lightning fast. YAYA!!
plot(sn.ms.r)

###try reading layers from file or dataframe
bcdc_layerlist_df <- read.csv2("C:/Users/dfilatow/Documents/R/WetlandsWatersheds/bcdc_layerlist_df2.csv", sep = ",")

for (i in 1:dim(bcdc_layerlist_df)[1]){
  assign(bcdc_layerlist_df$varname[i], bcdc_query_geodata(bcdc_layerlist_df$permalink[i], crs = 3005) %>% collect())
}#this works

##attempt to add the filter
#these work
assign(bcdc_layerlist_df$varname[2], bcdc_query_geodata(bcdc_layerlist_df$permalink[2], crs = 3005) %>% filter(WATERSHED_GROUP_CODE == "PARS") %>% collect())
assign(bcdc_layerlist_df$varname[1], bcdc_query_geodata(bcdc_layerlist_df$permalink[1], crs = 3005) %>% filter(MAP_TILE == "093J070") %>% collect())

##this does not work!!
for (i in 1:dim(bcdc_layerlist_df)[1]){
  assign(bcdc_layerlist_df$varname[i], bcdc_query_geodata(bcdc_layerlist_df$permalink[i], crs = 3005) %>% filter(bquote(bcdc_layerlist_df$filter[i])) %>% collect())
}
####################################################################

##############other stuff to work on#################

##make the bcdc read into a function with a map sheet and or a ws as input

##IMPROVEMENT: add in steps to make sure ww and dem.ws match origin, projection, pixel size, etc.
## consider stacking

##Read in area of interest aoi.file, transform to bcalb projection, union any internal polygons
##buffer with optional buffer distance
aoi <- read_sf(aoi.file,stringsAsFactors = TRUE) %>%
  st_transform(crs = bcalb) %>% 
  st_union() %>%
  st_buffer(dist = buffer) %>% 
  st_as_sf()
#plot(st_geometry(aoi))

dem.aoi <- mask(dem.bc, aoi)


