## function to fasterize bcdc stream network
##IMPROVEMENT: generalize the function to take any bcdcID or name. 
##Have to add an if statement for point and line features to add buffer.
##add bcdcID and buffer distance for point and line feateres to the input.
##or create two functions to fasterize lines and points with a buffer
## and a second to fasterize.bcdc.poly
#ms <- ms.ws %>% filter(MAP_TILE == mapsheet)
fasterize.sn <- function(reference.raster, mapsheet = "093J070", buf = (sum(res(raster.ref)/2))) {
  library (raster)
  library (fasterize)
  library (bcdata)
  
  ms <- bcdc_query_geodata("a61976ac-d8e8-4862-851e-d105227b6525", crs = bcalb) %>%
    filter( MAP_TILE == mapsheet) %>%
    collect()
  sn <- bcdc_query_geodata("92344413-8035-4c08-b996-65a9b3f62fca", crs = bcalb) %>%
    filter( INTERSECTS(ms)) %>%
    collect() %>% st_as_sf()
  sn.ms <- st_intersection(sn, ms)
  dem.ms <- dem.ws %>% raster::crop(ms) %>% raster::mask(ms)
  sn.ms.r <- sn.ms %>% st_buffer(dist=buf) %>% 
    st_union() %>% st_as_sf() %>% 
    fasterize(raster.ref)
  return(sn.ms.r)
}
sn.ms.r <- fasterize.sn(mapsheet, ww.ms, 25)
