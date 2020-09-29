library(bcdata)
library(sf)
library(rlang)
library(dplyr)
library(GetBCDataBatch)

e_sf<-sf::st_read('/home/hunter/Downloads/BowronAOI.geojson') %>%
  st_as_sf()

e_sf<-st_transform(e_sf,3005)

r1<-c("92344413-8035-4c08-b996-65a9b3f62fca",NA,"WATERSHED_GROUP_CODE == 'PARS'")
r2<-c("cb1e3aba-d3fe-4de1-a2d4-b8b6650fb1f6","e_sf",NA)

layers<-as.data.frame(rbind(r1,r2))

colnames(layers)<-c('code','geom_ext','filter_exp')



hmm<-get_bcdata_batch(layers)
