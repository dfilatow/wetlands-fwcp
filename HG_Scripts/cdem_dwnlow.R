library(sf)
library(ggplot2)
library(bcmaps)


cdem_tiles<-sf::read_sf("/home/hunter/Downloads/cdem_index_250k.kmz")


bc<-bcmaps::bc_bound()

iter<-sf::st_intersection(cdem_tiles,bc)

sf::sf_project(cdem_tiles,bc)
