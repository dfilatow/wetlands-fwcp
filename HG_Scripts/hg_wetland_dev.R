
#Fill sinks in a local SRTM DEM, name 'dem_fill'
RSAGA::rsaga.fill.sinks('/home/hunter/Downloads/dem.sgrd','/home/hunter/Downloads/dem_fill.sgrd',method="wang.liu.2006")


#Function return catchment scale statistics provided a DEM, x-y pour point coordinates (map units), grid(s) to derive stats from, and desired statistic(s). 
get_basin_stats<-function(dem_pth,x,y,agg_grids,stat)
{
  #Get up-slope catchment area for the point defined by 'x,y', derived from filled DEM using MFD algorithm in SAGA 
  RSAGA::rsaga.geoprocessor('ta_hydrology',4,param = list(TARGET_PT_X=x,TARGET_PT_Y=y,ELEVATION=dem_pth,AREA=paste(tempdir(),'/basin.sgrd',sep=""),METHOD=2))
  
  #Using output from up-slope catchment area, convert to binary grid where MFD contributing is > 50%
  RSAGA::rsaga.grid.calculus(in.grids = paste(tempdir(),'/basin.sgrd',sep=""), out.grid = paste(tempdir(),'/basin_bnry.sgrd',sep=""), formula = 'gt(a,50)')
  
  #Calculate zonal statistics using binary output grid as zones (1 corresponds to up-slope catchment), save stats to CSV
  RSAGA::rsaga.geoprocessor('statistics_grid',5,param = list(ZONES=paste(tempdir(),'/basin_bnry.sgrd',sep=""),STATLIST=agg_grids,OUTTAB=paste(tempdir(),'/basin_stats.csv',sep="")))
  
  #Read in zonal statistics results, return desired statistics. 
  rslt_tab<-read.csv(paste(tempdir(),'/basin_stats.csv',sep=""))
  if(stat=="tabl")
  {
    return(rslt_tab)
  }
  if(stat=="count")
  {
    return(rslt_tab[2,3])
  }
  if(stat=="min")
  {
    return(rslt_tab[2,4])
  }
  if(stat=="max")
  {
    return(rslt_tab[2,5])
  }
  if(stat=="mean")
  {
    return(rslt_tab[2,6])
  }
  if(stat=="sd")
  {
    return(rslt_tab[2,7])
  }
  if(stat=="sum")
  {
    return(rslt_tab[2,8])
  }
}

saga_create_basin<-function(dem_pth,x,y,UID,out_dir,envr)
{
  #Get up-slope catchment area for the point defined by 'x,y', derived from filled DEM using MFD algorithm in SAGA 
  RSAGA::rsaga.geoprocessor('ta_hydrology',4,param = list(TARGET_PT_X=x,TARGET_PT_Y=y,ELEVATION=dem_pth,AREA=paste(tempdir(),'/',UID,'_basin.sgrd',sep=""),METHOD=2),env=envr)
  
  #Using output from up-slope catchment area, convert to binary grid where MFD contributing is > 50%
  RSAGA::rsaga.grid.calculus(in.grids = paste(tempdir(),'/',UID,'_basin.sgrd',sep=""), out.grid = paste(tempdir(),'/',UID,'_basin_bnry.sgrd',sep=""), formula = 'gt(a,50)',env = envr)
  
  #Convert up-slope binary raster to vector. 
  RSAGA::rsaga.geoprocessor('shapes_grid',6,param = list(GRID=paste(tempdir(),'/',UID,'_basin_bnry.sgrd',sep=""),POLYGONS=paste(out_dir,'/',UID,'_basin.shp',sep=""),CLASS_ALL=0,SPLIT=0),env = envr)
  
}



#An example, for basin defined by point (EPSG:3005; x:1299974,y:954998)
start<-Sys.time()
hm<-get_basin_stats('/home/hunter/Downloads/dem_fill.sgrd',1299974,954998,'/home/hunter/Downloads/dem.sgrd',"tabl")
end<-Sys.time()


hmm<-saga_create_basin('/home/hunter/Downloads/dem_fill.sgrd',1299974,954998,103,'/home/hunter/Downloads')

