## batchUD  #####################################################################################################################

## Phil Taylor & Mark Miller, 2011

## batchUD calculates the Utilisation Distribution for groups of data in
## a larger dataset. The function relies on Adehabitat package for the
## calculations but returns the UD as SpatialPolygons so they can be exported
## as shapefiles and be more versatile in the script.

## DataGroup must be either a DataFrame or SpatialPointsDataFrame with Latitude,
## Longitude and ID as fields. UD will be calculated for each unique ID value and
## each row should be a location. For each UD to be comparable, the data should be
## regularly sampled or interpolated.
## Scale should be the smoothing factor to be used in the Kernel Density Estimation
## and should be provided in Km.
## UDLev should be the quantile to be used for the Utilisation Distribution.


batchUD <- function(DataGroup, Scale = 50, UDLev = 50)
{
  require(sp)
  require(maptools)
  require(rgdal)
  require(adehabitatHR)
  require(geosphere)
  
  if(!"Latitude" %in% names(DataGroup)) stop("Latitude field does not exist")
  if(!"Longitude" %in% names(DataGroup)) stop("Longitude field does not exist")
  if(!"ID" %in% names(DataGroup)) stop("ID field does not exist")
  
  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    mid_point<-data.frame(centroid(cbind(DataGroup$Longitude, DataGroup$Latitude)))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup$Longitude, DataGroup$Latitude), proj4string=CRS("+proj=longlat + datum=wgs84"))
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
  }else{DgProj<-DataGroup@proj4string}
  
  DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  
  UIDs <- unique(DataGroup$ID)
  note<-0
  KDE.Sp <- NULL
  for(i in 1:length(UIDs))
  {
    Trip <- DataGroup[DataGroup$ID == UIDs[i],]
    if(nrow(Trip@data)<6)
    {
      print(paste("ID =", UIDs[i], "has fewer than 6 points, too small to fit kernel. It will be excluded"))
      if(i == 1){note <- 1}
      next
    }
    TripCoords <- data.frame(Trip@coords)
    Temp <- data.frame(TripCoords[,1], TripCoords[,2])
    Ext <- (min(Temp[,1]) + 3 * diff(range(Temp[,1])))
    if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(Temp[,1]))))} else {BExt <- 3}
    KDE.Surface <- kernelUD(data.frame(TripCoords[,1], TripCoords[,2]), id=Trip$ID, h=(Scale * 1000), grid=100, extent=BExt, same4all=FALSE)
    KDE.UD <- getverticeshr(KDE.Surface, lev = UDLev)
    KDE.Sp1 <- kver2spol(KDE.UD)
    
    if(i==1 | note==1) {KDE.Sp <- KDE.Sp1} else
      KDE.Sp <- spRbind(KDE.Sp, KDE.Sp1)
    plot(KDE.Sp)
    note<-0
    if(i < length(UIDs)) {legend("bottomleft", paste(UIDs[i+1]))}
  }
  
  UIDs <- names(which(table(DataGroup$ID)>5))
  KDE.Sp@proj4string <- DgProj
  KDE.Wgs <- spTransform(KDE.Sp, CRS=CRS("+proj=longlat +ellps=WGS84"))
  Tbl <- data.frame(Name_0 = rep(1, length(UIDs)), Name_1 = 1:length(UIDs), ID = UIDs)
  row.names(Tbl) <- UIDs
  KDE.Spdf <- SpatialPolygonsDataFrame(KDE.Sp, data=Tbl)
  
  plot(KDE.Spdf, border=factor(UIDs))
  return(KDE.Spdf)
}
