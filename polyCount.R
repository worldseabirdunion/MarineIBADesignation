## polyCount  ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## polyCount calculates the number of overlapping Polygons. The calculation is
## done by overlapping the polygons within a grid and counting the number
## falling within each gridcell. The Res parameter sets the size of the grid
## (in decimal degrees) to be used. Smaller grids will allow for more detailed
## counts but will slow computing time. The function returns a raster object with
## values showing the proportion of Polys (i.e. nPolys/total nPolys) overlapping
## each cell. The raster extent is set at the bounding limits of the Polys or, when
## data crosses the dateline, set to the northern and southern most limits of the
## Polys but longitudinally crossing the circumference of the world.

## Polys must be a SpatialPolygonsDataFrame of the polygons to be counted.
## Res must be a numeric object indicating the resolution in decimal degrees.

## version 1.2    05-04-2012

polyCount <- function(Polys, Res = 0.1)
{
  
  require(raster)
  require(maps)
  
  if(!class(Polys) %in% c("SpatialPolygonsDataFrame", "SpatialPolygons")) stop("Polys must be a SpatialPolygonsDataFrame")
  if(is.na(projection(Polys))) stop("Polys must be projected")
  
  Poly.Spdf <- spTransform(Polys, CRS=CRS("+proj=longlat +ellps=WGS84"))
  DgProj <- Polys@proj4string
  
  DateLine <- Poly.Spdf@bbox[1,1] < -178 & Poly.Spdf@bbox[1,2] > 178
  if(DateLine == TRUE) {print("Data crosses DateLine")}
  
  UDbbox <- bbox(Poly.Spdf)
  if(DateLine == TRUE)  {UDbbox[1,] <- c(-180,180)}
  BL <- floor(UDbbox[,1]) + (Res/2)
  TR <- ceiling(UDbbox[,2])
  NRow <- ceiling(sqrt((BL[1] - TR[1])^2)/Res)
  NCol <- ceiling(sqrt((BL[2] - TR[2])^2)/Res) + (Res * 100)
  Grid <- GridTopology(BL, c(Res,Res), c(NRow, NCol))
  SpGrid <- SpatialPoints(Grid, proj4string = CRS("+proj=longlat + datum=wgs84"))
  SpdfGrid <- SpatialPointsDataFrame(SpGrid, data.frame(Longitude=SpGrid@coords[,1], Latitude=SpGrid@coords[,2]))
  SpGridProj <- spTransform(SpdfGrid, CRS=DgProj)
  GridIntersects <- over(SpGridProj, Polys)
  SpGridProj@data$Intersects$ID <- GridIntersects$ID
  SpGridProj <- SpGridProj[!is.na(SpGridProj@data$Intersects$ID),] 			###subset(SpGridProj, !is.na(SpGridProj@data$Intersects$ID))
  plot(SpGridProj)
  
  Count <- 0
  for(i in 1:length(Polys))
  {
    TempB <- Polys[i,]
    Temp <- over(SpGridProj, TempB)
    Temp[is.na(Temp)] <- 0
    Temp[Temp > 0] <- 1
    Count <- Count + Temp
    Prop <- Count/i
  }
  GridIntersects$inside<-as.numeric(as.character(GridIntersects$ID))
  GridIntersects$Prop <- 0
  GridIntersects$Prop[!is.na(GridIntersects$inside)] <- Prop[,1]
  SpdfGrid$Count <- 0
  SpGridVals <- SpatialPixelsDataFrame(SpGrid, data.frame(Values = GridIntersects$Prop))
  SGExtent <- extent(SpdfGrid)
  RT <- raster(SGExtent, as.double(NCol), as.double(NRow))
  WgsRas <- (rasterize(x=SpGridVals, y=RT, field = "Values"))
  
  plot(WgsRas, asp=1)
  map("world", add=T, fill=T, col="darkolivegreen3")
  projection(WgsRas) <- CRS("+proj=longlat + datum=wgs84")
  return(WgsRas)
}