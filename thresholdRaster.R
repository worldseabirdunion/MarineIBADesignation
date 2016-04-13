## thresholdRaster     ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## thresholdRaster applies a threshold to a raster, and isolates any areas above
## that threshold value. Converting the raster values to polygons is difficult and
## so this part takes some time. The function returns a SpatialPolygonsDataFrame
## containing the polygons that are above threshold, and with an attributes table
## holding each sites Maximum raster value.

## CountRas must be a raster object with values 0 - 1.
## Threshold must be a number indicating the percentage value to be used
## as the threshold.


thresholdRaster <- function(CountRas, Threshold = 10)
{
  
  require(raster)
  require(maps)
  require(geosphere)
  
  plot(CountRas, asp=1)
  map("world", add=T, fill=T, col="darkolivegreen3")
  Threshold <- Threshold/100
  RasSites <- CountRas >= Threshold
  plot(RasSites, asp=1, col=rev(heat.colors(25)))
  map("world", add=T, fill=T, col="darkolivegreen3")
  
  if(length(which(getValues(CountRas) > Threshold)) < 1)
  {
    Mid <- c(bbox(CountRas)[1,1]+(as.numeric(bbox(CountRas)[1,2])- as.numeric(bbox(CountRas)[1,1]))/2,  midPoint(bbox(CountRas)[,2], bbox(CountRas)[,1])[2])
    text(Mid, "No Site Identified", cex=1.25)
    stop("No cells were above the threshold value")
  }
  
  Cells <- rasterToPolygons(CountRas, fun=function(x) {x>Threshold})
  DateLine <- Cells@bbox[1,1] < -178 & Cells@bbox[1,2] > 178
  if(DateLine == TRUE)    {Cells <- spTransform(Cells, CRS=DgProj)}
  
  Sites <- dissolve(Cells)
  ifelse(DateLine == TRUE, projection(Sites) <- DgProj, projection(Sites) <- "+proj=longlat + datum=wgs84")
  Sites <- spTransform(Sites, CRS=CRS("+proj=longlat + datum=wgs84"))
  SiteTable <- data.frame(SiteID = names(Sites), MaxPerc = round(extract(CountRas, Sites, fun=max)*100,2))
  Sites <- SpatialPolygonsDataFrame(Sites, data=SiteTable)
  print(SiteTable)
  return(Sites)
}

dissolve <- function(Cells)
{
  require(sp)
  require(rgeos)
  require(maptools)
  CellsAv <- Cells
  plot(Cells)
  j <- 0
  while(length(CellsAv) > 0)
  {
    j <- j + 1
    
    Cell1 <- CellsAv[1,]
    CellsAv <- CellsAv[-1,]
    if(length(CellsAv) < 1)
    {
      CellsMerge <- spChFIDs(Cell1, as.character(j))
      if(j == 1) {Sites <- CellsMerge} else
        Sites <- spRbind(Sites, CellsMerge)
      next
    }
    CellsNr <- which(gTouches(Cell1, CellsAv, byid=T))
    if(length(CellsNr) < 1)
    {
      CellsMerge <- spChFIDs(Cell1, as.character(j))
      if(j == 1) {Sites <- CellsMerge} else
        Sites <- spRbind(Sites, CellsMerge)
      next
    }
    CellsSel <- CellsAv[as.double(CellsNr),]
    CellsMerge <- gUnion(Cell1, CellsSel)
    CellsAv <- CellsAv[-as.double(CellsNr),]
    if(length(CellsAv) < 1)
    {
      CellsMerge <- spChFIDs(CellsMerge, as.character(j))
      if(j == 1) {Sites <- CellsMerge} else
        Sites <- spRbind(Sites, CellsMerge)
      next
    }
    CellsNr <- which(gTouches(CellsMerge, CellsAv, byid=T))
    while(length(CellsNr) > 0)
    {
      CellsSel <- CellsAv[as.double(CellsNr),]
      CellsMerge <- gUnion(CellsMerge, CellsSel)
      plot(CellsMerge, add=T, col=2)
      CellsAv <- CellsAv[-as.double(CellsNr),]
      if(length(CellsAv) < 1) break
      CellsNr <- which(gTouches(CellsMerge, CellsAv, byid=T))
    }
    CellsMerge <- spChFIDs(CellsMerge, as.character(j))
    if(j == 1) {Sites <- CellsMerge} else
      Sites <- spRbind(Sites, CellsMerge)
    plot(Sites, add=T, col=2)
  }
  plot(Sites, col=names(Sites), add=T)
  return(Sites)
}

