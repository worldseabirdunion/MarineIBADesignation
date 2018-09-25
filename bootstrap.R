## bootstrap     ######################################################################################################

## Phil Taylor & Mark Miller, 2012

## this script Iteratively subsamples a dataset of tracking data, investigating the
## affect of sample size. At each iteration the data is split, one half is used as the
## 'training' data and the 50%UD is calculated from this. The second half is used as
## 'testing' and the proportion of it, captured within the 50%UD is calculated.
## A perfect dataset would tend towards 0.5.
## By fitting a trend line to this relationship we can establish the sample size at which
## any new data would simply add to the existing knowledge. This script indicates how close to
## this value the inputted data are.

## DataGroup must be a dataframe or SpatialPointsDataFrame with Latitude, Longitude and ID as fields.
## Scale determines the smoothing factor used in the kernel analysis.
## Iteration determines the number of times each sample size is iterated.


bootstrap <- function(DataGroup, Scale=100, Iteration=50)
{
  
  require(sp)
  require(geosphere)
  require(rgdal)
  require(adehabitatHR)
  require(foreach)
  require(doParallel)
  require(parallel)
  
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
  BoundBox <- bbox(DataGroup)
  UIDs <- unique(DataGroup$ID)
  Ntrips <- length(UIDs)
  Nloop<- seq(1,(Ntrips-1),ifelse(Ntrips>100,10,1))
  DoubleLoop <- data.frame(SampleSize = rep(Nloop,each=Iteration), Iteration=rep(seq(1:Iteration),length(Nloop)))
  LoopNr <- seq(1:dim(DoubleLoop)[1])	
  UDLev <- 50
  
  #setup parallel backend to use 4 processors
  cl<-makeCluster(detectCores())
  registerDoParallel(cl)
  Result<-data.frame()
  
  Result <- foreach(LoopN=LoopNr, .combine = rbind, .packages=c("sp","adehabitat","geosphere","rgdal")) %dopar% {
    
    N<-DoubleLoop$SampleSize[LoopN]
    i<-DoubleLoop$Iteration[LoopN]
    Coverage <- NULL
    Inclusion <- NULL
    History <- NULL
    
    Output <- data.frame(SampleSize = N, InclusionMean = 0,Iteration=i)
    
    RanNum <- sample(UIDs, N, replace=F)
    SelectedCoords <- as.data.frame(coordinates(DataGroup[DataGroup$ID %in% RanNum,]))
    NotSelected <- DataGroup[!DataGroup$ID %in% RanNum,]
    Temp <- data.frame(SelectedCoords[,1], SelectedCoords[,2])
    Ext <- (min(Temp[,1]) + 3 * diff(range(Temp[,1])))
    if(Ext < (Scale * 1000 * 2)) {BExt <- ceiling((Scale * 1000 * 3)/(diff(range(Temp[,1]))))} else {BExt <- 3}
    coordinates(SelectedCoords) <- ~DataGroup.Longitude + DataGroup.Latitude
    KDE.Surface <- kernelUD(SelectedCoords, h=Scale*1000, grid=70, extent=BExt,  same4all=FALSE)
    KDE.Spl <- getverticeshr(KDE.Surface, lev = UDLev)
    KDE.Spl@proj4string <- DgProj
    Overlain <- over(NotSelected, KDE.Spl)
    Output$InclusionMean <- length(which(Overlain == 1))/nrow(NotSelected)
    return(Output)
  }
  
  ## stop the cluster
  stopCluster(cl)
  
  par(mfrow=c(1,1), mai=c(1,1,1,1))
  #Result <- Output[1:nrow(Output)-1,]
  write.table(Result,"bootout_temp.csv", row.names=F, sep=",")
  M1 <- nls((Result$InclusionMean ~ (a*Result$SampleSize)/(1+b*Result$SampleSize)), data=Result, start=list(a=1,b=0.1))
  PredData <- data.frame(SampleSize = unique(Result$SampleSize))
  Result$pred<-predict(M1)
  P2 <- aggregate(pred~SampleSize, Result, FUN=mean)
  P2$sd <- aggregate(InclusionMean~SampleSize, Result, FUN=sd)[,2]
  plot(InclusionMean~SampleSize, data=Result, pch=16, cex=0.2, col="darkgray", ylim=c(0,1), xlim=c(0,nrow(PredData)), ylab="Inclusion", xlab="SampleSize")
  yTemp <- c((P2[,2] + P2[,3]), rev(P2[,2] - P2[,3]))
  xTemp <- c(1:nrow(P2), nrow(P2):1)
  polygon(x=xTemp, y=yTemp, col="gray93", border=F)
  points(InclusionMean~SampleSize, data=Result, pch=16, cex=0.2, col="darkgray")
  lines(P2, lty=1,lwd=2)
  Asymptote <- (summary(M1)$coefficients[1]/summary(M1)$coefficients[2])
  RepresentativeValue <- max(P2$pred)/Asymptote*100
  print(RepresentativeValue)
  text(x=0, y=1,paste(round(RepresentativeValue,2), "%", sep=""), cex=2, col="gray45", adj=0)
  Result$RepresentativeValue <- RepresentativeValue
  return(Result)
}

