rm(list=ls())
library(stringr)
library(terra)
library(viridis)
library(this.path)
setwd(this.path::here())

source("../inla_modelling/functions-lifting.R")

# Process the data by month (otherwise get memory errors)
files.month = list.files(file.path(getwd(),"data_month"),
                         full.names = TRUE)
for(file.month in files.month){
  
  # get the name (year+month)
  nm = substr(file.month, nchar(file.month) - 5, nchar(file.month))
  
  cat("Prepocessing month",nm,"\n")
  
  # get all *_RR.gtif for this month (hourly)
  files  = list.files(
    path = file.month,
    pattern = "_RR.gtif",
    recursive = TRUE,
    full.names = TRUE
  )
  
  # extract them as raster
  r = rast(files)
  ndays = dim(r)[3]

  # set out of range as NA
  r2 <- ifel(r == 65535, NA, r)
  
  # set corsica to NA
  m <- r2[[1]]
  p <- patches(m) # Identify connected patches
  f <- freq(p)    # Find largest patch (mainland)
  main_patch <- f[which.max(f$count), "value"]
  mask_mainland <- p == main_patch
  mask_mainland[mask_mainland == 0] <- NA
  r3 <- mask(r2, mask_mainland)
  
  # convert to Exp(1) margins
  r.exp = app(r3, fun = function(x) {
    if(all(is.na(x))){
      return(rep(NA,ndays))
    }
    ok <- !is.na(x)
    out <- rep(NA, length(x))
    
    # rank among valid values
    # out[ok] <- rank(x[ok], ties.method = "min") / (sum(ok)+1)
    # out[ok] = dplyr::dense_rank(x[ok]) / (sum(ok)+1)
    # out[ok] = get.uniform(x[ok],p=0.05)
    p0 <- mean(x[ok] == 0)
    pos <- x[ok][x[ok] > 0]
    # Fpos <- ecdf(pos)
    u <- numeric(length(x[ok]))
    u[x[ok] == 0] <- runif(sum(x[ok] == 0), 0, p0)
    u[x[ok] > 0] <- p0 + (1 - p0) * (rank(pos, ties.method="average")/(length(pos)+1))#Fpos(x[x > 0])
    out[ok] = qexp(u)
    return(out)
  })
  
  writeRaster(r.exp,
              paste0("exprank_data/raster_exprank_",nm,".tif"),
              datatype = "FLT4S",
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9"),
              overwrite = TRUE)
}

#############################################################

# check the exponential transformation
# load in the raster data, hourly rainfall in france for a 1km x 1km grid
files <- list.files(
  "exprank_data",
  pattern = "\\.tif$",
  full.names = TRUE
)
r <- rast(files)
dim(r)
object.size(r)

# define the conditioning sites
cities <- data.frame(
  name = c("paris", "bourges", "valence", "nice"),
  lon  = c(2.3514, 2.3980, 4.8967, 7.2620),
  lat  = c(48.8575, 47.0812, 44.9337, 43.7102)
)
pts <- terra::vect(cities, geom = c("lon", "lat"), crs = "EPSG:4326")
pts_proj <- terra::project(pts, crs(r))
cells <- terra::cellFromXY(r, crds(pts_proj))
vals.at.s0 <- t(terra::extract(r, pts_proj))
vals.at.s0 = vals.at.s0[-1,]

par(pty="m",mfrow=c(2,2))
hist(vals.at.s0[,1],freq=F,breaks=30,ylim=c(0,0.2))
lines(seq(0,9,length.out=100),dexp(seq(0,9,length.out=100)))
hist(vals.at.s0[,2],freq=F,breaks=30,ylim=c(0,0.2))
lines(seq(0,9,length.out=100),dexp(seq(0,9,length.out=100)))
hist(vals.at.s0[,3],freq=F,breaks=30,ylim=c(0,0.2))
lines(seq(0,9,length.out=100),dexp(seq(0,9,length.out=100)))
hist(vals.at.s0[,4],freq=F,breaks=30,ylim=c(0,0.2))
lines(seq(0,9,length.out=100),dexp(seq(0,9,length.out=100)))
