rm(list=ls())
library(stringr)
library(terra)
library(viridis)
library(this.path)
setwd(this.path::here())

terraOptions(memmax = 2)

# source("../inla_modelling/functions-lifting.R")

# # Process the data by month (otherwise get memory errors)
# files.month = list.files(file.path(getwd(),"data_month"),
#                          full.names = TRUE)
# 
# dir <- "data_month"

files <- list.files(
  "data_month",
  pattern = "_RR\\.gtif$",
  full.names = TRUE,recursive=TRUE
)

# rasters <- lapply(files, rast)
# 
# # get all *_RR.gtif for this month (hourly)
# files  = list.files(
#   path = file.month,
#   pattern = "_RR.gtif",
#   recursive = TRUE,
#   full.names = TRUE
# )

# extract them as raster
r = rast(files)
print(dim(r))
ndays = dim(r)[3]
print("Raster dimensions:")
print(dim(r))

# set out of range as NA
print("Set out of range to NA")
r2 <- ifel(r == 65535, NA, r)

# set corsica to NA
print("Setting Corsica to NA")
m <- r2[[1]]
p <- patches(m) # Identify connected patches
f <- freq(p)    # Find largest patch (mainland)
main_patch <- f[which.max(f$count), "value"]
mask_mainland <- p == main_patch
mask_mainland[mask_mainland == 0] <- NA
r3 <- mask(r2, mask_mainland)

print("r3 memory info")
print(mem_info(r3, print = TRUE))

# convert to Exp(1) margins via rank transform
# r.exp = app(r3, fun = function(x) {
#   if(all(is.na(x))){
#     return(rep(NA,ndays))
#   }
#   ok <- !is.na(x)
#   out <- rep(NA, length(x))
#   
#   out[ok] = rank(x[ok], ties.method="average")/(length(x[ok])+1)
#   out = qexp(out)
#   return(out)
# })
# writeRaster(r.exp,
#             paste0("exp_data_v2/raster_exprank.tif"),
#             datatype = "FLT4S",
#             gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9"),
#             overwrite = TRUE)

print("Rank transformation...")
r.exp <- app(r3, fun = function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  ok <- !is.na(x)
  out <- rep(NA_real_, length(x))
  out[ok] <- rank(x[ok], ties.method = "average") / (sum(ok) + 1)
  qexp(out)
}, 
filename = "exp_data_v2/raster_exprank.tif", overwrite = TRUE, wopt = list(datatype = "FLT4S"))

# convert to Exp(1) margins via TO smooth method
print("Smooth GPD transformation...")
source("../inla_modelling/functions-lifting.R")
r.exp2 = app(r3, fun = function(x) {
  if(all(is.na(x))){
    return(rep(NA,ndays))
  }
  ok <- !is.na(x)
  out <- rep(NA, length(x))
  
  out[ok] = get.uniform(x[ok],0.9)
  out = qexp(out)
  return(out)
})

writeRaster(r.exp2,
            paste0("exp_data_v2/raster_expsmoothGP.tif"),
            datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9"),
            overwrite = TRUE)
