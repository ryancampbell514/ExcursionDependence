# This code loads in raster files from the Comephore dataset and standardises the margins

rm(list=ls())
library(stringr)
library(terra)
library(viridis)
library(this.path)
setwd(this.path::here())

terraOptions(memmax = 2)

# read the raster file locations
files <- list.files(
  "PATH_TO_DATA",
  pattern = "_RR\\.gtif$",
  full.names = TRUE,recursive=TRUE
)

# extract them as raster
r = rast(files)
print(dim(r))
ndays = dim(r)[3]
print("Raster dimensions:")
print(dim(r))

# set out of range as NA (rather than NaN)
print("Set out of range to NA")
r2 <- ifel(r == 65535, NA, r)

# set Corsica to NA
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
print("Rank transformation...")
r.exp <- app(r3, fun = function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  ok <- !is.na(x)
  out <- rep(NA_real_, length(x))
  out[ok] <- rank(x[ok], ties.method = "average") / (sum(ok) + 1)
  qexp(out)
})
writeRaster(r.exp,
            paste0("raster_exprank.tif"),
            datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9"),
            overwrite = TRUE)
