rm(list=ls())
library(this.path)
setwd(this.path::here())
library(imager)
library(terra)
source("../Functions.R")
library(Rcpp)
sourceCpp("../excursion-sets/binary_fields.cpp")
source("../excursion-sets/excursions.R")

# load in the raster data, hourly rainfall in france for a 1km x 1km grid
files = "raster_exprank.tif"
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

# set the threshold values for the excursion sets
p.vals = seq(0.9,0.999,by=0.0025)

for (loc.idx in 1:3) {
  cat("Estimating perimeters for location",loc.idx,"\n")
  
  # target cell (row/col)
  cell0 <- terra::cellFromXY(r, crds(pts_proj[loc.idx, ]))
  rc0   <- terra::rowColFromCell(r, cell0)
  ix <- rc0[2]
  iy <- nrow(r) - rc0[1] + 1
  
  results <- vector("list", 0) #vector("list", dim(r)[3])  # or a growable list, appended with [[]]
  
  n.obs = dim(r)[3]
  for (obs.idx in 1:n.obs) {
    cat("\tProgress:",round(100*obs.idx / n.obs,2),"%\n")
    
    if (vals.at.s0[obs.idx, loc.idx] < qexp(p.vals[1])) next
    
    r.i <- r[[obs.idx]]
    r.i.mat <- terra::as.matrix(r.i, wide = TRUE)  # this is a bit slow
    
    for (p.val in p.vals) {
      u <- qexp(p.val)
      
      row <- data.frame(obs.idx=obs.idx, p=p.val,
                        perim_biased=NA_real_, perim_unbiased=NA_real_)
      
      if (vals.at.s0[obs.idx, loc.idx] < u) {
        row$perim_biased = 0
        row$perim_unbiased = 0
        results[[length(results) + 1]] <- row
        next
      }
      
      # compute the binary field and take the connected compontn touching s0
      bf <- (r.i.mat > u) * 1
      im <- as.cimg(t(bf[nrow(bf):1, ]))   
      lbl <- label(im, high_connectivity = FALSE)
      target_label <- lbl[ix, iy, 1, 1]
      bf.cc <- as.matrix(lbl == target_label) * 1
      
      # compute the perimeters in 2 ways
      row$perim_biased   <- perim.L1.na.rm(bf.cc)
      row$perim_unbiased <- perimBF(bf.cc, p = 2, incl_boundary = FALSE, box_size = "auto")

      results[[length(results) + 1]] <- row
    }
    res.df <- do.call(rbind, results)
    save(res.df, file = paste0("perims_loc", loc.idx, "_better_v2.RData"))
  }
}

stop("Finished.")