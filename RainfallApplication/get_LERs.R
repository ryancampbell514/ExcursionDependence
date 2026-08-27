rm(list=ls())

library(this.path)
setwd(this.path::here())
library(terra)
source("../Functions.R")
library(ExcursionSets)

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

# set the threshold values for the excursion sets
# x.axis = seq(log(-log(1-0.8)),log(-log(1-0.99)),length.out=30)
# p.vals = 1-exp(-exp(x.axis))[-c(1:6)]  # needed, because of gaps in data due to zero-inflation
p.vals = seq(0.95,0.999,by=0.0025)

# get seasonal information
start <- as.POSIXct("2024-08-01 00:00:00", tz = "UTC")
times <- start + (0:(nlyr(r) - 1)) * 3600
month <- as.integer(format(times, "%m"))
season <- ifelse(month %in% 4:9, 1, 0)

library(argparse)
parser <- ArgumentParser()
parser$add_argument("--loc", type = "integer", default=1, required = TRUE)
args <- parser$parse_args()
loc.idx <- args$loc

cond_site = vect(
  data.frame(lon=cities[loc.idx,2], lat=cities[loc.idx,3]),
  geom=c("lon","lat"),
  crs="EPSG:4326"
)
s0_proj <- project(cond_site, crs(r))

cat("Estimating perimeters for location",loc.idx,"\n")

results <- vector("list", 0) #vector("list", dim(r)[3])  # or a growable list, appended with [[]]

n.obs = dim(r)[3]
for (obs.idx in 1:n.obs) {
  cat("\tProgress:",100*obs.idx / n.obs,"%\n")
  
  if (vals.at.s0[obs.idx, loc.idx] < qexp(p.vals[1])) next
  r.i <- r[[obs.idx]]

  bf.all = r.i > qexp(p.vals)
  p.i = 1
  for (p.val in p.vals) {
    u <- qexp(p.val)
    if(vals.at.s0[obs.idx, loc.idx] < u){
      ler=0
    } else {
      bf = bf.all[[p.i]]  # this is slow
      ler=terra.ler(rf=NULL,bf=bf,thresh=NULL,half.width=NULL,s0=s0_proj)  # units=meters
    }
    row <- data.frame(obs.idx=obs.idx, season=season[obs.idx], p=p.val, ler=ler)
    results[[length(results) + 1]] <- row
    p.i = p.i+1
  }
  res.df <- do.call(rbind, results)
  save(res.df, file = paste0("LER_loc", loc.idx, "_better.RData"))
}

stop("Finished.")

#######################################################################

library(dplyr)
library(quantreg)
rm(list=ls())
# loc.idx = 1
for(loc.idx in 1:4){
  load(paste0("/scratch/rcampbell/nice/Comephore_data/LER_loc", loc.idx, "_better.RData"))
  # head(res.df)
  res.df = res.df %>% mutate(log.ler = log(ler), nl.nl.p = -log(-log(1-p))) %>% filter(is.finite(log.ler))
  fit1 = rq(log.ler~nl.nl.p, data=res.df)
  # fit2 = rq(log.ler~nl.nl.p+season, data=res.df)
  print(round(coefficients(fit1),2))
  # print(round(coefficients(fit2),2))
}
for(loc.idx in 1:4){
  load(paste0("/scratch/rcampbell/nice/Comephore_data/LER_loc", loc.idx, "_better.RData"))
  # head(res.df)
  res.df = res.df %>% mutate(log.ler = log(ler), nl.nl.p = -log(-log(1-p))) %>% filter(is.finite(log.ler))
  # fit1 = rq(log.ler~nl.nl.p, data=res.df)
  fit2 = rq(log.ler~nl.nl.p+season, data=res.df)
  # print(round(coefficients(fit1),2))
  print(round(coefficients(fit2),2))
}