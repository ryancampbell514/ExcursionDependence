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

files = "/scratch/rcampbell/nice/Comephore_data/exp_data_v2/raster_exprank.tif"
# files = "/scratch/rcampbell/nice/Comephore_data/exp_data_v2/raster_expsmoothGP.tif"
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

# # Check if the margins are Exp(1)  (in the right tails)
# par(mfrow=c(2,2))
# for(i in 1:ncol(vals.at.s0)){
#   hist(vals.at.s0[,i],freq=F)
#   lines(seq(0,30,by=0.1),dexp(seq(0,30,by=0.1)))
# }

# set the threshold values for the excursion sets
# p.vals1 = seq(0.7,0.95,by=0.05)
# p.vals2 = seq(0.99,0.999,by=0.001)
# p.vals = c(p.vals1,p.vals2)
# x.axis = seq(log(-log(1-0.8)),log(-log(1-0.99)),length.out=30)
# p.vals = 1-exp(-exp(x.axis))[-c(1:6)]  # needed, because of gaps in data due to zero-inflation
p.vals = seq(0.95,0.999,by=0.0025)
# u.vals = qexp(p.vals)
# p.vals = seq(0.8,0.99,by=0.01)

# window.sizes = c(seq(10,50,by=10),NA)  # crop radii (in km)
# finite.rads = window.sizes[is.finite(window.sizes)]

m="auto"

library(argparse)
parser <- ArgumentParser()
parser$add_argument("--loc", type = "integer", default=1, required = TRUE)
args <- parser$parse_args()
loc.idx <- args$loc

# for (loc.idx in 1:nrow(cities)) {
cat("Estimating perimeters for location",loc.idx,"\n")

# # distance (in metres) from every cell to this city, used instead of buffer()+mask()
# dist.r   <- terra::distance(r[[1]], pts_proj[loc.idx, ], unit = "m")
# dist.mat <- terra::as.matrix(dist.r, wide = TRUE)
# 
# # boolean "inside circle" mask per finite radius, reused for every obs.idx and every p
# circle.masks <- lapply(finite.rads, function(rad) dist.mat <= (rad * 1000))
# names(circle.masks) <- finite.rads

# target cell (row/col) for the "no crop" connected-component case
cell0 <- terra::cellFromXY(r, crds(pts_proj[loc.idx, ]))
rc0   <- terra::rowColFromCell(r, cell0)
ix <- rc0[2]
iy <- nrow(r) - rc0[1] + 1

# # precompute buffers once per location (independent of obs.idx and p)
# bufs <- lapply(window.sizes[!sapply(window.sizes, is.na)], function(rad) {
#   buffer(pts_proj[loc.idx, ], width = rad * 1000)
# })
# names(bufs) <- window.sizes[!sapply(window.sizes, is.na)]
# 
# # precompute the target cell (row/col) for the "no crop" connected-component case
# cell0 <- cellFromXY(r, crds(pts_proj[loc.idx, ]))
# rc0 <- rowColFromCell(r, cell0)
# ix <- rc0[2]; iy <- nrow(r) - rc0[1] + 1

# Crop 
RR = 100  # in km
r.obs = r[[5]]
buf <- buffer(pts_proj[loc.idx,], width = RR*1000)  # width input is in meters
r.obs.crop = mask(r.obs,buf)
plot(r.obs)
plot(r.obs.crop)

u.p = qexp(p.vals[1])
bf1 = (r.obs > u.p)*1
bf2 = (r.obs.crop > u.p)*1



results <- vector("list", 0) #vector("list", dim(r)[3])  # or a growable list, appended with [[]]

n.obs = dim(r)[3]
for (obs.idx in 1:n.obs) {
  cat("\tProgress:",round(100*obs.idx / n.obs,2),"%\n")
  
  if (vals.at.s0[obs.idx, loc.idx] < qexp(p.vals[1])) next
  
  # break
  r.i <- r[[obs.idx]]
  r.i.mat <- terra::as.matrix(r.i, wide = TRUE)  # this is a bit slow
  
  # # precompute cropped/masked rasters per rad, once per obs.idx (no p dependence)
  # # THIS IS SLOW:
  # r.i.by.rad <- lapply(window.sizes, function(rad) {
  #   if (is.na(rad)) r.i else mask(r.i, bufs[[as.character(rad)]])
  # })
  window.sizes = c(NA)
  for (p.val in p.vals) {
    u <- qexp(p.val)
    
    row <- data.frame(obs.idx=obs.idx, p=p.val, window.size=window.sizes,
                      perim_biased=NA_real_, perim_unbiased=NA_real_)
    
    if (vals.at.s0[obs.idx, loc.idx] < u) {
      row$perim_biased = rep(0,length(window.sizes))
      row$perim_unbiased = rep(0,length(window.sizes))
      
      results[[length(results) + 1]] <- row
      
      next
    }
    
    k=1
    # for (k in seq_along(window.sizes)) {
    #   rad <- window.sizes[[k]]

      # if (is.na(rad)) {
        bf <- (r.i.mat > u) * 1#terra::as.matrix(r.crop > u, wide = TRUE)  # this is a bit slow
        im <- as.cimg(t(bf[nrow(bf):1, ]))   
        lbl <- label(im, high_connectivity = FALSE)
        target_label <- lbl[ix, iy, 1, 1]
        bf.cc <- as.matrix(lbl == target_label) * 1
        # image(obs.idx)
        row$perim_biased[k]   <- perim.L1.na.rm(bf.cc)
        row$perim_unbiased[k] <- perimBF(bf.cc, p = 2, incl_boundary = FALSE, box_size = "auto")
      # } else {
      #   # bf <- terra::as.matrix(r.crop > u, wide = TRUE)
      #   # bf = t(bf[nrow(bf):1, ])
      #   # bf = bf*1
      #   # crop a circle surrounding the city via precomputed distance mask (no buffer()/mask() call)
      #   bf <- (r.i.mat > u) * 1
      #   bf[!circle.masks[[as.character(rad)]]] <- NA
      #   bf = t(bf[nrow(bf):1, ])  # re-orientate (not necessary)
      #   
      #   row$perim_biased[k]   <- perim.L1.na.rm(bf)
      #   bf.na0 = ifelse(is.na(bf),0,bf)  # this fixes the bug in the unbiased L2 perimeter
      #   row$perim_unbiased[k] <- perimBF(bf.na0, p = 2, incl_boundary = F, box_size = "auto")
      # }
    # }
    results[[length(results) + 1]] <- row
  }
  res.df <- do.call(rbind, results)
  save(res.df, file = paste0("perims_loc", loc.idx, "_better_v2.RData"))
}

stop("Finished.")
  
##################

library(dplyr)

beta.all = lapply(1:3, function(loc.idx) {
  load(paste0("/scratch/rcampbell/nice/Comephore_data/perims_loc",loc.idx,"_better_v2.RData"))
  # View(res.df)
  # head(res.df,10)
  
  # for each amount of cropping, find beta for all combination of values p1 and p2
  
  # # get a list of lists, each for a window size
  # res.split <- split(res.df, res.df$window.size, drop = FALSE)
  # res.split[[length(res.split)+1]] = res.df[is.na(res.df$window.size),]
  res.split = list(res.df)
  
  # 
  p.vals = sort(unique(res.df$p))
  p.idx.pairs = t(combn(1:length(p.vals),2))
  # p.idx.pairs = p.idx.pairs[apply(p.idx.pairs,1,function(vec) all(vec>5)),]
  # p.idx.pairs = p.idx.pairs[apply(p.idx.pairs,1,function(vec) all(vec>2 & vec < 7)),]
  beta.est.lst = lapply(res.split, function(lst){
    lst = lst %>% select(-window.size,-obs.idx)
    beta.med.res.biased = apply(p.idx.pairs,1,function(p.idx){
      p.idx1=p.idx[1]
      p.idx2=p.idx[2]
      perim1 = lst %>% filter(p==p.vals[p.idx1]) %>% select(perim_biased) %>% unlist() %>% unname() #%>% sample()
      perim2 = lst %>% filter(p==p.vals[p.idx2]) %>% select(perim_biased) %>% unlist() %>% unname() #%>% sample()
      beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
        log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
      # grp.len = 50
      # beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
      #                    function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)]))
      return(median(beta.est[is.finite(beta.est)]))
    })
    beta.med.res.unbiased = apply(p.idx.pairs,1,function(p.idx){
      p.idx1=p.idx[1]
      p.idx2=p.idx[2]
      perim1 = lst %>% filter(p==p.vals[p.idx1]) %>% select(perim_unbiased) %>% unlist() %>% unname() #%>% sample()
      perim2 = lst %>% filter(p==p.vals[p.idx2]) %>% select(perim_unbiased) %>% unlist() %>% unname() #%>% sample()
      beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
        log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
      # grp.len = 50
      # beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
      #                    function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)]))
      return(median(beta.est[is.finite(beta.est)]))
    })
    return(cbind(beta.bias=beta.med.res.biased,
                 beta.unbias=beta.med.res.unbiased))
  })
})

# beta.true.all.loc = c(1.41,1.19,0.4,0.9)
# for(rad.idx in 1:6){
rad.idx = 1
par(mfrow=c(1,3))
for(loc.idx in 1:3){
  beta.est.lst = beta.all[[loc.idx]]
  
  # pdf(paste0("beta_hist_loc",loc.idx,".pdf"),width=6,height=5)
  # par(mfrow=c(1,1),pty="m")
  hist(beta.est.lst[[rad.idx]][,1],breaks=20,
       main=NA,
       xlab=expression("Estimated"~beta),
       freq=F)
  # # abline(v=beta.true.all.loc[loc.idx],col="black",lwd=2)
  # abline(v=mean(beta.est.lst[[rad.idx]][,1],na.rm=T),col="green",lwd=2)
  abline(v=median(beta.est.lst[[rad.idx]][,1],na.rm=T),col="blue",lwd=2)
  # dev.off()
  print(round(median(beta.est.lst[[rad.idx]][,1],na.rm=T),2))
  
}
# }

# best is clearly CC, no cropping

## try median across all combinations of p1 and p2
rad.idx = 1
for(loc.idx in 1:4){
  beta.est.lst = beta.all[[loc.idx]]
  print(round(median(beta.est.lst[[rad.idx]][,1],na.rm=T),2))
}

# # Alternatively, try lowest absolute error to "true" values (not optimal)
# rad.idx = 6
# errors = NULL
# for(loc.idx in 1:4){
#   beta.est.p1p2 = beta.all[[loc.idx]][[rad.idx]][,1]
#   errors = cbind(errors,abs(beta.est.p1p2 - beta.true.all.loc[loc.idx]))
#   # best.idx = sort(errors,index.return=T)$ix
#   # print(best.idx)
# }
# total.error = apply(errors,1,sum)
# best.p1p2 = which.min(total.error)
# rad.idx = 6
# for(loc.idx in 1:4){
#   beta.est.lst = beta.all[[loc.idx]]
#   print(round(beta.est.lst[[rad.idx]][best.p1p2,1],2))
# }



# OLD, SLOWER VERSION
# for(loc.idx in c(1:nrow(cities))){
#   res.df = data.frame(obs.idx=numeric(),
#                       p=numeric(),
#                       window.size=numeric(),
#                       perim_biased=numeric(),
#                       perim_unbiased=numeric()
#                       )
#   row.idx = 1
#   for(obs.idx in c(1:dim(r)[3])){
#     if(vals.at.s0[obs.idx,loc.idx] < qexp(p.vals[1])){
#       next
#     }
#     
#     r.i = r[[obs.idx]]
#     
#     for(p in p.vals){
#       u = qexp(p)
#       if(vals.at.s0[obs.idx,loc.idx] < u){
#         break
#       }
#       
#       perims.biased = c()
#       perims.unbiased = c()
#       for(rad in window.sizes){
#         if(is.null(rad)){
#           # no cropping - instead get the perimeter of the connected component covering the location (city) of interest
#           bf = as.vector((r.i>u)*1)
#           bf <- matrix(bf, nrow = sqrt(length(bf)), ncol = sqrt(length(bf)), byrow = TRUE)
#           bf = t(bf)
#           bf <- bf[, nrow(bf):1]
#           im <- as.cimg(bf)
#           lon=cities[loc.idx,2]
#           lat=cities[loc.idx,3]
#           p <- vect(
#             data.frame(x = lon, y = lat),
#             geom = c("x", "y"),
#             crs = "EPSG:4326"
#           )
#           p_proj <- project(p, crs(r))
#           cell <- cellFromXY(r, crds(p_proj))
#           rc <- rowColFromCell(r, cell)
#           row <- rc[1]
#           col <- rc[2]
#           ix <- col
#           iy <- nrow(r) - row + 1
#           lbl <- label(im, high_connectivity = FALSE)
#           target_label <- lbl[ix, iy, 1, 1]
#           keep <- (lbl == target_label)
#           bf.cc = as.matrix(keep)
#           
#           perim.biased = perim.L1.na.rm(bf.cc)
#           perim.unbiased = perimBF(bf.cc, p = 2, incl_boundary = FALSE, box_size = "auto")
#           perims.biased = c(perims.biased,perim.biased)
#           perims.unbiased = c(perims.unbiased,perim.unbiased)
#         } else {
#           # crop a circle surrounding the city of interest, get perimeter of excursion sets
#           
#           buf <- buffer(pts_proj[loc.idx,], width = rad*1000)  # width input is in meters
#           r.i.crop = mask(r.i,buf)
#           
#           bf = as.vector((r.i.crop>u)*1)
#           bf <- matrix(bf, nrow = sqrt(length(bf)), ncol = sqrt(length(bf)), byrow = TRUE)
#           bf = t(bf)
#           bf <- bf[, nrow(bf):1]
#           
#           perim.biased = perim.L1.na.rm(bf)
#           perim.unbiased = perimBF(bf, p = 2, incl_boundary = FALSE, box_size = "auto")
#           perims.biased = c(perims.biased,perim.biased)
#           perims.unbiased = c(perims.unbiased,perim.unbiased)
#         }
#       }
#       
#       res.df = rbind(res.df,
#                      cbind(rep(obs.idx,length(window.sizes)),
#                            rep(p,length(window.sizes)),
#                            window.sizes,
#                            perims.biased,perims.unbiased)
#                      )
#       save(res.df,file=paste0("perims_loc",loc.idx,".RData"))
#     }
#   }
# }
