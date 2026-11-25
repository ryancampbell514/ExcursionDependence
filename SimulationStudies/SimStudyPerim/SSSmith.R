rm(list=ls())
library(this.path)
library(SpatialExtremes)
library(fields)
library(geoR)
library(spsUtil)
library(dplyr)
setwd(this.path::here())
source("../Functions.R")

library(Rcpp)
sourceCpp("../excursion-sets/binary_fields.cpp")
source("../excursion-sets/excursions.R")

#############################################################################

m="auto" # was previously 7, but doesn't make much a difference to perimeter esimtates
global.length = 30 #floor(10*m^(3/2))
half.width <- global.length/2
num.pixels = 2^10
pixel.width = global.length / num.pixels
x = y = seq(-half.width, half.width, length.out = num.pixels)  # c(-0.5,0.5), grid spacing 0.02
S.coord.grid = expand.grid(x,y)

# x.obs = SpatialExtremes::rmaxstab(1, cbind(x, y), "gauss", 
#                                   cov11 = 1, cov12 = 0, cov22 = 1, 
#                                   grid = TRUE)
# plotRandomField.v2(rf=x.obs,x=x,y=y)
# plotRandomField.v2(rf=x.obs,x=x,y=y,thresh=50)

# get p and u vals (using empirical distribution)
# p.vals =  c(seq(0.7,0.95,by=0.05),0.99)
# p.vals = seq(0.99,0.999,by=0.001)
p.vals = seq(0.8,0.99,by=0.01)
# up = function(p) -1/log(p)  # Frechet margins
# u.vals = up(p.vals)

# C1 = function(p,ell=range){
#   sqrt(pi/2)*(1/ell) *p*(-log(p))
# }
# p.low = seq(0.6,0.99,by=0.01)
# p.up = 0.5*(1+p.low)#p.low+0.005
# betap1p2.true = apply(cbind(p.low,p.up),1,function(p1p2){
#   (log(C1(p1p2[1]))-log(C1(p1p2[2])) -
#      log((1-p1p2[1])/(1-p1p2[2]))) /
#     log(log(1-p1p2[1])/log(1-p1p2[2]))
# })
# plot(p.low,betap1p2.true,type="l",xlab="p1 (p2=(p1+1)/2)",ylim=c(-0.05,0.7),
#      ylab=expression(beta(p1,p2)),lwd=2)
# abline(h=0,lty=3,col="red",lwd=2)

###################################################################################

# Run the study
N.rep = 200*50  # number of repetitions x number of betas
res.df = data.frame(matrix(NA,ncol = 3, nrow = 1))  # keep track of estimated perimeters
colnames(res.df) = c("mc.rep","p","perim")
row.idx=1
for(rep in 1:N.rep){
  cat("Monte Carlo simulation:",rep,"\n")
  
  rf = SpatialExtremes::rmaxstab(1, cbind(x, y), "gauss",
                                    cov11 = 1, cov12 = 0, cov22 = 1,
                                    grid = TRUE)
  # perims = NULL
  for(p in p.vals){
    u = -1/log(p)
    bf = (rf > u)*1
    perim.val = perimBF(bf,p = 2, incl_boundary = FALSE, box_size = m)
    perim.val = pixel.width * perim.val
    res.df[row.idx,] = c(rep,p,perim.val)
    row.idx = row.idx + 1
  }
  
  rm(bf,rf)
  save.image("SSSmith_lowerthresh.RData")
}
stop("End of simulation study.")

##########################

rm(list=ls())
library(dplyr)
setwd(this.path::here())
load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSSmith_lowerthresh.RData")
res.df = res.df[1:120260,]
head(res.df)  #120260

C1.smith = function(p,sd=1){
  (1/(sd*sqrt(2*pi)))*(pi/2)*p*(-log(p))
}

# # beta(p1,p2) plots
# p.low = 1-10^seq(log10(0.4), -14, length.out = 1000)
# p.up = 0.5*(p.low+1)
# betap1p2.smith.true = apply(cbind(p.low,p.up),1,function(p1p2){
#   (log(C1.smith(p1p2[1]))-log(C1.smith(p1p2[2])) -
#      log((1-p1p2[1])/(1-p1p2[2]))) /
#     log(log(1-p1p2[1])/log(1-p1p2[2]))
# })
# pdf("../figs/beta_p1p2_smith.pdf",width=6.5,height=6)
# par(mfrow=c(1,1),pty="m",mar = c(5.1, 5.5, 4.1, 2.1),cex.lab=1.5)
# plot(p.low,betap1p2.smith.true,type="l",xlab=expression(p),ylim=c(-0.05,0.7),
#      ylab=expression(beta(p,0.5(p+1))),lwd=2)
# abline(h=0.0,lty=3,col="red",lwd=2)
# dev.off()
# 
# Perimeter density estimation plots
res.df.perims = res.df %>% mutate(C1.est = 0.5*perim / global.length^2)
pp = seq(p.vals[1],1,length.out=1000)
C1.true = sapply(pp,C1.smith)
pdf("../figs/C1Smith.pdf",width=5,height=4.5)
plot(NULL,
     xlim=range(pp),ylim=range(c(C1.true,res.df.perims$C1.est),na.rm=T),
     xlab="p",ylab=expression(C*"*"[1](p)))
spsUtil::quiet(
  sapply(p.vals[c(1:length(p.vals))[c(1:length(p.vals))%%2 ==0]],function(p.val){
    res.df.perims %>% filter(p==p.val) %>% select(C1.est) %>%
      boxplot(add=T,at=p.val,boxwex = 0.03, axes = F, outline=F)
  })
)
lines(pp,C1.true,lwd=2,col="red")
dev.off()

# p.vals=c(0.990,0.991)#c(0.992,0.993)
p.vals = p.vals[-c(1:10)]
# p.vals = p.vals[-c(1:5)]
p.idx.pairs = t(combn(1:length(p.vals),2))
# p.idx.pairs = p.idx.pairs[1:9,]
# set.seed(12)
grp.len=30
beta.med.res = apply(p.idx.pairs,1,function(p.idx){
  p.idx1=p.idx[1]
  p.idx2=p.idx[2]
  # set.seed(13)
  perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
    log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
  # beta.med = beta.est
  beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
                     function(start.idx) {
                       median((beta.est[start.idx:(start.idx+grp.len)])[is.finite(beta.est[start.idx:(start.idx+grp.len)])])
                    })
  
  beta.true = (log(C1.smith(p.vals[p.idx1])/C1.smith(p.vals[p.idx2])) - 
                 log((1-p.vals[p.idx1])/(1-p.vals[p.idx2]))) / 
    log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
  return(list(beta.med=beta.med,beta.true=beta.true))
})

# beta.med.res.smith = beta.med.res
# save(beta.med.res.gauss, file="res_gauss.RData")

bias.skeweness <- data.frame(
  setting = seq_along(beta.med.res),
  mean = sapply(beta.med.res, function(lst) mean(lst$beta.med[is.finite(lst$beta.med)])),
  bias = sapply(beta.med.res, function(lst) {
    mean(lst$beta.med[is.finite(lst$beta.med)]) - lst$beta.true
  }),
  abs_bias = abs(sapply(beta.med.res, function(lst) {
    mean(lst$beta.med[is.finite(lst$beta.med)]) - lst$beta.true
  })),
  skewness = sapply(beta.med.res, function(lst) {
    moments::skewness(lst$beta.med[is.finite(lst$beta.med)])
  })
)

bias.skeweness$rank_bias <- rank(bias.skeweness$abs_bias)
bias.skeweness$rank_skew <- rank(abs(bias.skeweness$skewness))
bias.skeweness$total_rank <- bias.skeweness$rank_bias + bias.skeweness$rank_skew
bias.skeweness <- bias.skeweness[order(bias.skeweness$total_rank), ]

best.p = 25
# best.p = bias.skeweness[1,1]
# best.p = which.min(sapply(beta.med.res,function(lst) abs(median(lst$beta.med[is.finite(lst$beta.med)])-lst$beta.true)))
p1p2.smith = p.vals[p.idx.pairs[best.p,]]
beta.smith.est =  beta.med.res[[best.p]]$beta.med
beta.smith.est = beta.smith.est[is.finite(beta.smith.est)]
beta.smith.true = beta.med.res[[best.p]]$beta.true
pdf("../figs/beta_boxplot_Smith.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s")
boxplot(beta.smith.est,ylim=c(-1,1),ylab = "scaling index")
abline(h=beta.smith.true,lwd=2,col="red")
dev.off()

# test = (sqrt(global.length^2))*(((beta.smith.est)))/sd(beta.smith.est)  > qnorm(0.95)
test = (((beta.smith.est)))/sd(beta.smith.est) > qnorm(0.95)
Pr.reject = mean(test==1,na.rm=T)  # fail to reject at the 5% level
round(Pr.reject,3)

Z.smith = (beta.smith.est-beta.smith.true)/sd(beta.smith.est)
qq.smith = sort(Z.smith)
qq.true = qnorm(c(1:length(qq.smith))/(length(qq.smith)+1))

# bootstrap CIs
bs.reps=200
idx.by.rep <- split(seq_len(nrow(res.df)), res.df$mc.rep)
rep.names <- names(idx.by.rep)
n.rep <- length(rep.names)
set.seed(123)
beta.smith.bs = lapply(1:bs.reps, function(bs.rep){
  sampled.reps <- sample(rep.names, n.rep, replace = TRUE)
  idx.block <- unlist(
    idx.by.rep[sampled.reps],
    use.names = FALSE
  )
  res.df.bs = res.df[idx.block,]
  perim1.bs = res.df.bs %>% filter(p==p1p2.smith[1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  perim2.bs = res.df.bs %>% filter(p==p1p2.smith[2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  beta.est.bs = (log(perim1.bs) - log(perim2.bs) - log((1-p1p2.smith[1])/(1-p1p2.smith[2])))/
    log(log(1-p1p2.smith[1])/log(1-p1p2.smith[2]))
  # beta.med = beta.est
  beta.med.bs =  sapply((1:floor(length(beta.est.bs)/grp.len))*grp.len - grp.len+1,
                        function(start.idx) median(beta.est.bs[start.idx:(start.idx+grp.len)]))
  return(beta.med.bs)
})

Z.smith.bs = lapply(beta.smith.bs, function(vec) (vec-beta.smith.true)/sd(vec,na.rm=T))
qq.smith.bs = lapply(Z.smith.bs, function(vec) sort(vec))
qq.smith.bs.mat = do.call(rbind,qq.smith.bs)
qq.smith.low = apply(qq.smith.bs.mat,2,quantile,p=0.025)
qq.smith.upp = apply(qq.smith.bs.mat,2,quantile,p=0.975)
qq.true.bs = qnorm(c(1:length(qq.smith.low))/(length(qq.smith.low)+1))

CI.x = c(qq.true.bs,rev(qq.true.bs))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(qq.smith.low,rev(qq.smith.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0

pdf("../figs/QQbeta_smith.pdf",5,5)
par(mfrow=c(1,1),pty="s")
plot(NULL,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.smith,qq.true)),
     ylim=range(c(qq.smith,qq.true)))
polygon(x=CI.x,y=CI.y,col=scales::alpha("grey",0.6),border=NA)
points(qq.true,qq.smith,pch=20)
abline(0,1)
dev.off()

