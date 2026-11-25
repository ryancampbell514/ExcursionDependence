
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

k.df=3
# library(argparse)
# parser <- ArgumentParser()
# parser$add_argument("--df", type = "integer", default=3, required = TRUE)
# args <- parser$parse_args()
# k.df <- args$df

range=1#(num.pixels/100)/global.length#3  # smaller range, rougher field, higher \lambda_2
rho <- make_rho_squexp(l=range, sigma2=1)
lambda.true = 1/(range^2)
lambda.vals = NULL

# get p and u vals (using empirical distribution)
# p.vals = seq(0.99,0.999,by=0.001)
p.vals = seq(0.8,0.99,by=0.01)

###################################################################################

# Run the study
N.rep = 50*200  # number of repetitions x number of betas
res.df = data.frame(matrix(NA,ncol = 3, nrow = 1))  # keep track of estimated perimeters
colnames(res.df) = c("mc.rep","p","perim")
row.idx=1
for(rep in 1:N.rep){
  cat("Monte Carlo simulation:",rep,"\n")
  
  rf = lapply(1:k.df,function(ii){
    obs = stationary_Gaussian_process_EDiBcode(nx=num.pixels,ny=num.pixels,
                                               xlim=c(-global.length/2,global.length/2),
                                               ylim=c(-global.length/2,global.length/2),
                                               rho=rho)
    return(obs$field)
  })
  rf = simplify2array(rf)
  rf = apply(rf,c(1,2),function(vec)( sum(vec^2) - k.df) / sqrt(2*k.df))
  
  # perims = NULL
  for(p in p.vals){
    u = (qchisq(p,df=k.df)-k.df)/sqrt(2*k.df)  
    bf = (rf > u)*1
    perim.val = perimBF(bf,p = 2, incl_boundary = FALSE, box_size = m)
    perim.val = pixel.width * perim.val
    res.df[row.idx,] = c(rep,p,perim.val)
    row.idx = row.idx + 1
  }
  
  rm(bf,rf)
  base::save.image(paste0("SSChisq_df",k.df,"_lowerthresh_higherrange.RData"))
}

##########################

rm(list=ls())
library(dplyr)
setwd(this.path::here())
load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSChisq_df3_lowerthresh_higherrange.RData")
res.df = res.df[1:120260,]
head(res.df)

C1.chisq = function(p,k=k.df,lam=lambda.true){
  u = (qchisq(p,df=k)-k)/sqrt(2*k)
  (sqrt(pi*lam)/(2^(0.5*(k.df+1))*gamma(k.df/2)))*
    (k+ sqrt(2*k)*u)^(0.5*(k-1)) *
    exp(-0.5*(k+ sqrt(2*k)*u))
}

# # beta(p1,p2) plots
# p.low = 1-10^seq(log10(0.4), -14, length.out = 1000)
# p.up = 0.5*(p.low+1)
# betap1p2.chisq.true = apply(cbind(p.low,p.up),1,function(p1p2){
#   (log(C1.chisq(p1p2[1]))-log(C1.chisq(p1p2[2])) -
#      log((1-p1p2[1])/(1-p1p2[2]))) /
#     log(log(1-p1p2[1])/log(1-p1p2[2]))
# })
# pdf("../figs/beta_p1p2_chisq.pdf",width=6.5,height=6)
# par(mfrow=c(1,1),pty="m",mar = c(5.1, 5.5, 4.1, 2.1),cex.lab=1.5)
# plot(p.low,betap1p2.chisq.true,type="l",xlab=expression(p),ylim=c(-0.05,0.7),
#      ylab=expression(beta(p,0.5(p+1))),lwd=2)
# abline(h=0.5,lty=3,col="red",lwd=2)
# dev.off()
# 
# Perimeter density estimation plots
res.df.perims = res.df %>% mutate(C1.est = 0.5*perim / global.length^2)
pp = seq(p.vals[1],1,length.out=1000)
C1.true = sapply(pp,C1.chisq)
pdf("../figs/C1Chisq.pdf",width=5,height=4.5)
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
  
  beta.true = (log(C1.chisq(p.vals[p.idx1])/C1.chisq(p.vals[p.idx2])) - 
                 log((1-p.vals[p.idx1])/(1-p.vals[p.idx2]))) / 
    log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
  return(list(beta.med=beta.med,beta.true=beta.true))
})
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

best.p = 170  # for the Smith model
# best.p = bias.skeweness[1,1]#which.min(sapply(beta.med.res,function(lst) abs(median(lst$beta.med[is.finite(lst$beta.med)])-lst$beta.true)))
p1p2.chisq = p.vals[p.idx.pairs[best.p,]]
beta.chisq.est =  beta.med.res[[best.p]]$beta.med
beta.chisq.est = beta.chisq.est[is.finite(beta.chisq.est)]
beta.chisq.true = beta.med.res[[best.p]]$beta.true
pdf("../figs/beta_boxplot_Chisq.pdf",width=5,height=5)
par(mfrow=c(1,1),pty="s")
boxplot(beta.chisq.est,ylim=c(0,1),ylab = "scaling index")
abline(h=beta.chisq.true,lwd=2,col="red")
dev.off()

# test = (sqrt(global.length^2))*(((beta.chisq.est)))/sd(beta.chisq.est)  > qnorm(0.95)
test = (((beta.chisq.est)))/sd(beta.chisq.est)  > qnorm(0.95)
Pr.reject = mean(test==1,na.rm=T)  # fail to reject at the 5% level
round(Pr.reject,3)

Z.chisq = (beta.chisq.est-beta.chisq.true)/sd(beta.chisq.est)
qq.chisq = sort(Z.chisq)
qq.true = qnorm(c(1:length(qq.chisq))/(length(qq.chisq)+1))

# bootstrap CIs
bs.reps=200
idx.by.rep <- split(seq_len(nrow(res.df)), res.df$mc.rep)
rep.names <- names(idx.by.rep)
n.rep <- length(rep.names)
set.seed(123)
beta.chisq.bs = lapply(1:bs.reps, function(bs.rep){
  sampled.reps <- sample(rep.names, n.rep, replace = TRUE)
  idx.block <- unlist(
    idx.by.rep[sampled.reps],
    use.names = FALSE
  )
  res.df.bs = res.df[idx.block,]
  perim1.bs = res.df.bs %>% filter(p==p1p2.chisq[1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  perim2.bs = res.df.bs %>% filter(p==p1p2.chisq[2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  beta.est.bs = (log(perim1.bs) - log(perim2.bs) - log((1-p1p2.chisq[1])/(1-p1p2.chisq[2])))/
    log(log(1-p1p2.chisq[1])/log(1-p1p2.chisq[2]))
  # beta.med = beta.est
  beta.med.bs =  sapply((1:floor(length(beta.est.bs)/grp.len))*grp.len - grp.len+1,
                        function(start.idx) median(beta.est.bs[start.idx:(start.idx+grp.len)]))
  return(beta.med.bs)
})

Z.chisq.bs = lapply(beta.chisq.bs, function(vec) (vec-beta.chisq.true)/sd(vec,na.rm=T))
qq.chisq.bs = lapply(Z.chisq.bs, function(vec) sort(vec))
qq.chisq.bs.mat = do.call(rbind,qq.chisq.bs)
qq.chisq.low = apply(qq.chisq.bs.mat,2,quantile,p=0.025)
qq.chisq.upp = apply(qq.chisq.bs.mat,2,quantile,p=0.975)
qq.true.bs = qnorm(c(1:length(qq.chisq.low))/(length(qq.chisq.low)+1))

CI.x = c(qq.true.bs,rev(qq.true.bs))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(qq.chisq.low,rev(qq.chisq.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0

pdf("../figs/QQbeta_chisq.pdf",5,5)
par(mfrow=c(1,1),pty="s")
plot(NULL,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.chisq,qq.true)),
     ylim=range(c(qq.chisq,qq.true)))
polygon(x=CI.x,y=CI.y,col=scales::alpha("grey",0.6),border=NA)
points(qq.true,qq.chisq,pch=20)
abline(0,1)
dev.off()
