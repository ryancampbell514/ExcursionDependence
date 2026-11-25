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

m="auto"
global.length = 30 
half.width <- global.length/2
num.pixels = 2^10
pixel.width = global.length / num.pixels
x = y = seq(-half.width, half.width, length.out = num.pixels) 
S.coord.grid = expand.grid(x,y)


# make the anisotropy
aniso = c(2,1.5)
aniso.mat = matrix(c(aniso[1],1,1,aniso[2]),2,2)
S.coord.grid.anis = t(apply(S.coord.grid,1,function(vec){
  vec = matrix(vec,ncol=1)
  return(as.numeric(aniso.mat %*% vec))
}))

range=1
rho <- make_rho_squexp(l=range, sigma2=1)
lambda.true = 1/(range^2)
lambda.vals = NULL

p.vals = seq(0.8,0.99,by=0.01)

###################################################################################

N.rep = 50*200  # number of repetitions x number of betas
res.df = data.frame(matrix(NA,ncol = 3, nrow = 1))  # keep track of estimated perimeters
colnames(res.df) = c("mc.rep","p","perim")
row.idx=1
for(rep in 1:N.rep){
  cat("Monte Carlo simulation:",rep,"\n")
  
  rf = stationary_Gaussian_process_EDiBcode(nx=num.pixels,ny=num.pixels,
                                            xlim=c(-global.length/2,global.length/2),
                                            ylim=c(-global.length/2,global.length/2),
                                            rho=rho,anis=aniso.mat)
  rf = rf$field
  
  plotRandomField.v2(rf,x,y)
  
  # perims = NULL
  for(p in p.vals){
    u = qnorm(p)
    bf = (rf > u)*1
    perim.val = perimBF(bf,p = 2, incl_boundary = FALSE, box_size = m)
    perim.val = pixel.width * perim.val
    res.df[row.idx,] = c(rep,p,perim.val)
    row.idx = row.idx + 1
  }
  
  rm(bf,rf)
  base::save.image("SSGaussAnisotropic.RData")
}
stop("End of simulation study.")

##########################

rm(list=ls())
library(dplyr)
setwd(this.path::here())
load("SSGaussAnisotropic.RData")

eigen.vec = eigen(aniso.mat)$values
ellipse_perimeter <- function(a, b) {
  h <- ((a - b)^2) / ((a + b)^2)
  pi * (a + b) * (1 + (3 * h) / (10 + sqrt(4 - 3 * h)))
}
anis.perim.const = ellipse_perimeter(a = eigen.vec[1], b = eigen.vec[2]) / (2*pi)
C1.gaussaniso = function(p,lam=lambda.true,a.const=anis.perim.const){
  0.25*sqrt(2*pi*lam)*a.const *dnorm(qnorm(p))
}

# Perimeter density estimation plots
res.df.perims = res.df %>% mutate(C1.est = 0.5*perim / global.length^2)
pp = seq(p.vals[1],1,length.out=1000)
C1.true = sapply(pp,C1.gaussaniso)
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

p.idx.pairs = t(combn(1:length(p.vals),2))
grp.len=30
beta.med.res = apply(p.idx.pairs,1,function(p.idx){
  p.idx1=p.idx[1]
  p.idx2=p.idx[2]
  perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
    log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
  beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
                     function(start.idx) {
                       median((beta.est[start.idx:(start.idx+grp.len)])[is.finite(beta.est[start.idx:(start.idx+grp.len)])])
                     })
  
  beta.true = (log(C1.gaussaniso(p.vals[p.idx1])/C1.gaussaniso(p.vals[p.idx2])) - 
                 log((1-p.vals[p.idx1])/(1-p.vals[p.idx2]))) / 
    log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
  return(list(beta.med=beta.med,beta.true=beta.true))
})

best.p = 170  # (p1,p2)=(0.93,0.94)
p1p2.gaussaniso = p.vals[p.idx.pairs[best.p,]]
beta.gaussaniso.est =  beta.med.res[[best.p]]$beta.med
beta.gaussaniso.est = beta.gaussaniso.est[is.finite(beta.gaussaniso.est)]
beta.gaussaniso.true = beta.med.res[[best.p]]$beta.true
par(mfrow=c(1,1),pty="s")
boxplot(beta.gaussaniso.est,ylim=c(0,1),ylab = "scaling index")
abline(h=beta.gaussaniso.true,lwd=2,col="red")

# test = (sqrt(global.length^2))*(((beta.gaussaniso.est)))/sd(beta.gaussaniso.est)  > qnorm(0.95)
test = (((beta.gaussaniso.est)))/sd(beta.gaussaniso.est)  > qnorm(0.95)
Pr.reject = mean(test==1,na.rm=T)  # fail to reject at the 5% level
round(Pr.reject,3)

Z.gaussaniso = (beta.gaussaniso.est-beta.gaussaniso.true)/sd(beta.gaussaniso.est)
qq.gaussaniso = sort(Z.gaussaniso)
qq.true = qnorm(c(1:length(qq.gaussaniso))/(length(qq.gaussaniso)+1))

# bootstrap CIs
bs.reps=200
idx.by.rep <- split(seq_len(nrow(res.df)), res.df$mc.rep)
rep.names <- names(idx.by.rep)
n.rep <- length(rep.names)
set.seed(123)
beta.gaussaniso.bs = lapply(1:bs.reps, function(bs.rep){
  sampled.reps <- sample(rep.names, n.rep, replace = TRUE)
  idx.block <- unlist(
    idx.by.rep[sampled.reps],
    use.names = FALSE
  )
  res.df.bs = res.df[idx.block,]
  perim1.bs = res.df.bs %>% filter(p==p1p2.gaussaniso[1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  perim2.bs = res.df.bs %>% filter(p==p1p2.gaussaniso[2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  beta.est.bs = (log(perim1.bs) - log(perim2.bs) - log((1-p1p2.gaussaniso[1])/(1-p1p2.gaussaniso[2])))/
    log(log(1-p1p2.gaussaniso[1])/log(1-p1p2.gaussaniso[2]))
  # beta.med = beta.est
  beta.med.bs =  sapply((1:floor(length(beta.est.bs)/grp.len))*grp.len - grp.len+1,
                        function(start.idx) median(beta.est.bs[start.idx:(start.idx+grp.len)]))
  return(beta.med.bs)
})

Z.gaussaniso.bs = lapply(beta.gaussaniso.bs, function(vec) (vec-beta.gaussaniso.true)/sd(vec,na.rm=T))
qq.gaussaniso.bs = lapply(Z.gaussaniso.bs, function(vec) sort(vec))
qq.gaussaniso.bs.mat = do.call(rbind,qq.gaussaniso.bs)
qq.gaussaniso.low = apply(qq.gaussaniso.bs.mat,2,quantile,p=0.025)
qq.gaussaniso.upp = apply(qq.gaussaniso.bs.mat,2,quantile,p=0.975)
qq.true.bs = qnorm(c(1:length(qq.gaussaniso.low))/(length(qq.gaussaniso.low)+1))

CI.x = c(qq.true.bs,rev(qq.true.bs))
CI.x[length(CI.x)+1] = CI.x[1]
CI.y = c(qq.gaussaniso.low,rev(qq.gaussaniso.upp))
CI.y[length(CI.y)+1] = CI.y[1]
CI.y[is.na(CI.y)] = 0

par(mfrow=c(1,1),pty="s")
plot(NULL,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.gaussaniso,qq.true)),
     ylim=range(c(qq.gaussaniso,qq.true)))
polygon(x=CI.x,y=CI.y,col=scales::alpha("grey",0.6),border=NA)
points(qq.true,qq.gaussaniso,pch=20)
abline(0,1)



