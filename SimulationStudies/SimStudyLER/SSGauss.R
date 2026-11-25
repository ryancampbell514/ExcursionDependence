rm(list=ls())
library(this.path)
library(SpatialExtremes)
library(fields)
library(geoR)
library(this.path)
library(plotly)
setwd(this.path::here())

library(ExcursionSets)

library(quantreg)

source("../Functions.R")

##############################################################################

HALF_WIDTH <- 15#0.5
SIDE_LENGTH <- 121
x = y = seq(-HALF_WIDTH, HALF_WIDTH, length.out = SIDE_LENGTH)  # c(-0.5,0.5), grid spacing 0.02
S.coord.grid = expand.grid(x,y)
cond_site = data.frame("x" = 0,"y" = 0)

# Generate GRF
range=9#0.3
sigma=1
cov.fn = function(h) {
  h = h/range
  return(sigma^2 * exp(-(h^2)))
}
gen = random.field.generator(grid=S.coord.grid,
                             cov_fun=cov.fn,
                             cond_sites = cond_site)
x.axis = seq(1.097189,2.443470,length.out=50) # to have matching x axis on all plots on the log(-log(1-p))
p.vals = 1-exp(-exp(x.axis))
u.vals = qnorm(p.vals)  # the threshold values, marginal quantile function applied

N.beta = 200
N.rep = 30
beta.est.vec = c()

for(rep in 1:N.beta){
  print(paste("Simulation Study repetition",rep))
  
  QR.data.mat = NULL
  
  for(u in u.vals){
    epys = c()
    lers = c()
    
    for (simulation in 1:N.rep){
      rf = generateGaussian.truncnorm(generator=gen, cond_val=u, exceedance = TRUE)
      ep = extent.profile(rf, u, x0 = cond_site)
      epys = cbind(epys, ep$y)
      ler = proportion.range(ep, 1)
      ler = ifelse(is.infinite(ler),1/sqrt(2),ler)
      QR.data.mat = rbind(QR.data.mat,c(log(ler),-log(-log(1-pnorm(u)))))
    }
  }
  
  QR.data.df = data.frame(X=QR.data.mat[,2],Y=QR.data.mat[,1])
  fit = rq(Y~X, data=QR.data.df)
  beta.est = coefficients(fit)[2]
  
  beta.est.vec = c(beta.est.vec,beta.est)
  save(beta.est.vec,file="SSGauss.RData")
  print(beta.est.vec)
}

