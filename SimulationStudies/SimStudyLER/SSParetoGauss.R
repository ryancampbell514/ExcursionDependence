rm(list=ls())
library(this.path)
library(SpatialExtremes)
library(fields)
library(geoR)
library(this.path)
library(plotly)
library(rmutil)
setwd(this.path::here())

# install.packages("~/nice/Local-Excursion-Set-Coefficients/ExcursionSets_0.0.4.tar.gz")
library(ExcursionSets)

library(quantreg)

source("../Functions.R")

##############################################################################

HALF_WIDTH <- 0.5 
SIDE_LENGTH <- 121
x = y = seq(-HALF_WIDTH, HALF_WIDTH, length.out = SIDE_LENGTH)  # c(-0.5,0.5), grid spacing 0.02
S.coord.grid = expand.grid(x,y)
cond_site = data.frame("x" = 0,"y" = 0)
s0.idx = which(apply(S.coord.grid,1,function(vec) all(vec==0)))

# Generate GRF
# range=0.3
# nu=2.5
# cov.fn = cov.model("Matern", range=range, nu = nu)
range=0.3
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
u.vals = (2*(1-p.vals))^(-1/2)  # the threshold values, marginal quantile function applied

cdf.approx = function(x){ 1-(1/(2*x^2)) }
quantile.approx = function(p){ 1/sqrt(2*(1-p)) }

N.beta = 200
N.rep = 50
beta.est.vec = c()

for(rep in 1:N.beta){
  print(paste("Simulation Study repetition",rep))
  
  QR.data.mat = NULL
  
  for(u in u.vals){
    epys = c()
    lers = c()
    
    for (simulation in 1:N.rep){
      # Lam.obs = sqrt(rexp(1,0.5))
      # G.obs = generateGaussian.truncnorm(generator=gen, 
      #                                    cond_val=u/Lam.obs, 
      #                                    exceedance = TRUE)$values
      # Z = Lam.obs * G.obs
      
      x0 = quantile.approx(1-runif(1)*(1-cdf.approx(u)))
      nloc = dim(S.coord.grid)[1]
      log.target = function(x,thresh=x0) {
        log(2)-3*log(x)+ dnorm(thresh/x,log=T) - nloc*log(x)
      }
      lambda = logMetropolisAlg(lambda0=quantile.approx(0.5),
                              log.target.dens=log.target,
                              prop.fn = function(loc) exp(rnorm(1,mean=loc)))
      w0 = x0/lambda
      G.obs = generateGaussian(generator=gen, cond_val=w0, exceedance = FALSE)$values
      Z = G.obs*lambda
      
      if(Z[s0.idx] <= u){
        stop("Conditional sample is not larger than the threshold at the reference site.")
      }
      
      rf = random.field(Z, S.coord.grid)
      # plotRandomField.v2(rf,x,y,thresh=u)
      ep = extent.profile(rf, u, x0 = cond_site)
      epys = cbind(epys, ep$y)
      ler = proportion.range(ep, 1)
      ler = ifelse(is.infinite(ler),1/sqrt(2),ler)
      QR.data.mat = rbind(QR.data.mat,c(log(ler),-log(-log(1-cdf.approx(u)))))
    }
  }
  
  QR.data.df = data.frame(X=QR.data.mat[,2],Y=QR.data.mat[,1])
  fit = rq(Y~X, data=QR.data.df)
  beta.est = coefficients(fit)[2]
  # print(beta.est)
  beta.est.vec = c(beta.est.vec,beta.est)
  save(beta.est.vec,file="SSParetoGauss_SqExp.RData")
  print(beta.est.vec)
}

# hist(beta.est.vec)
# boxplot(beta.est.vec)
stop()

load("/scratch/rcampbell/nice/SimStudyLER_directsim/SSParetoGauss.RData")
print(beta.est.vec)
