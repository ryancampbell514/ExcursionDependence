rm(list=ls())
library(this.path)
library(SpatialExtremes)
library(fields)
library(evd)
# library(mev)
library(geoR)
library(this.path)
library(plotly)
setwd(this.path::here())
#setwd("/home/rcampbell/nice/SimStudy")

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

# Generate GRF
# k.df=3
library(argparse)
parser <- ArgumentParser()
parser$add_argument("--df", type = "integer", default=3, required = TRUE)
args <- parser$parse_args()
k.df <- args$df
# range=0.3
# nu=2.5
# cov.fn = cov.model("Matern", range=range, nu = nu)
range=0.3
# range = 0.341
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
u.vals = (qchisq(p.vals,df=k.df)-k.df)/sqrt(2*k.df)  # the threshold values

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
      rf = generateChisq(generator=gen, k=k.df, thresh=(sqrt(2*k.df)*u)+k.df)
      rf$values = (rf$values - k.df)/sqrt(2*k.df)
      ep = extent.profile(rf, u, x0 = cond_site)
      epys = cbind(epys, ep$y)
      ler = proportion.range(ep, 1)
      ler = ifelse(is.infinite(ler),1/sqrt(2),ler)
      QR.data.mat = rbind(QR.data.mat,c(log(ler),-log(-log(1-pchisq((sqrt(2*k.df)*u)+k.df,df=k.df)))))
    }
  }
  
  QR.data.df = data.frame(X=QR.data.mat[,2],Y=QR.data.mat[,1])
  fit = rq(Y~X, data=QR.data.df)
  beta.est = coefficients(fit)[2]
  
  beta.est.vec = c(beta.est.vec,beta.est)
  print(beta.est.vec)
  save(beta.est.vec,file=paste0("SSChisq_df",k.df,"_SqExp.RData"))
}

# hist(beta.est.vec)

