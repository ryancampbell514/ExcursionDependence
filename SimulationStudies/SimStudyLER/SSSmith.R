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

# In this script, we use the custom R function to sample the conditional Smith model,
# described in Section 3.2 of the manuscript
rSmith_exceeding_storms <- function(u, sigma = 1) {
  
  stopifnot(u > 0, sigma > 0)
  
  lambda <- 1 / u
  
  # Zero-truncated Poisson number of storms ( to assure no 0 draws... need at least one storm )
  N <- qpois(
    runif(1, min = exp(-lambda), max = 1),
    lambda = lambda
  )
  
  # Storm centres
  U <- matrix(
    rnorm(2 * N, mean = 0, sd = sigma),
    ncol = 2
  )
  
  # Gaussian storm density evaluated at the origin
  f0 <- exp(-rowSums(U^2) / (2 * sigma^2)) /
    (2 * pi * sigma^2)
  
  # Storm magnitudes conditional on exceeding u at the origin
  # xi <- u / (f0 * rexp(N))
  xi <- u / (f0 * runif(N))
  
  data.frame(
    xi = xi,
    U1 = U[, 1],
    U2 = U[, 2]
  )
}

##############################################################################

HALF_WIDTH <- 15
SIDE_LENGTH <- 121
x = y = seq(-HALF_WIDTH, HALF_WIDTH, length.out = SIDE_LENGTH)  # c(-0.5,0.5), grid spacing 0.02
S.coord.grid = expand.grid(x,y)
cond_site = data.frame("x" = 0,"y" = 0)
cond_site_idx = which(apply(S.coord.grid,1,function(coord) all(coord==0)))

x.axis = seq(1.097189,2.443470,length.out=50) # to have matching x axis on all plots on the log(-log(1-p))
p.vals = 1-exp(-exp(x.axis))
u.vals = -1/log(p.vals)

N.beta = 200
N.rep = 50
beta.est.vec = c()

sigma.smith = 3

for(rep in 1:N.beta){
  print(paste("Simulation Study repetition",rep, "of", N.beta))
  
  QR.data.mat = NULL
  
  for(u in u.vals){
    epys = c()
    lers = c()
    
    for (simulation in 1:N.rep){
      
      storms <- rSmith_exceeding_storms(
        u = u,
        sigma = sigma.smith
      )
      Z <- matrix(0, nrow = length(y), ncol = length(x))
      
      for (i in seq_len(nrow(storms))) {
        
        d2 <- outer(
          y - storms$U2[i],
          x - storms$U1[i],
          function(a, b) a^2 + b^2
        )
        
        contribution <- storms$xi[i] *
          exp(-d2 / (2 * sigma.smith^2)) / (2 * pi * sigma.smith^2)
        
        Z <- pmax(Z, contribution)
      }
      
      if(sum(Z>u)==0){next}  # occasionally my sampling method fails
      
      rf = random.field(Z,S.coord.grid)
      
      rf.mat = matrix(rf$values,length(x),length(y))
      bf.mat = (rf.mat>=u) *1
      idx.0 <- which(bf.mat == 0, arr.ind = TRUE)
      nr <- nrow(bf.mat)
      nc <- ncol(bf.mat)
      origin <- c((nr + 1) / 2, (nc + 1) / 2)
      d <- sqrt((idx.0[, 1] - origin[1])^2 +
                  (idx.0[, 2] - origin[2])^2)
      min.dist = min(d[d>0])
      ler = min.dist * (2*HALF_WIDTH) / SIDE_LENGTH
      
      QR.data.mat = rbind(QR.data.mat,c(log(ler),-log(-log(1-exp(-1/u)))))
    }
  }
  
  QR.data.df = data.frame(X=QR.data.mat[,2],Y=QR.data.mat[,1])
  QR.data.df = QR.data.df %>% filter(Y>min(Y))

  fit = rq(Y~X, data=QR.data.df)
  beta.est = coefficients(fit)[2]
  
  beta.est.vec = c(beta.est.vec,beta.est)
  print(beta.est.vec)
  save(beta.est.vec,file="SSSmith.RData")
  
}

