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
pixel.width = 30/1023  # reduce this to (30/1023) or (30/1023)*2
num.pixels = ceiling(global.length/pixel.width)
if(num.pixels %% 2 !=0){num.pixels=num.pixels+1}
x = y = seq(-half.width, half.width, length.out = num.pixels)  # c(-0.5,0.5), grid spacing 0.02
S.coord.grid = expand.grid(x,y)

range=(num.pixels/100)/global.length#3  # smaller range, rougher field, higher \lambda_2
rho <- make_rho_squexp(l=range, sigma2=1)
lambda.true = 1/(range^2)
lambda.vals = NULL

# get p and u vals (using empirical distribution)
p.vals1 = seq(0.7,0.95,by=0.05)
p.vals2 = seq(0.99,0.999,by=0.001)
p.vals = c(p.vals1,p.vals2)

###################################################################################

# Run the study
N.rep = 200*50  # number of repetitions x number of betas
res.df = data.frame(matrix(NA,ncol = 3, nrow = 1))  # keep track of estimated perimeters
colnames(res.df) = c("mc.rep","p","perim")
row.idx=1
for(rep in 1:N.rep){
  cat("Monte Carlo simulation:",rep,"\n")
  
  rf = stationary_Gaussian_process_EDiBcode(nx=num.pixels,ny=num.pixels,
                                            xlim=c(-global.length/2,global.length/2),
                                            ylim=c(-global.length/2,global.length/2),
                                            rho=rho)
  rf = rf$field
  
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
  save.image("SSGauss_betterthresh.RData")
}
stop("End of simulation study.")

###########################################

rm(list=ls())
library(dplyr)
setwd(this.path::here())
# load("SSGauss.RData")
load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSGauss.RData")
head(res.df)
dim(res.df)
p.idx1=3  #3
p.idx2=4  #4
# p.len = floor(dim(res.df)[1]/2)
# perim1 = res.df[1:p.len,] %>% filter(p==p.vals[1]) %>% select(perim) %>% unlist() %>% unname()
# perim2 = res.df[(p.len+1):dim(res.df)[1],] %>% filter(p==p.vals[2]) %>% select(perim) %>% unlist() %>% unname()
perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
  log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
grp.len = 50
beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
                   function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)],na.rm=T))
beta.true = 0.5
pdf("SSGauss_beta.pdf",width=5,height=4.5)
par(mfrow=c(1,1),pty="m")
boxplot(beta.med)
abline(h=beta.true,col="red")
dev.off()

Z = (beta.med-beta.true)/sd(beta.med)
pdf("SSGauss_PPQQ.pdf",width=8,height=5)
par(mfrow=c(1,2),pty="s")
qq.model = sort(Z)
qq.theoretical = qnorm(c(1:length(qq.model))/(length(qq.model)+1))
plot(qq.theoretical,qq.model,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.theoretical,qq.model)),
     ylim=range(c(qq.theoretical,qq.model)),pch=20)
abline(0,1)

pp.model = pnorm(qq.model)
pp.theoretical = c(1:length(qq.model))/(length(qq.model)+1)
plot(pp.theoretical,pp.model,
     xlab="Theoretical",ylab="Model",
     xlim=c(0,1),
     ylim=c(0,1),pch=20)
abline(0,1)
dev.off()