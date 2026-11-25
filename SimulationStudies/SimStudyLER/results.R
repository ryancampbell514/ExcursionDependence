library(this.path)
setwd(this.path::here())

load("SSGauss.RData")
gauss.res = beta.est.vec
load("SSGaussAnisotropic.RData")
gaussaniso.res = beta.est.vec
load("SSChisq.RData")
chisq.res = beta.est.vec
load("SSSmith.RData")
smith.res = beta.est.vec

rm(beta.est.vec)

boxplot(gauss.res,
        gaussaniso.res,
        chisq.res,
        smith.res,
        ylab = "scaling index",
        names = NA)
x_labs = c("Gauss.",
           "Gauss. Aniso.",
           "Chi Sq.",
           "Smith"
)
axis(1, at=1:length(x_labs), labels = x_labs)    # add tick marks
abline(h=c(0,0.5),lty=2,lwd=1,col="red")
