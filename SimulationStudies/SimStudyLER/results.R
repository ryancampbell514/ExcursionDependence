library(this.path)
setwd(this.path::here())

load("/scratch/rcampbell/nice/SimStudyLER_directsim/SSGauss_SqExp.RData")
gauss.res = beta.est.vec
load("/scratch/rcampbell/nice/SimStudyLER_directsim/SSGaussAnisotropic_SqExp.RData")
gaussaniso.res = beta.est.vec
load("/scratch/rcampbell/nice/SimStudyLER_directsim/SSChisq_SqExp.RData")
chisq.res = beta.est.vec
load("/scratch/rcampbell/nice/SimStudyLER_directsim/SSParetoGauss_SqExp.RData")
paretoguss.res = beta.est.vec
# load("SSStudent.RData")
# student.res = beta.est.vec
# load("SSEllipticalStudent.RData")
# ellipticstudent.res = beta.est.vec
# load("SSLaplace.RData")
# laplace.res = beta.est.vec
# load("SSGaussPareto.RData")
# gausspareto.res = beta.est.vec
# # load("SSStable.R")
# stable.res = NA
# load("SSrPareto.RData")
# rpareto.res = beta.est.vec
# load("SSMSdisks.RData")
# msdisks.res = beta.est.vec
# load("SSShotNoise.RData")
# shotnoise.res = beta.est.vec

# rm(beta.est.vec)

pdf("resultsLERbetas.pdf",width=8,height=6)
boxplot(gauss.res,
        gaussaniso.res,
        chisq.res,
        paretoguss.res,
        # laplace.res,
        # ellipticstudent.res,
        # stable.res,
        # rpareto.res,
        ylab = "scaling index",
        names = NA)
x_labs = c("Gauss.",
           "Gauss. Aniso.",
           "Chi Sq.",
           "Pareto * Gauss."
)
axis(1, at=1:length(x_labs), labels = x_labs)    # add tick marks
abline(h=c(0,0.5),lty=2)
dev.off()
