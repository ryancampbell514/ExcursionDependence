rm(list=ls())
setwd(this.path::here())

grp.len = 50

load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSGauss.RData")
p.idx1=3  #3
p.idx2=4  #4
perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
  log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
# grp.len = 50
beta.gauss =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
                   function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)],na.rm=T))

load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSGaussAniso.RData")
p.idx1=4  #3
p.idx2=5  #4
perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
  log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
# grp.len = 50
beta.gaussaniso =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
                     function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)],na.rm=T))
# load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSGaussAniso_betterthresh.RData")
# p.vals = p.vals[1:7]
# # p.vals = p.vals[7:16]
# p.idx.pairs = t(combn(1:length(p.vals),2))
# beta.med.res = apply(p.idx.pairs,1,function(p.idx){
#   p.idx1=p.idx[1]
#   p.idx2=p.idx[2]
#   perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
#   perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
#   beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
#     log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
#   beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
#                      function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)]))
#   return(beta.med)
# })
# best.p = which.min(apply(beta.med.res,2,function(vec) abs(median(vec[is.finite(vec)])-0.5)))
# p.vals[p.idx.pairs[best.p,]]
# # beta.med.res[best.p]
# beta.gaussaniso =  beta.med.res[,best.p]


load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSChisq_df3.RData")
# p.vals = p.vals[1:7]
# p.vals = p.vals[7:16]
p.idx.pairs = t(combn(1:length(p.vals),2))
beta.med.res = apply(p.idx.pairs,1,function(p.idx){
  p.idx1=p.idx[1]
  p.idx2=p.idx[2]
  perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
  beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
    log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
  beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
                     function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)]))
  return(beta.med)
})
best.p = which.min(apply(beta.med.res,2,function(vec) abs(median(vec[is.finite(vec)])-0.5)))
p.vals[p.idx.pairs[best.p,]]
# beta.med.res[best.p]
beta.chisq =  beta.med.res[,best.p]

load("/scratch/rcampbell/nice/SimStudyPerimeter_FFT/SSParetoGaussian.RData")
p.idx1=4  #3
p.idx2=5  #4
perim1 = res.df %>% filter(p==p.vals[p.idx1]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
perim2 = res.df %>% filter(p==p.vals[p.idx2]) %>% select(perim) %>% unlist() %>% unname() #%>% sample()
beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
  log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
# grp.len = 50
beta.paretogauss =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
                     function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)],na.rm=T))

############################################################################

# boxplots
# pdf("resultsPerimbetas.pdf",width=8,height=6)
par(mfrow=c(1,1),pty="m")
boxplot(beta.gauss,
        beta.gaussaniso,
        beta.chisq,
        beta.paretogauss,
        ylab = "scaling index",
        names = NA
        )
x_labs = c("Gauss.",
           "Gauss. Aniso.",
           "Chi Sq.",
           "Pareto * Gauss."
)
axis(1, at=1:length(x_labs), labels = x_labs)    # add tick marks
abline(h=c(0,0.5),lty=2)
# dev.off()

############################################################################

# QQ-PP plots, individual, then all together

Z.gauss = (beta.gauss-0.5)/sd(beta.gauss)
beta.prop=0
test = (((beta.gauss-beta.prop)))/sd(beta.gauss)  > qnorm(0.95)
Pr.accept = mean(test==0,na.rm=T)
round(Pr.accept,3)
Pr.reject = mean(test==1,na.rm=T)
round(Pr.reject,3)

pdf("QQbetaMedianPerim_Gauss.pdf",5,5)
par(mfrow=c(1,1),pty="s")
qq.model = sort(Z.gauss)
qq.theoretical = qnorm(c(1:length(qq.model))/(length(qq.model)+1))
plot(qq.theoretical,qq.model,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.theoretical,qq.model)),
     ylim=range(c(qq.theoretical,qq.model)),pch=20)
abline(0,1)
dev.off()

pdf("PPbetaMedianPerim_Gauss.pdf",5,5)
par(mfrow=c(1,1),pty="s")
pp.model = pnorm(qq.model)
pp.theoretical = c(1:length(qq.model))/(length(qq.model)+1)
plot(pp.theoretical,pp.model,
     xlab="Theoretical",ylab="Model",
     xlim=c(0,1),
     ylim=c(0,1),pch=20)
abline(0,1)
dev.off()

Z.gauss.aniso = (beta.gaussaniso-0.5)/sd(beta.gaussaniso)
beta.prop=0
test = (((beta.gaussaniso-beta.prop)))/sd(beta.gaussaniso)  > qnorm(0.95)
Pr.accept = mean(test==0,na.rm=T)
Pr.accept
Pr.reject = mean(test==1,na.rm=T)
Pr.reject

pdf("QQbetaMedianPerim_GaussAnisotropic.pdf",5,5)
par(mfrow=c(1,1),pty="s")
qq.model = sort(Z.gauss.aniso)
qq.theoretical = qnorm(c(1:length(qq.model))/(length(qq.model)+1))
plot(qq.theoretical,qq.model,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.theoretical,qq.model)),
     ylim=range(c(qq.theoretical,qq.model)),pch=20)
abline(0,1)
dev.off()

pdf("PPbetaMedianPerim_GaussAnisotropic.pdf",5,5)
par(mfrow=c(1,1),pty="s")
pp.model = pnorm(qq.model)
pp.theoretical = c(1:length(qq.model))/(length(qq.model)+1)
plot(pp.theoretical,pp.model,
     xlab="Theoretical",ylab="Model",
     xlim=c(0,1),
     ylim=c(0,1),pch=20)
abline(0,1)
dev.off()


Z.chisq = (beta.chisq-0.5)/sd(beta.chisq,na.rm=T)
Z.chisq = Z.chisq[is.finite(Z.chisq)]
beta.prop=0
test = (((beta.chisq[is.finite(beta.chisq)]-beta.prop)))/sd(beta.chisq[is.finite(beta.chisq)])  > qnorm(0.95)
Pr.accept = mean(test==0,na.rm=T)
Pr.accept
Pr.reject = mean(test==1,na.rm=T)
Pr.reject

pdf("QQbetaMedianPerim_Chisq.pdf",5,5)
par(mfrow=c(1,1),pty="s")
qq.model = sort(Z.chisq)
qq.theoretical = qnorm(c(1:length(qq.model))/(length(qq.model)+1))
plot(qq.theoretical,qq.model,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.theoretical,qq.model)),
     ylim=range(c(qq.theoretical,qq.model)),pch=20)
abline(0,1)
dev.off()

pdf("PPbetaMedianPerim_Chisq.pdf",5,5)
par(mfrow=c(1,1),pty="s")
pp.model = pnorm(qq.model)
pp.theoretical = c(1:length(qq.model))/(length(qq.model)+1)
plot(pp.theoretical,pp.model,
     xlab="Theoretical",ylab="Model",
     xlim=c(0,1),
     ylim=c(0,1),pch=20)
abline(0,1)
dev.off()

Z.paretogauss = (beta.paretogauss-0)/sd(beta.paretogauss)
beta.prop=0
test = (((beta.paretogauss-beta.prop)))/sd(beta.paretogauss)  > qnorm(0.95)
Pr.accept = mean(test==0,na.rm=T)
Pr.accept
Pr.reject = mean(test==1,na.rm=T)
Pr.reject

pdf("QQbetaMedianPerim_ParetoGauss.pdf",5,5)
par(mfrow=c(1,1),pty="s")
qq.model = sort(Z.paretogauss)
qq.theoretical = qnorm(c(1:length(qq.model))/(length(qq.model)+1))
plot(qq.theoretical,qq.model,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.theoretical,qq.model)),
     ylim=range(c(qq.theoretical,qq.model)),pch=20)
abline(0,1)
dev.off()

pdf("PPbetaMedianPerim_ParetoGauss.pdf",5,5)
par(mfrow=c(1,1),pty="s")
pp.model = pnorm(qq.model)
pp.theoretical = c(1:length(qq.model))/(length(qq.model)+1)
plot(pp.theoretical,pp.model,
     xlab="Theoretical",ylab="Model",
     xlim=c(0,1),
     ylim=c(0,1),pch=20)
abline(0,1)
dev.off()

####################

# plot all on 1 figure for main body

qq.gauss = sort(Z.gauss)
qq.gauss.aniso = sort(Z.gauss.aniso)
qq.chisq = sort(Z.chisq)
qq.paretogauss = sort(Z.paretogauss)
qq.model.all = list(qq.gauss,
                    qq.gauss.aniso,
                    qq.chisq,
                    qq.paretogauss)
qq.true.all = lapply(qq.model.all,function(qq.model) c(1:length(qq.model))/(length(qq.model)+1))

pdf("QQbetaMedianPerim_ALL.pdf",6,6)
par(mfrow=c(1,1),pty="s")
plot(NULL,
     xlab="Theoretical",ylab="Model",
     xlim=range(c(qq.model.all,qq.true.all)),
     ylim=range(c(qq.model.all,qq.true.all)))
for(plt.i in 1:length(qq.model.all)){
  points(qnorm(qq.true.all[[plt.i]]),qq.model.all[[plt.i]],pch=plt.i,col=plt.i+1)
}
abline(0,1)
legend('bottomright',legend=c("Gauss.",
                              "Gauss. Aniso.",
                              "Chi Sq.",
                              "Pareto * Gauss."
),
       pch=1:length(qq.model.all),col=1+c(1:length(qq.model.all)))
dev.off()

pdf("PPbetaMedianPerim_ALL.pdf",6,6)
par(mfrow=c(1,1),pty="s")
plot(NULL,
     xlab="Theoretical",ylab="Model",
     xlim=c(0,1),
     ylim=c(0,1))
for(plt.i in 1:length(qq.model.all)){
  points(qq.true.all[[plt.i]],pnorm(qq.model.all[[plt.i]]),pch=plt.i,col=plt.i+1)
}
abline(0,1)
legend('bottomright',legend=c("Gauss.",
                              "Gauss. Aniso.",
                              "Chi Sq.",
                              "Pareto * Gauss."
),
       pch=1:length(qq.model.all),col=1+c(1:length(qq.model.all)))
dev.off()