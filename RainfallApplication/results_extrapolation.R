# Gather results of the LER and perimeters

rm(list=ls())
library(this.path)
library(SpatialExtremes)
library(fields)
library(geoR)
library(this.path)
setwd(this.path::here())
library(ExcursionSets)
library(terra)
library(dplyr)
library(quantreg)

beta.vals.qr = c()
par(mfrow=c(1,3))
for(loc.idx in 1:3){
  load(paste0("LER_loc", loc.idx, ".RData"))
  QR.dat.df = QR.dat.df %>% filter(p>0.9)
  fit1 = rq(log.lers~nl.nl.p, data=QR.dat.df)
  beta.vals.qr = c(beta.vals.qr,coefficients(fit1)[2])
}

print("LER-based values")
print(round(beta.vals.qr,2))

beta.vals.perim = c()
for(loc.idx in 1:3){
  load(paste0("perims_loc",loc.idx,".RData"))
  p.vals = sort(unique(res.df$p))
  p.idx.pairs = t(combn(1:length(p.vals),2))
  res.split = list(res.df)
  
  beta.est.lst = lapply(res.split, function(lst){
    lst = lst %>% select(-window.size,-obs.idx)
    beta.med.res = apply(p.idx.pairs,1,function(p.idx){
      p.idx1=p.idx[1]
      p.idx2=p.idx[2]
      perim1 = lst %>% filter(p==p.vals[p.idx1]) %>% select(perim_biased) %>% unlist() %>% unname() #%>% sample()
      perim2 = lst %>% filter(p==p.vals[p.idx2]) %>% select(perim_biased) %>% unlist() %>% unname() #%>% sample()
      beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
        log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
      return(median(beta.est[is.finite(beta.est)]))
    })
    return(cbind(beta.val=beta.med.res))
  })
  beta.vals.perim = c(beta.vals.perim,median(unlist(beta.est.lst)))
}

print("Perimeter-based values")
print(round(beta.vals.perim,2))

# block bootstrap (LER)
n.bs = 200
beta.ler.bs = matrix(NA,nrow=n.bs,ncol=3)
seed.val=2112
for(bs.rep in 1:n.bs){
  beta.vals = c()
  for(loc.idx in 1:3){
    load(paste0("LER_loc", loc.idx, ".RData"))
    
    block_length <- 72
    
    # Unique hourly indices
    hours <- sort(unique(res.df$obs.idx))
    n_hours <- length(hours)
    boot_hours <- numeric(0)
    set.seed(seed.val)
    seed.val = seed.val+1
    while (length(boot_hours) < n_hours) {
      # Random starting hour
      start <- sample(hours, 1)
      # Hours in the block: start through 72 hours later
      block <- hours[
        hours >= start &
          hours < start + block_length
      ]
      # Append block
      boot_hours <- c(boot_hours, block)
    }
    
    # Keep exactly the same number of hours
    boot_hours <- boot_hours[seq_len(n_hours)]
    
    # Generate bootstrap sequence of hours
    boot.df <- bind_rows(
      lapply(seq_along(boot_hours), function(i) {
        res.df %>%
          filter(obs.idx == boot_hours[i]) %>%
          mutate(boot.idx = i)
      })
    )
    
    boot.df = boot.df %>% mutate(log.ler = log(ler), nl.nl.p = -log(-log(1-p))) %>% filter(is.finite(log.ler))
    fit1 = quantreg::rq(log.ler~nl.nl.p, data=boot.df)
    beta.ler.bs[bs.rep,loc.idx]=unname(coefficients(fit1)[2])
  }
}

print("LER block bootstrap CI:")
print(round(apply(beta.ler.bs,2,quantile,p=c(0.025,0.975)),2))

# block bootstrap (perimeter)
n.bs = 200
beta.perim.bs = matrix(NA,nrow=n.bs,ncol=3)
seed.val=2112
for(bs.rep in 1:n.bs){
  beta.vals = c()
  for(loc.idx in 1:3){
    load(paste0("perims_loc",loc.idx,".RData"))

    p.vals = sort(unique(res.df$p))
    p.idx.pairs = t(combn(1:length(p.vals),2))
    res.split = list(res.df)
    
    block_length <- 72
    
    # Unique hourly indices
    hours <- sort(unique(res.df$obs.idx))
    n_hours <- length(hours)
    boot_hours <- numeric(0)
    set.seed(seed.val)
    seed.val = seed.val+1
    while (length(boot_hours) < n_hours) {
      # Random starting hour
      start <- sample(hours, 1)
      # Hours in the block: start through 72 hours later
      block <- hours[
        hours >= start &
          hours < start + block_length
      ]
      # Append block
      boot_hours <- c(boot_hours, block)
    }
    
    # Keep exactly the same number of hours
    boot_hours <- boot_hours[seq_len(n_hours)]
    
    # Generate bootstrap sequence of hours
    boot.df <- bind_rows(
      lapply(seq_along(boot_hours), function(i) {
        res.df %>%
          filter(obs.idx == boot_hours[i]) %>%
          mutate(boot.idx = i)
      })
    )
    
    res.split = list(boot.df)
    
    # 
    p.vals = sort(unique(boot.df$p))
    p.idx.pairs = t(combn(1:length(p.vals),2))
    # p.idx.pairs = p.idx.pairs[apply(p.idx.pairs,1,function(vec) all(vec>5)),]
    # p.idx.pairs = p.idx.pairs[apply(p.idx.pairs,1,function(vec) all(vec>2 & vec < 7)),]
    beta.est.lst = lapply(res.split, function(lst){
      lst = lst %>% select(-window.size,-obs.idx)
      beta.med.res.biased = apply(p.idx.pairs,1,function(p.idx){
        p.idx1=p.idx[1]
        p.idx2=p.idx[2]
        perim1 = lst %>% filter(p==p.vals[p.idx1]) %>% select(perim_biased) %>% unlist() %>% unname() #%>% sample()
        perim2 = lst %>% filter(p==p.vals[p.idx2]) %>% select(perim_biased) %>% unlist() %>% unname() #%>% sample()
        beta.est = (log(perim1) - log(perim2) - log((1-p.vals[p.idx1])/(1-p.vals[p.idx2])))/
          log(log(1-p.vals[p.idx1])/log(1-p.vals[p.idx2]))
        # grp.len = 50
        # beta.med =  sapply((1:floor(length(beta.est)/grp.len))*grp.len - grp.len+1,
        #                    function(start.idx) median(beta.est[start.idx:(start.idx+grp.len)]))
        return(median(beta.est[is.finite(beta.est)]))
      })
      
      return(beta.med.res.biased)
    })
    
    beta.perim.bs[bs.rep,loc.idx]=median(unlist(beta.est.lst)[is.finite(unlist(beta.est.lst))])
  }
}

print("Perimeter block bootstrap CI:")
print(round(apply(beta.perim.bs,2,quantile,p=c(0.025,0.975)),2))

#############################################################

par(mfrow=c(1,3))
for(loc.idx in 1:3){
  load(paste0("LER_loc",loc.idx,".RData"))
  p.vals = sort(unique(res.df$p))
  p.vals = p.vals[p.vals>=0.95]  # we want to extrapolate to the far tails
  
  p2 = p.vals[1]
  p1 = p.vals[2:length(p.vals)]
  
  x.axis.emp = p1/p2
  theta.emp = res.df %>% group_by(p) %>% summarise(theta=median(ler[ler>0])) %>% select(theta) %>% unlist() %>% unname()
  LER.ratio.emp = theta.emp[2:length(theta.emp)] / theta.emp[1]
  p1.seq = seq(min(p1),0.9999,length.out=1000)
  beta.val.qr = beta.vals.qr[loc.idx]
  beta.val.perim = beta.vals.perim[loc.idx]
  LER.ratio.est.qr = (log(1-p1.seq)/log(1-p2))^(-beta.val.qr)
  LER.ratio.est.perim = (log(1-p1.seq)/log(1-p2))^(-beta.val.perim)
  
  LER.ratio.est.qr.bs = lapply(beta.ler.bs[,loc.idx], function(bb){
    (log(1-p1.seq)/log(1-p2))^(-bb)
  })
  LER.ratio.est.qr.bs = do.call(rbind,LER.ratio.est.qr.bs)
  LER.ratio.est.qr.bs.low = apply(LER.ratio.est.qr.bs,2,quantile,p=0.025)
  LER.ratio.est.qr.bs.up = apply(LER.ratio.est.qr.bs,2,quantile,p=0.975)
  
  LER.ratio.est.perim.bs = lapply(beta.perim.bs[,loc.idx], function(bb){
    (log(1-p1.seq)/log(1-p2))^(-bb)
  })
  LER.ratio.est.perim.bs = do.call(rbind,LER.ratio.est.perim.bs)
  LER.ratio.est.perim.bs.low = apply(LER.ratio.est.perim.bs,2,quantile,p=0.025)
  LER.ratio.est.perim.bs.up = apply(LER.ratio.est.perim.bs,2,quantile,p=0.975)
  
  plot(NULL,
       xlim=range(c(x.axis.emp,p1.seq/p2)),
       ylim=range(c(LER.ratio.emp,LER.ratio.est.qr,LER.ratio.est.perim)),
       xlab=bquote(p[2]/p[1]*","~p[1]==.(format(p1,digits=2))),
       ylab=expression(q[alpha](R[s[0]]^{(p[2])}) / q[alpha](R[s[0]]^{(p[1])})),
       cex.axis=1.2,cex.lab=1.6)

  CI.x = c(p1.seq/p2,rev(p1.seq/p2))
  CI.x[length(CI.x)+1] = CI.x[1]
  CI.y = c(LER.ratio.est.qr.bs.low,rev(LER.ratio.est.qr.bs.up))
  CI.y[length(CI.y)+1] = CI.y[1]
  CI.y[is.na(CI.y)] = 0
  polygon(x=CI.x,y=CI.y,col=scales::alpha("blue",0.2),border=NA)
  lines(p1.seq/p2,LER.ratio.est.qr,col="blue",lty=2,lwd=2)
  
  CI.x = c(p1.seq/p2,rev(p1.seq/p2))
  CI.x[length(CI.x)+1] = CI.x[1]
  CI.y = c(LER.ratio.est.perim.bs.low,rev(LER.ratio.est.perim.bs.up))
  CI.y[length(CI.y)+1] = CI.y[1]
  CI.y[is.na(CI.y)] = 0
  polygon(x=CI.x,y=CI.y,col=scales::alpha("red",0.2),border=NA)
  lines(p1.seq/p2,LER.ratio.est.perim,col="red",lty=3,lwd=2)
  # 
  points(x.axis.emp,LER.ratio.emp)
  legend(x="topright",legend=c("LER","perim"),lty=c(2,3),col=c("blue","red"),lwd=rep(2,2))
}

###############################################################

# Seasonal model beta values

beta.vals.qr.season.bs = matrix(0,200,3)
for(loc.idx in 1:3){
  load(paste0("LER_loc", loc.idx, ".RData"))
  fit1 = rq(log.ler~nl.nl.p+season, data=res.df)
  print(round(coefficients(fit1),2))
  
  for(i in 1:dim(beta.vals.qr.season.bs)[1]){
    bs.idx = sample(1:dim(res.df)[1],dim(res.df)[1],res.df=T)
    fit.bs = rq(log.lers~nl.nl.p+season, data=res.df[bs.idx,])
    beta.vals.qr.season.bs[i,] = coefficients(fit.bs)
  }
  print(round(apply(beta.vals.qr.season.bs,2,quantile,p=c(0.025,0.975)),2))
}
