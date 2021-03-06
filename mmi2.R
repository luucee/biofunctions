cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)     
  cd  <- md/csd        
  return(cd)
}

aracne2 = function(mexp,from,to,nperm=1000,DPI=TRUE) {
  require(parmigene)
  pb = txtProgressBar(min=1,max=nperm,style=3)
  mi = knnmi.cross(mexp[from,],mexp[to,])
  mi.perm = array(0,dim=c(nrow(mi),ncol(mi)))
  rates = array(0,dim=c(nrow(mi),ncol(mi)))
  for (i in 1:nperm) {
    mi.tmp = knnmi.cross(mexp[from,],mexp[to,sample(1:ncol(mexp))])
    mi.perm = mi.perm + (mi < mi.tmp)
    rates = rates + mi.tmp
    setTxtProgressBar(pb, i)
  }
  
  pval = mi.perm/nperm # higher tail
  rates = nperm/rates
  
  mi0 = mi[pval==0]
  ra0 = rates[pval==0]
  pval[pval==0] = sapply(1:length(mi0),function(i) pexp(mi0[i],rate=ra0[i],lower.tail = F))
  
  #pval.fdr = p.adjust(pval,method = "fdr")
  #mi[pval.fdr<=0.05] = 0
  
  if (DPI) {
    commong = intersect(from,to)
    onlyfrom = setdiff(from,to)
    onlyto = setdiff(to,from)
    extn = length(commong)+length(onlyfrom)+length(onlyto)
    extmi = matrix(0,nrow=extn,ncol=extn)
    rownames(extmi) = c(commong,onlyto,onlyfrom)
    colnames(extmi) = c(commong,onlyto,onlyfrom)
    extmi[commong,commong] = mi[commong,commong,drop=FALSE]
    extmi[onlyto,commong] = t(mi[commong,onlyto,drop=FALSE])
    extmi[onlyfrom,commong] = mi[onlyfrom,commong,drop=FALSE]
    extmi[commong,onlyfrom] = t(mi[onlyfrom,commong,drop=FALSE])
    extmi[onlyto,onlyfrom] = t(mi[onlyfrom,onlyto,drop=FALSE])
    extmi[commong,onlyto] = mi[commong,onlyto,drop=FALSE]
    extmi[onlyfrom,onlyto] = mi[onlyfrom,onlyto,drop=FALSE]
    extmi = aracne.a(extmi)
    mi = extmi[from,to]
  }
  close(pb)
  return(list(MI=mi,PVAL=pval))
}



cross.cor <- function(x, y, verbose = TRUE, ncore="all", met="pearson", ...){
  require(doParallel)
  if(ncore=="all"){
    ncore = parallel:::detectCores()
    doParallel:::registerDoParallel(cores=ncore)
  } else{
    doParallel:::registerDoParallel(cores=ncore)
  }
  N <- nrow(x)
  M <- nrow(y)
  
  ntasks <- ncore
  TM<-round(sqrt(ntasks*M/N))
  TN<-ntasks %/% TM
  if (TM==0) TM=1
  if (TN==0) TN=1
  
  Nsize <- N %/% TN
  Msize <- M %/% TM
  corMAT<-foreach(i = 1:TN, .combine='rbind') %:% 
    foreach(j = 1:TM, .combine='cbind') %dopar% {
      s1<-(i-1)*Nsize+1
      e1<-s1+Nsize-1
      s2<-(j-1)*Msize+1
      e2<-s2+Msize-1     
      if(i==TN) {
        e1<-N
      }
      if(j==TM) {
        e2<-M
      }
      #cat(s1,e1,s2,e2,"\n")
      cor(t(x[s1:e1,]), t(y[s2:e2,]), method = met, use="pairwise.complete.obs", ...)
    }
  gc()
  return(corMAT)
}



mi.ksampled = function(m,r1,r2,k,nboot=100,method="MI") {
  ptm = proc.time()[3]
  l1=length(r1)
  l2=length(r2)
  xnboot = round(sqrt(nboot))
  m1 = array(0,dim=c(l1*xnboot,k))
  m2 = array(0,dim=c(l2*xnboot,k))
  for (i in 1:xnboot) {
    ksample = sample(1:ncol(m),k)
    m1[1:l1 + l1*(i-1),] = m[r1,ksample]
    m2[1:l2 + l2*(i-1),] = m[r2,ksample]
  }
  if (method=="MI") {
    require(parmigene)
    mi = knnmi.cross(m1,m2)
  } else if (method=="pearson") {
    mi = abs(cross.cor(m1,m2,met="pearson"))
    mi[is.na(mi)] = 0
  } else if (method=="spearman") {
    mi = abs(cross.cor(m1,m2,met="spearman"))
    mi[is.na(mi)] = 0
  }
  print(paste0(" mi.ksampled ",proc.time()[3]-ptm," sec."))
  return(mi)
}

mi.split = function(mi,r1,r2) {
  ptm = proc.time()[3]
  l1=length(r1)
  l2=length(r2)
  nboot = nrow(mi) %/% l1
  mi.boot = array(0,dim=c(l1,l2,nboot),dimnames=list(r1,r2))
  mi.perm = array(0,dim=c(l1,l2,nboot**2-nboot),dimnames=list(r1,r2))
  ni =1
  for(i in 1:(nboot-1)) {
    mi.boot[,,i] = mi[1:l1 +l1*(i-1),1:l2 + l2*(i-1)]
    for(j in (i+1):nboot) {
      mi.perm[,,ni] = mi[1:l1 +l1*(i-1),1:l2 + l2*(j-1)]
      ni = ni + 1
      mi.perm[,,ni] = mi[1:l1 +l1*(j-1),1:l2 + l2*(i-1)]
      ni = ni + 1
    }
  }
  mi.boot[,,nboot] = mi[1:l1 +l1*(nboot-1),1:l2 + l2*(nboot-1)]
  mi = apply(mi.boot,c(1,2),mean)

  bb = mean(mi.perm[1:length(mi.perm)])
  cc = sqrt(var(mi.perm[1:length(mi.perm)]))

  pval1 = apply(mi[1:length(mi)] < mi.perm,c(1,2),sum)/(nboot**2-nboot) # higher tail
  pval2 = apply(mi[1:length(mi)] > mi.perm,c(1,2),sum)/(nboot**2-nboot) # lower tail

  pval1[pval1==0] = pnorm(mi[pval1==0], mean = bb, sd = cc, lower.tail = F)
  pval2[pval2==0] = pnorm(mi[pval2==0], mean = bb, sd = cc, lower.tail = T)        

  pval2[pval2>pval1] = pval1[pval2>pval1] 
  #pval1[pval1>pval2] = pval2[pval1>pval2]
  
  return(list(DELTA=mi,PVAL=pval2))
}

mindy2 = function(mexp,mod,tf,target,nbins=5,h=0,nboot=100,perm=F,siglev=0.05,method="MI",verbose=T,scale=F) {
  if(verbose) {
    print(paste0("Exp Matrix: ",paste0(dim(mexp),collapse="x")))
    print(paste0("Modulator: ",mod))
    print(paste0("N° targets: ",length(target)))
    print(paste0("N° TF: ",length(tf)))
    print(paste0("N° bins: ",nbins))
    print(paste0("N° intra-bins to leave: ",h))
    print(paste0("N° boots: ",nboot))
    print(paste0("Sig level: ",siglev))
    print(paste0("Method: ",method))
  }
  
  dbin = ncol(mexp) %/% nbins
  kordering = order(mexp[mod,],decreasing = T)
  if (perm) {
    kordering = sample(kordering)
  }
  
  mii=list()
  for(i in 1:(nbins-1-h)) {
    #ksmpli = 1:(dbin*i)
    ksmpli = (dbin*(i-1)+1):(dbin*i)
    ks=dbin
    if (method=="pearson") {
      ks = length(ksmpli)
    }
    mii[[i]] = mi.ksampled(mexp[,kordering[ksmpli]],tf,target,k = ks,nboot = nboot,method=method)
  }
  mij=list()
  for(j in (2+h):nbins) {
    #ksmplj = (dbin*(j-1)+1):ncol(mexp)
    ksmplj = (dbin*(j-1)+1):(dbin*j)
    ks=dbin
    if(method=="pearson") {
      ks = length(ksmplj)
    }
    mij[[j-1]] = mi.ksampled(mexp[,kordering[ksmplj]],tf,target,k = ks,nboot = nboot,method=method)
  }
  
  retval=NULL
  for(i in 1:(nbins-1-h)) {
    for(j in (i+1+h):nbins) {
      dij = mi.split(mii[[i]] - mij[[j-1]],tf,target)
      trgnames = rep(target,times=rep(length(tf),length(target)))
      tfnames = rep(tf,length(target))
      retval=rbind(retval,data.frame(MOD=mod,TF=tfnames,TRG=trgnames,I=i,J=j,
                              DELTA=dij$DELTA[1:length(dij$DELTA)],PVAL=dij$PVAL[1:length(dij$PVAL)]))
    }
  }
  if(is.numeric(siglev)) {
    #retval$PVAL = p.adjust(retval$PVAL,method = "fdr")
    retval = subset(retval,PVAL<siglev)
  }
  return(retval)
}


cut.outliers = function(mexp) {
  for(i in 1:nrow(mexp)){
    mexp[i,] = (mexp[i,]-mean(mexp[i,]))/sd(mexp[i,],na.rm = T)
  }
  mexp[mexp <= quantile(mexp, 0.05,na.rm = T)] = quantile(mexp, 0.05,na.rm = T)
  mexp[mexp >= quantile(mexp, 0.95,na.rm = T)] = quantile(mexp, 0.95,na.rm = T)
  return(mexp)
}

plot.mod = function(mexp,mod,tf,target,nettarget,fus="",high=TRUE,ordmod=order(mexp[mod,],decreasing=T)) {
  require(amap)
  # ordino per tf
  mexp=mexp[,order(mexp[tf,])]
  crp = colorRampPalette(c("red","white","green"))
  colortf = crp(100)[cut(mexp[tf,],100,labels=F)]
  
  # ordino per mod
  colortf = colortf[ordmod]
  mexp = mexp[,ordmod]
  sbinhigh = ncol(mexp) %/% 3
  sbinlow =  ncol(mexp) %/% 3

  # righe da mostrare
  righe = unique(c(tf,target,nettarget))
  
  # prendo le code ordinate per tf
  mhigh = mexp[righe,1:(sbinhigh)]
  mlow = mexp[righe,(ncol(mexp)-sbinlow):ncol(mexp)]
  colhigh = colortf[1:(sbinhigh)]
  collow = colortf[(ncol(mexp)-sbinlow):ncol(mexp)]
  ordtfhigh = order(mhigh[tf,])
  ordtflow = order(mlow[tf,])
  
  colrighe = rep("blue",length(righe))
  colrighe[righe %in% nettarget] = "red"
  colrighe[righe %in% tf] = "yellow"
  colrighe = colrighe[!righe %in% tf]
  
  
  require(gplots)
  colcolhigh1 = rep("not Fused",ncol(mhigh))
  colcolhigh1[colnames(mhigh) %in% fus] = "Fused"
  colcollow1 = rep("not Fused",ncol(mlow))
  colcollow1[colnames(mlow) %in% fus] = "Fused"
  colcolhigh1=as.factor(colcolhigh1)
  colcollow1=as.factor(colcollow1)
  if (fus!="") {
  tfhigh=data.frame(colcolhigh1[ordtfhigh],mhigh[tf,ordtfhigh])
  colnames(tfhigh)=c("TACC3-FGFR3",tf)
  tflow=data.frame(colcollow1[ordtflow],mlow[tf,ordtflow])
  colnames(tflow)=c("TACC3-FGFR3",tf)
  } else {
    tfhigh=data.frame(mhigh[tf,ordtfhigh])
    colnames(tfhigh)=c(tf)
    tflow=data.frame(mlow[tf,ordtflow])
    colnames(tflow)=c(tf)
  }
  mlow = mlow[!righe %in% tf,]
  mhigh = mhigh[!righe %in% tf,]
  print(dim(mhigh))
  print(dim(mlow))
  cnn=colnames(mlow[,ordtflow])
  print(cnn[colcollow1[ordtflow]=="Fused"])
  require(heatmap3)
  if(high) {
    rownames(mhigh)=rep("",nrow(mhigh))
    colnames(mhigh)=rep("",ncol(mhigh))
    heatmap3(mhigh[,ordtfhigh],col=bluered(100),Colv=NA,
              #RowSideColors=colrighe,
              #ColSideColors=colcolhigh1[ordtfhigh],
              ColSideWidth=1,
              ColSideFun=function(x) showAnn(x),ColSideAnn = tfhigh,
             showRowDendro=F,method="ward.D2",legendfun=function(x){image(matrix(0),col=0,axes=F)},
              scale="none",main="Highest")
  } else {
    rownames(mlow)=rep("",nrow(mlow))
    colnames(mlow)=rep("",ncol(mlow))
    heatmap3(mlow[,ordtflow],col=bluered(100),Colv=NA,
             #RowSideColors=colrighe,
             #ColSideColors=colcollow1[ordtflow],
             ColSideWidth=1,
             ColSideFun=function(x) showAnn(x),ColSideAnn = tflow,
             showRowDendro=F,method="ward.D2",legendfun=function(x){image(matrix(0),col=0,axes=F)},
             scale="none",main="Lowest")
  }
  #return(ic)
}
