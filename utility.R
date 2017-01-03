consolidate.genenames <- function(M,genenames,t=0.2) {
  MM_r<-NULL
  nuovinomi<-c()
  for(g in unique(genenames)) {
    print(g)
    #righe<-which(annotations$X.Gene.Symbol==g)
    righe <- which(genenames == g)
    if(length(righe)==1) {
      newrow<-M[righe,]
      nuovinomi<-c(nuovinomi,g)
      MM_r<-rbind(MM_r,newrow)
      next
    }
    cor.dist <- as.dist(1 - cor(t(M[righe, ])))
    ht <- hclust(cor.dist,method = "complete")
    clusters<-cutree(ht,h=t)  # 1- corr
    cl<-unique(clusters)
    if (length(cl)==1) {
      newrow<-M[righe,]
      if (length(righe)>1) {
        newrow<-apply(M[righe,],2,mean)
      }
      nuovinomi<-c(nuovinomi,g)
      MM_r<-rbind(MM_r,newrow)
    } else {
      cat(g,"aggrego",length(righe),"probes in",length(cl),"clusters\n")
      for(i in cl) {
        subrighe<-righe[clusters==i]
        newrow<-M[subrighe,]
        if (length(subrighe)>1) {
          newrow<-apply(M[subrighe,],2,mean)
        }
        nuovinomi<-c(nuovinomi,paste(g,i,sep="!"))
        MM_r<-rbind(MM_r,newrow)
      }
    }
  }
  rownames(MM_r) <- nuovinomi
  return(MM_r)
}

