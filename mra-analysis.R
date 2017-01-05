

activity = function(mexp,cormat,tflist=NULL,tau=0.6) {
  
  mexp.s = scale(mexp)
  if (is.null(tflist)) {
    tflist=rownames(cormat)
  }
  actmat = mexp[tflist,]
  actmat[1:length(actmat)]=0
  pb = txtProgressBar(min=1,max=length(tflist),style=3)
  for(tfi in tflist) {
    postrg = names(cormat[tfi,cormat[tfi,] > tau])
    negtrg = names(cormat[tfi,cormat[tfi,] < -tau])
    apos = 1
    if (length(postrg)>0) {
      apos = apply(mexp.s[postrg,,drop=F],2,sum)/length(postrg)
    }
    aneg = 1
    if (length(negtrg)>0) {
      aneg = apply(mexp.s[negtrg,,drop=F],2,sum)/length(negtrg)
    }
    actmat[tfi,] =  apos - aneg
    setTxtProgressBar(pb, i)
  }
  return(actmat)
}