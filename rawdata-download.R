# this file is loaded by each Rnw
# all common variables/functions declarations etc. go here
require(GenomicFeatures)
require(biomaRt)
require(BSgenome)

# workaround to get previous version of mart
mymart79 = function(biomart="ensembl",
                    dataset="hsapiens_gene_ensembl",
                    host="www.biomart.org",
                    port=80) {
  ### Redefine GenomicFeatures useMart2 to use a previous version of assembly
  message("use previous mart version ... ")
  useEnsembl(biomart=biomart, dataset=dataset,version=79)
}
mymartCurrent = function(biomart="ensembl",
                         dataset="hsapiens_gene_ensembl",
                         host="www.biomart.org",
                         port=80) {
  message("use current mart version ... ")
  useEnsembl(biomart=biomart, dataset=dataset)
}
tmpfun = NULL
tmpfun <- get(".useMart2", envir = asNamespace("GenomicFeatures"))
environment(mymart79) = environment(tmpfun)
attributes(mymart79) = attributes(tmpfun)
environment(mymartCurrent) = environment(tmpfun)
attributes(mymartCurrent) = attributes(tmpfun)

# set database genome nomenclature
ens.org = c("hsapiens_gene_ensembl","drerio_gene_ensembl", "mmusculus_gene_ensembl")
names(ens.org) = c("Human","Zebrafish","Mouse")
ucsc.org = c("hg38","danRer10","mm10")
names(ucsc.org) = names(ens.org)

# define some tx classes
longclass=c("lincRNA","sense_intronic","sense_overlapping","antisense")
small = c("rRNA","snRNA","miRNA","snoRNA")

# detect installed genomes
instGenomes=ens.org
for(i in names(ens.org)) {
  assemblyver = ucsc.org[i] #sub("(\\w+)\\_(\\w+)\\_(\\w+)", "\\1", ens.org[i],perl=T)
  inst = grep(assemblyver,installed.genomes(),ignore.case = T)
  if (length(inst)>0) {
    instGenomes[i] = installed.genomes()[inst]
  } else {
    instGenomes[i] = ""
  }
}

# Le coppie di specie omologhe
species = list(HM=c("Human","Mouse"),MZ=c("Mouse","Zebrafish"),HZ=c("Human","Zebrafish"))

# prob random seqs
RNDprobs=list(Human=c(0.41/2,0.41/2,0.59/2,0.59/2),
              Mouse=c(0.42/2,0.42/2,0.58/2,0.58/2),
              Zebrafish=c(0.38/2,0.38/2,0.62/2,0.62/2))


# Kozak consensus sequences PWM
# Natural Variability of Kozak Sequences Correlates with Function in a Zebrafish Model
# Steven J. Grzegorski,  Estelle F. Chiari,  Amy Robbins,  Phillip E. Kish,  Alon Kahana
# PlosONE 2014
kozak=list()
kozak[["Mouse"]] = rbind(A=c(21,19,21,20,18,25,48,29,17,1,0,0,23,28,15),
                         C=c(27,34,31,23,32,38,9,39,47,0,0,0,14,40,26),
                         G=c(34,28,27,39,29,25,36,20,28,0,0,1,49,18,37),
                         T=c(19,19,21,18,20,12,6,13,8,0,1,0,14,14,22))
kozak[["Human"]] = rbind(A=c(20,20,21,21,19,24,46,29,19,1,0,0,22,28,16),
                         C=c(27,33,32,23,32,38,10,38,45,0,0,0,15,40,26),
                         G=c(35,29,28,39,30,26,37,20,28,0,0,1,49,18,37),
                         T=c(18,18,19,17,19,12,7,13,8,0,1,0,14,15,21))
kozak[["Zebrafish"]] = rbind(A=c(25,24,23,24,21,33,58,37,25,1,0,0,32,28,16),
                             C=c(24,31,28,18,32,35,7,29,43,0,0,0,8,34,24),
                             G=c(30,23,19,33,25,17,29,18,25,0,0,1,42,19,38),
                             T=c(22,22,30,26,23,15,6,16,7,0,1,0,18,18,21))
for(k in names(kozak)) {
  s=apply(kozak[[k]],2,sum)
  kozak[[k]]=t(t(kozak[[k]])/s)
  kozak[[k]]=log2(kozak[[k]]/0.25)
}


# Fickett Score look-up tables 
# Recognition of protein coding regions in DNA sequences 
# Fickett, James W.
# Nucleic acids research 1982

Bases = c("A","C","G","T") 

posPfrom = c(0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9)
posPto = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,Inf)
posP = cbind(posPfrom,posPto)

contPfrom = c(0,0.17,0.19,0.21,0.23,0.25,0.27,0.29,0.31,0.33)
contPto = c(0.17,0.19,0.21,0.23,0.25,0.27,0.29,0.31,0.33,0.99)
contP = cbind(contPfrom,contPto)

ProbCoding.posP = rbind(c(0.22,0.23,0.08,0.09),
                        c(0.20,0.30,0.08,0.09),
                        c(0.34,0.33,0.16,0.20),
                        c(0.45,0.51,0.27,0.54),
                        c(0.68,0.48,0.48,0.44),
                        c(0.58,0.66,0.53,0.69),
                        c(0.93,0.81,0.64,0.68),
                        c(0.84,0.70,0.74,0.91),
                        c(0.68,0.70,0.88,0.97),
                        c(0.94,0.80,0.90,0.97))
colnames(ProbCoding.posP) = Bases

ProbCoding.contP = rbind(c(0.21,0.31,0.29,0.58),
                         c(0.81,0.39,0.33,0.51),
                         c(0.65,0.44,0.41,0.69),
                         c(0.67,0.43,0.41,0.56),
                         c(0.49,0.59,0.73,0.75),
                         c(0.62,0.59,0.64,0.55),
                         c(0.55,0.64,0.64,0.40),
                         c(0.44,0.51,0.47,0.39),
                         c(0.49,0.64,0.54,0.24),
                         c(0.28,0.82,0.40,0.28))
colnames(ProbCoding.contP) = Bases

posW = c(0.26,0.18,0.31,0.33)
names(posW) = Bases

contW = c(0.11,0.12,0.15,0.14)
names(contW) = Bases

getProm = function(upseq,gseq,up=2000,dw=500) {
  inizio = width(upseq) - up
  inizio[inizio<=0]=1
  fine = rep(dw,length(gseq))
  fine[width(gseq)<dw]=width(gseq)[width(gseq)<dw]
  rseq=gseq
  rupseq=subseq(upseq,start=inizio,end=width(upseq))
  rgseq=subseq(gseq,start=1,end=fine)
  for(h in 1:length(rseq)) {
    rseq[[h]] = c(rupseq[[h]],rgseq[[h]])
                  
  }
  return(rseq)
}

downloadRaw = function(i="Mouse",mart,txseq=T,geneseq=T,upseq=1000, 
                       promReg=c(2000,500), orthologs="",
                       suppcol=T,skipIfFileExists=F) {
  o=ucsc.org[i]
  message("Starting download ",i," raw data (assembly ",o,")")
  fn=paste0(o,".Rdata")
  
  assignInNamespace(".useMart2", mart, ns="GenomicFeatures")
  txdb <- makeTxDbFromBiomart(biomart="ensembl",dataset=ens.org[i])
  print(txdb)
  seqlevels(txdb) = grep("^\\d+|^X|^Y",seqlevels0(txdb),value=T,perl=T)
  seqlevels(txdb) = paste0("chr",seqlevels(txdb))
    
  genome = NULL
  if (instGenomes[i]!="") {
    require(instGenomes[i],character.only = TRUE)
    genome = get(instGenomes[i])
  } else {
    warning(paste(i,"genome not installed seqs will not be downloaded"))
  }
  
  if (!(skipIfFileExists & file.exists(paste0("txSeqs-",fn)))) {
    if (!is.null(genome) & txseq) {
      txSeqs <- extractTranscriptSeqs(genome, txdb,use.names=T)
      save(txSeqs,file=paste0("txSeqs-",fn))
    }
  }

  if (!(skipIfFileExists & file.exists(paste0("geneRanges-",fn)))) {
    geneRanges = exonsBy(txdb, by="gene",use.names=F)
    save(geneRanges,file=paste0("geneRanges-",fn))
  } else {
    load(paste0("geneRanges-",fn))
  }
  
  
  if (!(skipIfFileExists & file.exists(paste0("txRanges-",fn)))) {
    txRanges <- exonsBy(txdb, by = "tx",use.names=T)
    save(txRanges,file=paste0("txRanges-",fn))
  } else {
    load(paste0("txRanges-",fn))
  }
  
  if (!is.null(genome) & geneseq) {
    if (!(skipIfFileExists & file.exists(paste0("geneSeqs-",fn)))) {
      geneSeqs = extractTranscriptSeqs(genome, geneRanges)
      save(geneSeqs,file=paste0("geneSeqs-",fn))
    }
  }
  
  if (!is.null(genome) & length(upseq)>0) {
    for (us in upseq) {
      if (!(skipIfFileExists & file.exists(paste0("up",us,"Seqs-",fn)))) {
        upSeqs = extractUpstreamSeqs(genome, txdb,width=us)
        save(upSeqs,file=paste0("up",us,"Seqs-",fn))
      }
    }
  }
  
  if (!is.null(genome) & length(promReg)>0) {
    if (!(skipIfFileExists & file.exists(paste0("promSeqs-",fn)))) {
      promRanges = trim(promoters(txdb,upstream = promReg[1],downstream = promReg[2]))
      save(promRanges,file=paste0("promRanges-",fn))
      promSeqs = getSeq(genome,promRanges)
      names(promSeqs) = promRanges$tx_name
      save(promSeqs,file=paste0("promSeqs-",fn))
    }
  }


  txClass=select(txdb, keys=names(txRanges), columns=c("TXNAME","TXTYPE","GENEID"), keytype="TXNAME")
  txClass$Class = "other"
  txClass$Class[grep("pseudogene",txClass$TXTYPE)] = "pseudogene"
  txClass$Class[txClass$TXTYPE %in% longclass] = txClass$TXTYPE[txClass$TXTYPE %in% longclass]
  txClass$Class[txClass$TXTYPE=="protein_coding"] = "PCT"
  
  if(suppcol) {
    # download supplemental columns vegaid
    mym = mart(biomart="ensembl",dataset=ens.org[i])
    #la = listAttributes(mym)
    #orthologs=c("drerio","hsapiens")
    if(length(orthologs)>0) {
      orthocol = paste0(orthologs,"_homolog_ensembl_gene")
      orthotable = getBM(attributes = c("ensembl_transcript_id",orthocol), mart = mym)
      orthotable=aggregate(orthotable,by=list(orthotable$ensembl_transcript_id),paste,collapse=",")
      rownames(orthotable) = as.character(orthotable$ensembl_transcript_id)
      txpresenti = txClass$TXNAME %in% orthotable$ensembl_transcript_id
      for(oi in 1:length(orthologs)) {
        txClass[txpresenti,orthologs[oi]] = orthotable[txClass$TXNAME[txpresenti],orthocol[oi]]
      }
    }
    
    vegatable = getBM(attributes = c("ensembl_transcript_id","ottt","transcript_status"), mart = mym)
    rownames(vegatable) = as.character(vegatable$ensembl_transcript_id)
    txpresenti = txClass$TXNAME %in% vegatable$ensembl_transcript_id
    txClass$VEGAID = ""
    txClass$VEGAID[txpresenti] = vegatable[txClass$TXNAME[txpresenti],"ottt"]
    txClass$STATUS = ""
    txClass$STATUS[txpresenti] = vegatable[txClass$TXNAME[txpresenti],"transcript_status"]
  }
  txClass$Species = i
  
  save(txClass,file=paste0("txClass-",fn))
}


