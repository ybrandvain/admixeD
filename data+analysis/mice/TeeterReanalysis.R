# Reanalysis of Teeter et al with Admixe'd 
# Yaniv Brandvain and Molly Schumer
# 6/5/15
setwd("/Users/ybrandva/Desktop/projects/admixeD/data+analysis/mice")
#Make sure you're directory is the mouse data
#teeter.data              <- read.csv("EVO_846_sm_TableS1.csv",head=TRUE)
teeter.data              <- read.csv("TeeterSuppTable4.csv",head=TRUE)
geno.cols                <- grep("chr", names(teeter.data))
teeter.data <- teeter.data[teeter.data$Locality == "Neufahrn, Germany",]

teeter.data[, geno.cols] <-apply(teeter.data[,geno.cols],2,function(locus){
  locus[locus %in% c("MD","DM")]       <- .5
  locus[locus == "DD"]                 <-  0 
  locus[locus == "MM"]                 <-  1
  locus[locus %in% c(".","","M","D")]  <- NA # Note, I am treating X chroms as missing data for dudes
  return(as.numeric(locus))
})

hybridindex <- rowMeans(teeter.data[, geno.cols] ,na.rm =T) # caclulating individual admixture props
hist(hybridindex)

#extreme.acest <- hybridindex > .5
#middle.acest <- hybridindex > .3 &  hybridindex < .7
#hybridindex  <- hybridindex[middle.acest] 
#teeter.data  <- teeter.data[middle.acest,]
#hist(hybridindex)

loc.chr <- sapply(strsplit( names(teeter.data)[geno.cols], ".", fixed=T),function(X){X[1]})
all.chr <- as.list(unique(loc.chr))
names(all.chr) <- unlist(all.chr)



source("../../scripts/admixedPrimaryFunctions.R")

prepChr <- function(this.chr.genos, hybridindex){
  apply(this.chr.genos, 2,  getCline, hybridindex,  reps = 1000, return.processed = TRUE)
}

all.processes.chr <- lapply(all.chr, function(CHR){
  print(CHR)
  these.cols <- geno.cols[loc.chr == CHR]
  these.loci <- teeter.data[,these.cols,drop = FALSE]
  prepChr(these.loci, hybridindex = hybridindex)
})


n.chr <- names(all.processes.chr)
chr.combos <- t(combn( x = names(all.chr), 2 ))
rownames(chr.combos) <-  paste(chr.combos[,1],chr.combos[,2],sep="_")
  
all.ld <- do.call(rbind,apply(chr.combos,1,function(COMBO){
  print(COMBO)
  all.pw.comps = do.call(rbind,lapply( all.processes.chr[[COMBO[1]]] ,function(LA){
    do.call(rbind,lapply( (all.processes.chr[[COMBO[2]]]) ,function(LB){
      c(LDcalcs(l1 = LA, l2 = LB, a = hybridindex))      
    }))
  }))
  rownames(all.pw.comps) <- 
    paste(rep(names(all.processes.chr[[COMBO[1]]]),each = length(names(all.processes.chr[[COMBO[2]]]))),
          rep(names(all.processes.chr[[COMBO[2]]]),times = length(names(all.processes.chr[[COMBO[1]]]))),sep="_")
  return(all.pw.comps)
}))


all.ld.filtered <- all.ld[rowSums(all.ld[,c("A","B","a","b")]<10) == 0,]

