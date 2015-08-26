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

all.ld <- data.frame(all.ld)
with( all.ld, min(A,a) )
with( all.ld, plot(  apply(cbind(A,a),1,min) / (A+a) ,  apply(cbind(B,b),1,min) / (B+b)  , xlim = c(0,.5), ylim = c(0,.5) , 
                    col = ifelse(p.value > .05, "black", ifelse(estimate > 0 , "red", " blue")  ),
                    pch=ifelse(p.value > .05,NA,1)))

high.sig <- all.ld[which(all.ld [,"p.value" ]< 0.001),]
high.sig <- data.frame(high.sig, sign = sign(high.sig[,"estimate"]), do.call(rbind,strsplit(rownames(high.sig),"_")))
gene.names <- unlist(do.call(c, sapply( all.processes.chr , names)))
tmp<-rev(seq_along(gene.names)); names(tmp) =gene.names
high.sig$y0 <- tmp[high.sig[,"X1"]]
high.sig$y1 <- tmp[high.sig[,"X2"]]
plot(0,xlim = c(0,3), ylim = range(seq_along(gene.names)), type = "n")
text(x = .5, y = rev(seq_along(gene.names)), gene.names , cex = .5, adj = c(1,0)  )
text(x = 2, y = rev(seq_along(gene.names)), gene.names , cex = .5, adj = c(0,0)  )
with(high.sig,segments(x0 = .52, y0 = y0 + .25, x1 = 1.98, y1 =y1+.25 , col = ifelse(sign == -1, "red", "blue") ))





moreLDcalcs <- function(l1,l2,a, method = "pearson"){
  ok <- !is.na(l1$l + l2$l)
  these.genos <- cbind(A = l1$l, B = l2$l)[ok,]
  p.cor = pcor.test( these.genos[,1], these.genos[,2],a[ok])[1:3]
  my1 <- l1$sim.l[ok,sample(seq_along( l1$sim.l[1,]),replace = T)]
  my2 <- l2$sim.l[ok,sample(seq_along( l2$sim.l[1,]),replace = T)]
  mya <- a[ok]
  sim.dist <- lapply( seq_along(my1[1,]), function(A){
    pcor.test( my1[,A], my2[,A],mya,method = method)[1:3]
  }  )
  sim.dist <- do.call(rbind,sim.dist)
  list( obs = p.cor,   sim.dist  =  sim.dist)
}

all.ld.sim <- do.call(c,lapply(seq_along(chr.combos[,1]),function(THIS){
  COMBO <- chr.combos[THIS,]
  print(COMBO)
  all.pw.comps <- do.call(c,lapply( all.processes.chr[[COMBO[1]]] ,function(LA){
      lapply( (all.processes.chr[[COMBO[2]]]) ,function(LB){
      moreLDcalcs(l1 = LA, l2 = LB, a = hybridindex, method = "spearman")
    })
  }))
  return(all.pw.comps)
}))


head(all.ld.sim)
