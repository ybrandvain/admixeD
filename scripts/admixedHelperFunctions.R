# Yaniv work on admixed simulations Aug 21 2015

rm(list=ls())
ls()
source("admixedPrimaryFunctions.R") # MAKE SURE THIS PATH IS CORRECT

summarizeSim <- function(geno.data,trim.method, trim.intense = .2){
  # geno.data can either be a path to a file, or can be a matrix
  if( class(geno.data)  == "character"){    geno.data <- as.matrix(read.csv(geno.data, sep = "\t"))/2 }
  locus.pairs <-t(combn(seq_along(geno.data[1,]),2))
  if(trim.method == "rm.sel.chr" ){
      sel.chr <- c(grep("group1.",colnames(geno.data),fixed=T), grep("group2.",colnames(geno.data),fixed=T))
      admixture.prop <- rowMeans(geno.data[,-sel.chr])
      removed = NA
  }
  if(trim.method == "none" ){    admixture.prop <- rowMeans(geno.data) ; removed = NA}
  if(trim.method != "none" & trim.method != "rm.focal" & trim.method  != "rm.sel.chr"){
    rem <- trimAncestryProp(all.loci=geno.data, trim = trim.intense, method = trim.method)
    admixture.prop <- rem$new.alpha
    removed <- paste(rem$weirdos,collapse="_")
    rm(rem)
  }
  if(trim.method == "rm.focal"){
    chrs <- do.call(cbind,strsplit(colnames(geno.data),".",fixed=T))[1,]
    removed.focal <- t(apply(locus.pairs,1,function(PAIR){ rowMeans(geno.data[,!chrs%in%chrs[PAIR]]) }))
    rownames(removed.focal) <- paste(locus.pairs[,1],locus.pairs[,2])
    removed <- NA
  }
  pw.sum <- apply(locus.pairs,1,function(PAIR){
    if(trim.method == "rm.focal"){admixture.prop <- removed.focal[paste(PAIR,collapse=" "),] }
    LDcalcs(geno.data[,PAIR[1]], geno.data[,PAIR[2]], admixture.prop)
  })
  pw.sum <-cbind(do.call(rbind,strsplit(colnames(geno.data)[locus.pairs[,1]],".",fixed =T)),
      do.call(rbind,strsplit(colnames(geno.data)[locus.pairs[,2]],".",fixed =T)),
      data.frame(t(pw.sum)))
  colnames(pw.sum)[1:4]  <-  paste( rep(c("chrom","loc"),2) ,  rep(c("A","B"), each = 2) , sep = "_")
  pw.sum <- pw.sum[ pw.sum$chrom_A !=pw.sum$chrom_B, ]
  #summarize performance
  p.val.rank.true.pos <- sum(with(pw.sum,p.value <= p.value[1]))
  cor.rank.true.pos <- sum(with(pw.sum,abs(estimate) >= abs(estimate)[1]))
  ld.rank.true.pos <- sum(with(pw.sum,abs(reg.D) >= abs(reg.D)[1]))
  r2.rank.true.pos <- sum(with(pw.sum,abs(reg.R) >= abs(reg.R)[1]))
  means <- colMeans(pw.sum[,c("reg.D","reg.R","estimate")])
  focal <- with(pw.sum,which(chrom_A == "group4" & loc_A == 10 & chrom_B == "group5" & loc_B == 5))
  focal <- unlist(pw.sum[focal,c("reg.D","reg.R","estimate")] )
  quantiles <- apply(pw.sum[,c("reg.D","reg.R","estimate")],2,quantile)
  quant <- c(quantiles)
  names(quant) <- paste( rep(colnames(quantiles) , each = nrow(quantiles)) , rep(rownames(quantiles) , ncol(quantiles)),sep=".")
  #note this type two error quant is 'cooked' for this sim, removing selected chroms
  typeII <- with(pw.sum[pw.sum$chrom_A!="group1" & pw.sum$chrom_B!="group2" ,], 
       sapply (c(0.001,0.01,.05),function(X){sum( p.value < X)  / length(p.value)
  }))
  names(typeII) <- c(0.001,0.01,.05)
  #print("done")
  pw.sum <- pw.sum[with(pw.sum,as.numeric(loc_A)<11 & as.numeric(loc_B)<11),]
  pw.sum <- pw.sum[order(with(pw.sum,paste(chrom_A,chrom_B,loc_A,loc_B))),]
  return(
    c(est.rank.true = cor.rank.true.pos, 
    ld.rank.true = ld.rank.true.pos,
    r2.rank.true = r2.rank.true.pos,
    means = means,
    focal = focal,
    quant = quant,
    typeII = typeII,
    removed = removed)
  )
}
near.path <- "/Users/ybrandva/Dropbox/LD_lowess_simulations_with_Yaniv/Admix\'em_simulation_results/yb/pulse/msg_format_subsample" 


runAdmixeD <- function(in.file){
  none <- summarizeSim(in.file,trim.method = "none")
  rm.focal <- summarizeSim(in.file,trim.method = "rm.focal")
  rm.put.sel <- summarizeSim(in.file,trim.method = "ancestry")
  
  print(in.file)
  
  # Formating output
  all.init <- none[c( "ld.rank.true", "r2.rank.true", "means.reg.D", "means.reg.R", "focal.reg.D", "focal.reg.R" , "quant.reg.D.0%" , "quant.reg.D.25%", "quant.reg.D.50%" ,"quant.reg.D.75%","quant.reg.D.100%", "quant.reg.R.0%" , "quant.reg.R.25%", "quant.reg.R.50%" ,"quant.reg.R.75%","quant.reg.R.100%" )]
  all.others <- do.call(c,lapply( list( none = none,  rm.focal= rm.focal,rm.put.sel=rm.put.sel), function(X){
    X[!names(X)%in%names(all.init)]
  }))
  return(c(all.init,all.others))
}


# change here for mutlicore
sim.sum <- t(sapply(paste(near.path,1:5,sep="_"), runAdmixeD))
  rownames(sim.sum) <- NULL
  sim.sum <- data.frame(apply(sim.sum[,-ncol(sim.sum)],2,as.numeric), rm.put.sel.removed = sim.sum[,ncol(sim.sum)]) 
  






# I tried a few different approaches. 
  # "rm.sel.chr" removes the first two chroms when calculating ancestry proportions 
  # "rm.focal" removes focal chromsomes when calculating ancestry proportions
  # "none" calculates admixture propostions from every marker
  # "lowcor" removes markers that poorly predict ancestry proportion when calculating ancestry proportion
  # "cline" removes markers whos genotypes are poorly predicted by the 'genomic cline' fit for that locus  when calculating ancestry proportion
  # "hit" removes markers with too many homozygotes given ancetry prop  when calculating ancestry proportion
  # "ancestry" removes markers that are poorly predicted by ancestry prop when calculating ancestry proportion

# In practice, "lowcor" &  "cline" are about the same, and "hit" and "ancestry" are about the same. with the latter doing a much better job of tagging selected sites. So I am now just looking at none, rm.focal, and ancestry 
# In this implementation (on admixsimul.cfg) "rm.focal" appears to do best.






