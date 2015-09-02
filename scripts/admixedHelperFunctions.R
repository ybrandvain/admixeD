# Yaniv work on admixed simulations Aug 27 2015
library(parallel)

rm(list=ls())
ls()
source('~/Desktop/ad/admixedPrimaryFunctions.R')


summarizeSim <- function(geno.data,trim.method, sel.loci = NA, just.estimate = FALSE, trim.intense = .2){
  # geno.data can either be a path to a file, or can be a matrix
  if( class(geno.data)  == "character"){    
    geno.data <- as.matrix(read.csv(geno.data, sep = "\t"))/2 
  }
  if( nrow(geno.data) == 0){return(NULL)}
  locus.pairs <-t(combn(seq_along(geno.data[1,]),2))
  if(trim.method == "none" | trim.method == "regularR2"){    
    admixture.prop <- rowMeans(geno.data)
    removed <- NA
    if(trim.method == "none"){
      pw.sum <- apply(locus.pairs,1,function(PAIR){
        LDcalcs(geno.data[,PAIR[1]], geno.data[,PAIR[2]], admixture.prop, to.return = "partial")
      })
    }
    if(trim.method == "regularR2"){
      pw.sum <- apply(locus.pairs,1,function(PAIR){
        LDcalcs(geno.data[,PAIR[1]], geno.data[,PAIR[2]], admixture.prop, to.return = "regular")
      })}
  }
  if(trim.method == "rm.sel.chr" ){
    sel.chr <- c(grep("group1.",colnames(geno.data),fixed=T), grep("group2.",colnames(geno.data),fixed=T))
    admixture.prop <- rowMeans(geno.data[,-sel.chr])
    removed = NA
    pw.sum <- apply(locus.pairs,1,function(PAIR){
      LDcalcs(geno.data[,PAIR[1]], geno.data[,PAIR[2]], admixture.prop, to.return = "partial")
    })
  }
  if(trim.method != "none" & trim.method != "rm.focal" & trim.method  != "rm.sel.chr" & trim.method != "regularR2"){
    rem <- trimAncestryProp(all.loci=geno.data, trim = trim.intense, method = trim.method)
    admixture.prop <- rem$new.alpha
    removed <- paste(rem$weirdos,collapse="_")
    rm(rem)
    pw.sum <- apply(locus.pairs,1,function(PAIR){
      LDcalcs(geno.data[,PAIR[1]], geno.data[,PAIR[2]], admixture.prop, to.return = "partial")
    })
  }
  if(trim.method == "rm.focal"){
    chrs <- do.call(cbind,strsplit(colnames(geno.data),".",fixed=T))[1,]
    removed.focal <- t(apply(locus.pairs,1,function(PAIR){ rowMeans(geno.data[,!chrs%in%chrs[PAIR]]) }))
    rownames(removed.focal) <- paste(locus.pairs[,1],locus.pairs[,2])
    removed <- NA
    pw.sum <- apply(locus.pairs,1,function(PAIR){
      admixture.prop <- removed.focal[paste(PAIR,collapse=" "),] 
      LDcalcs(geno.data[,PAIR[1]], geno.data[,PAIR[2]], admixture.prop, to.return = "partial")
    })
  }
  c1 <- do.call(rbind,strsplit(colnames(geno.data)[locus.pairs[,1]],".",fixed =T))
  c2 <- do.call(rbind,strsplit(colnames(geno.data)[locus.pairs[,2]],".",fixed =T))
  if(class(pw.sum) == "list"){ pw.sum <- data.frame( do.call(rbind,pw.sum)   )}
  if(class(pw.sum) == "matrix"){ pw.sum <- data.frame(t(pw.sum)) }
  pw.sum <-cbind(c1,c2,pw.sum)
  colnames(pw.sum)[1:4]  <-  paste( rep(c("chrom","loc"),2) ,  rep(c("A","B"), each = 2) , sep = "_")
  pw.sum <- pw.sum[ pw.sum$chrom_A !=pw.sum$chrom_B, ]
  if(just.estimate) {return(pw.sum)}
  #summarize performance  
  if(is.na(sel.loci)){foc <- "none" }
  if(!is.na(sel.loci)){foc <- sapply(sel.loci, function(X){  paste(unlist(strsplit(X,".",fixed = T))[c(1,3)],collapse =" ")  }) }  
  not.our.chr.pairs <- !with(pw.sum,paste(chrom_A,chrom_B) %in% foc)
  arb.loc <- which(with(pw.sum[not.our.chr.pairs,], chrom_A == "group4" & loc_A == 10 & chrom_B == "group5" & loc_B == 5))
  perform <- with(pw.sum[not.our.chr.pairs,],c(
      false.pos.alpha.001 = sum(p.value < 0.001,na.rm=T) / sum(!is.na(p.value)),
      false.pos.alpha.01  = sum(p.value < 0.01,na.rm=T) / sum(!is.na(p.value)),
      false.pos.alpha.05  = sum(p.value < 0.05,na.rm=T) / sum(not.our.chr.pairs,na.rm=T),
      quantiles = quantile(estimate,na.rm=T),
      focal.p = ifelse(length(arb.loc)==0,NA, p.value[arb.loc]),
      focal.estimate = ifelse(length(arb.loc)==0,NA, estimate[arb.loc])
  ))  
  if(!is.na(sel.loci)){
      true.pos <- lapply( sel.loci , function(PAIR){  
        focal.pair <- with(pw.sum, paste(chrom_A,loc_A,sep = ".")%in%PAIR & paste(chrom_B,loc_B,sep = ".")%in%PAIR)
        if(sum(focal.pair) == 0){ return( c(p.val.rank.true.pos = NA, cor.rank.true.pos = NA))}
        with(pw.sum,
             c( p.val.rank.true.pos = sum(p.value <= p.value[focal.pair],na.rm=T) , 
                cor.rank.true.pos   = sum(abs(estimate) >= abs(estimate)[focal.pair],na.rm=T)   ))
      })
      tmp.names <- names(true.pos[[1]])
      true.pos <- Reduce(paste,true.pos)
      names(true.pos) <-  tmp.names
      perform <- c(perform,true.pos, n.comps = sum(!is.na(pw.sum$estimate)))
    }
  return(list( perform=perform, removed = removed ))
}

#setwd("~/Dropbox/LD_lowess_simulations_with_Yaniv/Admix'em_simulation_results/molly/neutral")
#	sapply(list.files(),function(p){
#		runAdmixeD <- function(in.file){
#	  		print(in.file)
#	  		approach <- list(regularR2="regularR2", none = "none", rm_focal = "rm.focal", rm_put_sel = "ancestry")
#	  		sims <- lapply( approach, function(TRIM){ summarizeSim(in.file,TRIM) }) #,sel.loci = list(c("group1.1","group2.1"))) })
#	  		list( perform = do.call(c,lapply(sims,function(X){X[[1]]})) , removed = sims$rm_put_sel$removed) 
#		}
#		path <- paste( p, "/msg_format_subsample_",sep="")
		# change here for mutlicore		 
#		sim.sum <- mclapply(paste( path, 1:500, sep = ""), runAdmixeD,mc.cores=7)
#		sim.perform <- t(sapply(sim.sum,function(X){X[["perform"]]}))
#		IDing.sel <- sapply(sim.sum,function(X){X[["removed"]]})
#		write.csv( sim.perform , file = paste(p,"/","perform.csv",sep="") )
#		write.csv( IDing.sel , file = paste(p,"/","IDingsel.csv",sep="") )
#		p
#	}
#)




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





setwd("~/Dropbox/LD_lowess_simulations_with_Yaniv/Admix'em_simulation_results/molly/selection")
	sapply(list.files()[-c(1:2)],function(p){
		runAdmixeD <- function(foc,sel.loci){
			if(sel.loci == "rand") {
				temp.sel <- read.csv(paste(p, "/selected_pairs_log_",foc,sep=""),sep = "\t",stringsAsFactors=FALSE)
				temp.sel$chr1name <- as.numeric(sapply( strsplit( gsub("chr","",temp.sel$chr1name) , "_") ,function(X){X[1]}))
				temp.sel$chr2name <- as.numeric(sapply( strsplit( gsub("chr","",temp.sel$chr2name) , "_") ,function(X){X[1]}))
				sel.loci <- lapply( 1:nrow(temp.sel),function(X){
					with(temp.sel[X,],c( paste("group",chr1name,".",coordinate1,sep=""), paste("group",chr2name,".",coordinate2,sep=""))[order(c(chr1name,chr2name))])
				})
			}
			in.file <- paste( path, foc, sep = "")
	  		print(in.file)
	  		approach <- list(regularR2="regularR2", none = "none", rm_focal = "rm.focal", rm_put_sel = "ancestry")
	  		sims <- lapply( approach, function(TRIM){ summarizeSim(in.file,TRIM,sel.loci = sel.loci ) })
	  		list( perform = do.call(c,lapply(sims,function(X){X[[1]]})) , removed = sims$rm_put_sel$removed) 
		}
		path <- paste( p, "/msg_format_subsample_",sep="")
		sel <- ifelse( length(grep("selected_pairs_log",list.files(p))) == 0, list(c("group1.1","group2.1")),"rand")
		# change here for mutlicore		 
		sim.sum <- mclapply(1:500, runAdmixeD,sel ,mc.cores=8)
		nulls <- sapply(sim.sum,function(Z){is.null(Z[[1]])})
		sim.perform <- do.call(rbind,sapply(sim.sum[!nulls],function(X){X[["perform"]]},simplify=FALSE))
		IDing.sel <- do.call(rbind,sapply(sim.sum[!nulls],function(X){X[["removed"]]},simplify=FALSE))
		write.csv( sim.perform , file = paste(p,"/","perform.csv",sep="") )
		write.csv( IDing.sel , file = paste(p,"/","IDingsel.csv",sep="") )
		p
	}
)