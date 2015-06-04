# Yaniv Brandvain & Molly Schumer 6/4/15
# Helpfull functions for admixe'D
# Used to load and write large files, loop across all pairwise comparisons etc. 
# YB is still thinking about them. perhaps not suitable for a formal "package" 
  #   because too much file making, naming, etc.. but has made our lives easiers... thoughts?

chromPrep <- function(CHR, genome, a){
  v           <- seq_along(genome[[CHR]][1,])
  fname       <- sprintf("prep_%s.txt",CHR)
  chr.file    <- this.lapply(v, function(L){getCline(l = genome[[CHR]][,L] , a=a)} )
  chr.file    <- do.call(rbind, chr.file)
  write.table(chr.file, file = fname, quote=F,row.names=F,col.names=F )
}

allChromPreps = function(genome,a){ 
  lapply(names(genome), chromPrep, genome = genome, a = a)  
}


loadProcessChr <- function(chr.filename, chr.start, chr.end){
  skip.chr   <- (chr.start - 1)
  length.chr <- countLines(chr.filename)[[1]]
  if(length(chr.end) != 0){chr.end <- min(chr.end,length.chr)}
  nlines.chr <- ifelse(length(chr.end)==0, 0, chr.end - skip.chr)
  nrow.chr   <- ifelse(nlines.chr == 0, length.chr - skip.chr, nlines.chr)
  chr.dat    <- apply(matrix(
      scan(chr.filename, nlines = nlines.chr, skip = skip.chr), 
      nrow = nrow.chr, byrow = TRUE),1,processLocus)
  return(chr.dat)
}

betweenChrCompSplit <- function(CHR_A.FILENAME, CHR_B.FILENAME, a, s1 = 1, e1 = numeric(0), s2 = 1, e2 = numeric(0), return.table = FALSE){
  #note these should be the names of the chrs; ie prep prep_'CHR_A'.txt
  #added s for start (must be >=1)
  #added e for end (must be <= number of markers in that prep file)
  a.dat    <- loadProcessChr(CHR_A.FILENAME, s1, e1)
  b.dat    <- loadProcessChr(CHR_B.FILENAME, s2, e2)
  along.a  <- seq_along(a.dat)
  along.b  <- seq_along(b.dat)
  all.comps = do.call(rbind,lapply( along.a ,function(LA){
    do.call(rbind,lapply(along.b ,function(LB){
      c(la=LA,lb=LB,LDcalcs(l1 = a.dat[[LA]], l2 = b.dat[[LB]], a = a))      
    }))
  }))
  colnames(all.comps)[1] = unlist(strsplit("CHR_A.FILENAME",".",fixed=T))[1]
  colnames(all.comps)[2] = unlist(strsplit("CHR_B.FILENAME",".",fixed=T))[1]
  resultsfilename=paste("results",colnames(all.comps)[1],"start",s1,"end",e1,colnames(all.comps)[2],"start",s2,"end",e2,"comp.txt",sep="_")
  #write.table(all.comps,file = resultsfilename, quote=F,row.names=F )
  if(return.table){return(all.comps)}
}


# Example usage
setwd("/Users/ybrandva/Dropbox/LD_lowess_simulations_with_Yaniv/forYaniv_troubleshoot_mem") 
arrArgs <- unlist(strsplit("genotypes_CALM_updatemsg.txt hybrid_index_CALM_msgupdate 1 2", " "))
hybridindex<-unlist(read.table(arrArgs[2],header=TRUE,as.is=T))
focalchrom1<-as.numeric(arrArgs[3])
focalchrom2<-as.numeric(arrArgs[4])
betweenChrCompSplit("prep_1.txt","prep_2.txt",hybridindex, s1 = 1, e1 =2, return.table = TRUE)

