
pcor <- function (x, method = c("pearson", "kendall", "spearman")) {
  method <- match.arg(method)
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x)) 
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  n <- dim(x)[1]
  gp <- dim(x)[2] - 2
  cvx <- cov(x, method = method)
  if(det(cvx) == 0){
    return(list(estimate = matrix(rep(NA,9),nrow=3), p.value = matrix(rep(NA,9),nrow=3), statistic = matrix(rep(NA,9),nrow=3), 
                n = n, gp = gp, method = method)) }
  icvx <- solve(cvx)
  pcor <- -cov2cor(icvx)
  diag(pcor) <- 1
  if (method == "kendall") {
    statistic <- pcor/sqrt(2 * (2 * (n - gp) + 5)/(9 * (n -  gp) * (n - 1 - gp)))
    p.value <- 2 * pnorm(-abs(statistic))
  }
  else {
    statistic <- pcor * sqrt((n - 2 - gp)/(1 - pcor^2))
    p.value <- 2 * pnorm(-abs(statistic))
  }
  diag(statistic) <- 0
  diag(p.value) <- 0
  list(estimate = pcor, p.value = p.value, statistic = statistic, 
       n = n, gp = gp, method = method)
}

pcor.test <- function (x, y, z, method = c("pearson", "kendall", "spearman")) {
  method <- match.arg(method)
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  xyz <- data.frame(x, y, z)
  if( sum(abs(with(xyz,x-y))) == 0){
    return(data.frame(estimate = NA, p.value = NA, 
               statistic = NA, n = NA, gp = NA, 
               Method = method)  )
  }
  pcor = pcor(xyz, method = method) 
  data.frame(estimate = pcor$est[1, 2], p.value = pcor$p.value[1, 2], 
             statistic = pcor$statistic[1, 2], n = pcor$n, gp = pcor$gp, 
             Method = method)
}

LDcalcs <- function(l1,l2,a){
  ok <- !is.na(l1 + l2)
  these.genos <- cbind(A = l1, B = l2)[ok,]
#  allele.counts <- c(colSums(these.genos), colSums(abs(1-these.genos)))
#  names(allele.counts) <- c("A","B","a","b")
#  geno.table <- table(data.frame(A=factor(these.genos[,1], levels = c(1,.5,0)),B=factor(these.genos[,2], levels = c(1,.5,0))))
#  geno.counts <- c(  
#    AB = sum(geno.table[1,] * c(1,1/2,0)) + sum(geno.table[2,]/2 * c(1,1/2,0)),
#    Ab = sum(geno.table[1,] * c(0,1/2,1)) + sum(geno.table[2,]/2 * c(0,1/2,1)),
#    aB = sum(geno.table[3,] * c(1,1/2,0)) + sum(geno.table[2,]/2 * c(1,1/2,0)),
#    ab = sum(geno.table[3,] * c(0,1/2,1)) + sum(geno.table[2,]/2 * c(0,1/2,1))
#  )
  p.cor = pcor.test( these.genos[,1], these.genos[,2],a[ok])[1:3]
  results <-c(
    #allele.counts, 
    #geno.counts,
    reg.D         = mean(l1*l2,na.rm=T) - mean(l1[ok]) * mean(l2[ok]),
    reg.R         = (mean(l1*l2,na.rm=T) - mean(l1[ok]) * mean(l2[ok])) / sqrt(var(l1[ok])*var(l2[ok])),
    adx.D         = mean(l1*l2-a^2,na.rm=T),
    p.cor
  )
  return( unlist(results) )
}






















LL.fn <- function(par, SD = SD, x, n, y) {
  # calculate the log-likelihood of the data given the cline
  u            <- par[1]
  v            <- par[2]
  p            <- x^v/(x^v + (1 - x)^v * exp(u))
  p            <- replace(p, x == 1, 1)
  p            <- replace(p, x == 0, 0)
  p            <- replace(p, p <= 0, .Machine$double.xmin)
  p            <- replace(p, p >= 1, 1 - .Machine$double.neg.eps)
  this.LL      <- sum(dbinom(y, n, p, log = TRUE))
}

gcline <- function(x, n, y){
  # Modified from the HIest package (http://cran.r-project.org/web/packages/HIest/index.html)
  # Generates a maximum likelihood 'genomic cline'
  # n = ploidy
  start        <- c(u = 0, v = 1)
  k            <- 2
  SD           <- rep(0.1, 2)
  lower        <- c(-Inf, 0)
  upper        <- c(Inf, Inf)
  null.LL      <- LL.fn(par = c(u = 0, v = 1), SD, x, n, y)
  est          <- optim(par = start, fn = LL.fn, x = x, n = n, y = y, method = "L-BFGS-B", lower = lower, upper = upper, control = list(fnscale = -1))
  ests         <- est$par
  convergence  <- est$convergence
  names(ests)  <- names(start)
  u            <- ests["u"]
  v            <- ests["v"]
  this.gc      <- (x^v/(x^v+(1-x)^v *exp(u)))
  return(this.gc)
}


getCline <- function(l, a, return.processed = FALSE){
  # summarizes the "genomic cline" for a locus, includes raw data, and data generated for the null distribution
  old.l        <- l
  na.inds      <- which(is.na(l))
  ok.inds      <- !seq_along(l)%in%na.inds
  l            <- l[ok.inds]
  cline        <- try(gcline(a[ok.inds], n = 2, 2*l))
  tryerror     <- class(cline) == "try-error"
  if(class(cline) == "try-error" ){
    return(lklhd = NA)
  }
  return( sum(dbinom(2*l,2,cline,log=T)) / length(l))
}



#trimAncestryProp (geno.data)

trimAncestryProp <- function(all.loci, trim = .20, method = "lowcor" ){
  #all.loci genotypes for everybody [rows = inds, columns = loci]
  # "lowcor",  "ancestry", "cline", "hit"
  doTrim <- function(this.stat){
    normals <- this.stat > quantile( this.stat, prob = trim)
    list(new.alpha = rowMeans(all.loci[,normals]), weirdos = which(!normals) )
  }
  initial.alpha <- rowMeans(all.loci,na.rm=T)
  allele.freqs <- colMeans(all.loci,na.rm=T) 
  if(method == "lowcor"){
    # Remove loci with the lowest correlation with alpha
    this.tmp <- c( cor(initial.alpha,  all.loci, use = "pairwise.complete.obs") )
    return( doTrim(  this.tmp  ) )
    rm(this.tmp)
  }
  if(method == "ancestry"){
    # ameasure of how poorly the model fits given ancestry propotion, 
    #standadized by ample size and expected variance.
    this.tmp <- colSums(apply(
      2*all.loci,2,  dbinom,  size = 2, prob = initial.alpha , log = T),
      na.rm = T)
    this.tmp <- this.tmp  / ( colSums(!is.na(all.loci))) 
    return( doTrim(  this.tmp  ) )
    rm(this.tmp)
  }
  if(method == "cline"){
    this.tmp <- apply(all.loci,2,function(X){as.numeric(getCline(l = X, a = initial.alpha))  })
    return( doTrim(  this.tmp  ) )
    rm(this.tmp)
  }
  if(method == "hit"){
    this.tmp <-colSums(all.loci == .5,na.rm=T) / 
      colSums(2 * initial.alpha  * (1-initial.alpha ) * !is.na(all.loci))
    return( doTrim(  this.tmp  ) )
    rm(this.tmp)
  }
}


