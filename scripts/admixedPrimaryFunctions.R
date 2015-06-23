# Yaniv Brandvain & Molly Schumer 6/4/15
# Primary functions for admixe'D
# Used to calculate LD in admixed population [controlling for heterogeneity in ancestry proportions]
# If we where to release a packge these would probably be the primary functions.

#install.packages("ppcor")
#library(ppcor)
pcor <- function (x, method = c("pearson", "kendall", "spearman")) 
{
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

pcor.test <- function (x, y, z, method = c("pearson", "kendall", "spearman")) 
{
  method <- match.arg(method)
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  xyz <- data.frame(x, y, z)
  pcor = pcor(xyz, method = method)
  data.frame(estimate = pcor$est[1, 2], p.value = pcor$p.value[1, 2], 
             statistic = pcor$statistic[1, 2], n = pcor$n, gp = pcor$gp, 
             Method = method)
}

# Here we make class "admx" (probably not necessary), but why not?
setClass(Class="admx", 
         slots = c( 
           n.inds       = "numeric",
           n.reps       = "numeric",
           l            = "numeric",
           sim.l        = "matrix",
           cline        = "numeric",
           l.rnk.lklhd  = "numeric",
           sim.cline    = "matrix"
         )
)

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

processLocus <- function(l.vals){
  # This takes in, as data, a long vector that is the output of getCline (below)
  # Returns a well formated verstion of this data as an object of class "admx"  
  n.inds       <-   l.vals[1]
  n.reps       <-   l.vals[2]
  l            <-   l.vals[3:(2+n.inds)]
  last         <-   2+n.inds
  sim.l        <-   matrix(l.vals[(1+last):(last+n.reps*n.inds)], nrow = n.inds)
  last         <-   last+n.reps*n.inds
  cline        <-   l.vals[(1+last):(last+n.inds)]
  last         <-   last+n.inds
  sim.cline    <-   matrix(l.vals[(1+last):(last+n.reps*n.inds)], nrow = n.inds)
  l.rnk.lklhd  <-   l.vals[length(l.vals)] 
  return(new(Class="admx", 
             n.inds = n.inds, n.reps = n.reps, l = l, sim.l = sim.l, cline = cline, l.rnk.lklhd = l.rnk.lklhd, sim.cline = sim.cline
  ))
}

getCline <- function(l, a, reps = 1000, return.processed = FALSE){
  # summarizes the "genomic cline" for a locus, includes raw data, and data generated for the null distribution
  old.l        <- l
  na.inds      <- which(is.na(l))
  ok.inds      <- !seq_along(l)%in%na.inds
  l            <- l[ok.inds]
  cline        <- try(gcline(a[ok.inds], n = 2, 2*l))
  tryerror     <- class(cline) == "try-error"
  if(class(cline) == "try-error" ){
    chr.admx       <- c(length(old.l), reps,old.l,rep(NA,reps*length(old.l)),
                    rep(NA,length(old.l)), rep(NA,reps*length(old.l)), NA )
    if(return.processed){return(  processLocus(chr.admx)  )}
    return(chr.admx)
  }
  sim.l        <- replicate(reps,rbinom(length(cline),2,cline)/2)
  sim.vals     <- apply(sim.l,2,function(p){
    sim.cline      <- try(gcline(a[ok.inds], n = 2, 2*p))
    if(class(sim.cline)=="try-error") {sim.cline <- rep(NA,sum(ok.inds))}
    return(sim.cline)
  })
  p.l.given.cline <- sum(dbinom(2*l,2,cline,log=T))
  p.l.sims.given.cline <- apply(sim.l,2,function(X){ sum(dbinom(2*X,2,cline,log=T)) })
  l.rnk.lklhd  <- sum(p.l.given.cline  > p.l.sims.given.cline)/ reps
  #dealing with NAs
  tmp.sim.l    <- matrix(nrow= length(a),ncol=reps); tmp.sim.l[ok.inds,] =  sim.l; sim.l = tmp.sim.l
  tmp.cline    <- rep(NA,length(old.l)); tmp.cline[ok.inds] = cline; cline = tmp.cline
  tmp.sim.vals <- matrix(nrow= length(a),ncol=reps); tmp.sim.vals[ok.inds,] = sim.vals; sim.vals = tmp.sim.vals
  chr.admx     <- c( length(old.l), reps, old.l, c(sim.l), cline, c(sim.vals), l.rnk.lklhd )
  if(return.processed){return(  processLocus(chr.admx)  )}
  return(chr.admx)
}

LDcalcs <- function(l1,l2,a){
  ok <- !is.na(l1@l + l2@l)
  these.genos <- cbind(A = l1@l, B = l2@l)[ok,]
  allele.counts <- c(colSums(these.genos), colSums(abs(1-these.genos)))
  names(allele.counts) <- c("A","B","a","b")
  geno.table <- table(data.frame(A=factor(these.genos[,1], levels = c(1,.5,0)),B=factor(these.genos[,2], levels = c(1,.5,0))))
  geno.counts <- c(  
    AB = sum(geno.table[1,] * c(1,1/2,0)) + sum(geno.table[2,]/2 * c(1,1/2,0)),
    Ab = sum(geno.table[1,] * c(0,1/2,1)) + sum(geno.table[2,]/2 * c(0,1/2,1)),
    aB = sum(geno.table[3,] * c(1,1/2,0)) + sum(geno.table[2,]/2 * c(1,1/2,0)),
    ab = sum(geno.table[3,] * c(0,1/2,1)) + sum(geno.table[2,]/2 * c(0,1/2,1))
  )
  cln.D        <- mean(  (l1@l - l1@cline) * (l2@l - l2@cline), na.rm=T)
  cln.R        <- cln.D / (sqrt(mean((l1@l-l1@cline)[ok]^2))*sqrt(mean((l2@l-l2@cline)[ok]^2)))
  # sim.cln.D    <- colMeans(l1@sim.l*l2@sim.l - l1@sim.cline*l2@sim.cline,na.rm=T)
  sim.cln.D    <- colMeans(  (l1@sim.l- l1@sim.cline) * (l2@sim.l- l2@sim.cline),na.rm=T)
  sim.cln.R    <- sim.cln.D / 
    (sqrt(colMeans((l1@sim.l-l1@sim.cline)[ok,]^2)) * sqrt(colMeans((l2@sim.l-l2@sim.cline)[ok,]^2)))
  no.var <- length(unique(l1@l[ok])) == 0 |  length(unique(l2@l[ok]) )== 0
  if(no.var){pcor = c(estimate= NA, p.value = NA, statistic = NA)}
  if(!no.var){pcor = pcor.test( these.genos[,1], these.genos[,2],a[ok])[1:3]}
  results <-c(
    allele.counts, 
    geno.counts,
    reg.D         = mean(l1@l*l2@l,na.rm=T) - mean(l1@l,na.rm=T) * mean(l2@l,na.rm=T),
    adx.D         = mean(l1@l*l2@l-a^2,na.rm=T),
    cln.D         = cln.D,
    p.cln.D       = 2 * min(sum(cln.D > sim.cln.D,na.rm=T) , sum(cln.D < sim.cln.D,na.rm=T)) / sum(!is.na(sim.cln.D)),
    cln.R         = cln.R,
    p.cln.R       = 2 * min(sum(cln.R > sim.cln.R,na.rm=T) , sum(cln.R < sim.cln.R,na.rm=T)) / sum(!is.na(sim.cln.D)),
    n.chances     = sum(!is.na(sim.cln.D)),
    pcor
  )
  return( unlist(results) )
}
