##
#
# Notice: This code is here such that pyrace can test against it.
#   The original code can be found at
#        http://dl.dropboxusercontent.com/u/15624220/Loganetal.zip
#
# Please note the copyright-notice below.
#

################################################################################
############### UNIFORM CRITERION VARIABILITY WALD  ACCUMULATOR ################
################################################################################
# Comments and changes added by Andrew Heathcote

unit.pars <- c("ter","A","b","v")

pdf1=function(dt,pa,i) {
  # pdf for a single unit
  if ( length(dt)>dim(pa)[3] )
    stop("dt longer than parameters in pdf1") else
      pa <- pa[,,1:length(dt),drop=FALSE]

  # reparameterize boundrary at k +/- U(0,a) to
  # boundary b = k + A/2 and start point noise U(0,A)
  # hence k = b-A/2 and a=A/2 and v=l (i.e., mean rate)
  digt(dt,k=pa["b",i,]-pa["A",i,]/2,
       l=pa["v",i,],a=pa["A",i,]/2)
}

cdf1=function(dt,pa,i) {
  # cdf for a single unit

  if ( length(dt)>dim(pa)[3] )
    stop("dt longer than parameters in cdf1") else
      pa <- pa[,,1:length(dt),drop=FALSE]
  
  # reparameterize boundrary at k +/- U(0,a) to
  # boundary b = k + A/2 and start point noise U(0,A)
  # hence k = b-A/2 and a=A/2 and v=l (i.e., mean rate)
  pigt(dt,k=pa["b",i,]-pa["A",i,]/2,
         l=pa["v",i,],a=pa["A",i,]/2)
}

rfun1=function(pa) {
  # returns dim(pa)[2:3] matrix

  nsamp <- prod(dim(pa)[2:3])
  matrix(rwald(n=nsamp,l=pa["v",,],
      k=pa["b",,] - runif(nsamp,0,pa["A",,])) +
    pa["ter",,],nrow=dim(pa)[2]) 
}


########################## pdf1/cdf1 SUPPORT FUNCTONS #########################
# Copyright (C) 2013  Trisha Van Zandt distributed with: Psychological Review.
# Logan, Van Zandt, Verbruggen, and Wagenmakers (in press).  On the ability to 
# inhibit thought and action: General and special theories of an act of control.
# Comments and changes added by Andrew Heathcote

rwald <- function(n,k,l) {
  # random sample of n from a Wald (or Inverse Gaussian)
  # k = criterion, l = rate, assumes sigma=1 Browninan motion
  
  rlevy <- function(n=1, m=0, c=1) {
    if (any(c<0)) stop("c must be positive")
    c/qnorm(1-runif(n)/2)^2+m
  }
  
  if(length(k)==1) k <- rep(k,times=n)
  if(length(l)==1) l <- rep(l,times=n)
  
  tiny <- 1e-6
  flag <- l>tiny
  x <- rep(NA,times=n)
  
  x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
  mu <- k/l
  lambda <- k^2
  
  y <- rnorm(sum(flag))^2
  mu.0 <- mu[flag]
  lambda.0 <- lambda[flag]
  
  x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
    sqrt(4*mu.0*lambda.0*y +
           mu.0^2*y^2)*mu.0/(2*lambda.0)
  
  z <- runif(length(x.0))
  test <- mu.0/(mu.0+x.0)
  x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
  x[flag] <- x.0
  x[x<0] <- max(x)
  x
}

digt.0 <- function(t,k=1,l=1) {
  # pdf of inverse gaussian at t with no k variability
  
  options(warn=-1)
  if(length(k)==1) k <- rep(k,times=length(t))
  if(length(l)==1) l <- rep(l,times=length(t))
  
  mu <- k/l
  lambda <- k^2
  
  e <- -(lambda/(2*t)) * (t^2/mu^2 - 2*t/mu  + 1)
  
  e[mu==Inf] <- -.5*lambda/t
  con <- .5*log(lambda) - .5*log(2*t^3*pi)
  
  x <- exp(e+con)
  x[t<=0] <- 0
  x
}

pigt.0 <- function(t,k=1,l=1) {
  # cdf of inverse gaussian at t with no k variability
  
  options(warn=-1)
  if(length(k)==1) k <- rep(k,times=length(t))
  if(length(l)==1) l <- rep(l,times=length(t))
  
  mu <- k/l
  lambda <- k^2
  
  e <- exp(log(2*lambda) - log(mu))
  add <- sqrt(lambda/t) * (1 + t/mu)
  sub <- sqrt(lambda/t) * (1 - t/mu)
  
  p.1 <- 1 - pnorm(add)
  p.2 <- 1 - pnorm(sub)
  
  x <- exp(e + log(p.1)) + p.2
  
  x[t<0] <- 0
  x
  
}

digt <- function(t,k=1,l=1,a=1) {
  # pdf of inverse gaussian at t with k +/-a/2 uniform variability
  # returns digt.0 if a<1e-10
  
  options(warn=-1)
  if(length(k)==1) k <- rep(k,times=length(t))
  if(length(l)==1) l <- rep(l,times=length(t))
  if(length(a)==1) a <- rep(a,times=length(t))
  tiny <- 10^(-10)
  
  a[a<=tiny] <- 0
  l[l<=tiny] <- 0
  
  sqr.t <- sqrt(t)
  log.t <- log(t)

  ## formula 16
  term.1a <- -(a-k+t*l)^2/(2*t)
  term.1b <- -(a+k-t*l)^2/(2*t)
  term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t)
  
  term.2a <- log(.5)+log(l)
  term.2b <- 2*pnorm((-k+a)/sqr.t+sqr.t*l)-1
  term.2c <- 2*pnorm((k+a)/sqr.t-sqr.t*l)-1
  term.2d <- term.2b+term.2c
  term.2 <- exp(term.2a)*term.2d
  
  term.3 <- term.1+term.2
  term.4 <- log(term.3)-log(2)-log(a)
  x <- exp(term.4)
  
  ## formula 17
  term.1 <- -.5*(log(2)+log(pi)+log.t)
  term.2 <- (k-a)^2/(2*t)
  term.3 <- (k+a)^2/(2*t)
  term.4 <- (exp(-term.2[l==0])-exp(-term.3[l==0]))
  term.5 <- term.1[l==0]+log(term.4) - log(2) - log(a[l==0])
  x[l==0] <- exp(term.5)

  ## for a<1e-10
  x[a<=tiny] <- digt.0(t[a<=tiny],k[a<=tiny],l[a<=tiny])
  
  x[t<=0] <- 0
  x[x<0 | is.nan(x) ] <- 0
  
  x
}

pigt <- function(t,k=1,l=1,a=1) {
  # cdf of inverse gaussian at t with k +/-a/2 uniform variability
  # returns pigt.0 if a<=0
  
  options(warn=-1)
  if(length(k)==1) k <- rep(k,times=length(t))
  if(length(l)==1) l <- rep(l,times=length(t))
  if(length(a)==1) a <- rep(a,times=length(t))
  tiny <- 10^(-10)
  
  l[l<=tiny] <- 0
  a[a<=tiny] <- 0
  
  sqr.t <- sqrt(t)
  log.t <- log(t)


  ## normal
  term.1a <- .5*log.t-.5*log(2*pi)
  term.1b <- exp(-((k-a-t*l)^2/t )/2)
  term.1c <- exp(-((k+a-t*l)^2/t )/2)
  term.1 <- exp(term.1a)*(term.1b-term.1c)
  
  term.2a <- exp(2*l*(k-a)+ log(pnorm(-(k-a+t*l)/sqr.t)))
  term.2b <- exp(2*l*(k+a)+ log(pnorm(-(k+a+t*l)/sqr.t)))
  term.2 <- a + (term.2b-term.2a)/(2*l)
  
  term.4a <- 2*pnorm((k+a)/sqr.t-sqr.t*l)-1
  term.4b <- 2*pnorm((k-a)/sqr.t-sqr.t*l)-1
  term.4c <- .5*(t*l - a - k + .5/l)
  term.4d <- .5*(k - a - t*l - .5/l)
  term.4 <- term.4c*term.4a + term.4d*term.4b
  
  term <- (term.4 + term.2 + term.1)/(2*a)
  x <- (term)


  ## l==0 case
  term.5a <- 2*pnorm( (k+a)/sqr.t)-1
  term.5b <- 2*pnorm(-(k-a)/sqr.t)-1
  term.5 <- (-(k+a)*term.5a - (k-a)*term.5b)/(2*a)
  
  term.6a <- -.5*(k+a)^2/t - .5*log(2) -.5*log(pi) + .5*log.t - log(a)
  term.6b <- -.5*(k-a)^2/t - .5*log(2) -.5*log(pi) + .5*log.t - log(a)
  term.6 <- 1 + exp(term.6b) - exp(term.6a)
  
  term.7 <- term.5 + term.6
  
  x[l<=0] <- term.7[l<=0]


  ## a==0 case
  x[a<=0] <- pigt.0(t[a<=0],k[a<=0],l[a<=0])


  ## error-catching
  x[t<0] <- 0
  x[x<0 | is.nan(x) ] <- 0
  x
  
}

################################# Checks #######################################

# n=1000 # points on graphs, /100 to give seconds
# 
# # Note: ter=0 as not used by pdf1 and cdf1
# pars1 <- c(ter=0,A=1,b=2,v=1)   # "correct" error accumulator
# pars2 <- c(ter=0,A=1,b=2,v=.5)  # "error" accumulator
# 
# # For drawing graphs
# pa <- array(0,dim=c(length(unit.pars),2,n),
#             dimnames=list(unit.pars,NULL,NULL))
# pa[,1,] <- pars1; pa[,2,] <- pars2
# 
# # Check density area (shoudl be 1)
# integrate(pdf1,0,Inf,pa=pa,i=1) 
# integrate(pdf1,0,Inf,pa=pa,i=2) 
# 
# # Large random sample
# pa1 <- array(0,dim=c(length(unit.pars),2,n*100),
#              dimnames=list(unit.pars,NULL,NULL))
# pa1[,1,] <- pars1; pa1[,2,] <- pars2
# samp1 <- rfun1(pa1)
# 
# # Plot pdf1 and cdf1 and overlapy random sample on former
# x=c(1:n)/100; probs <- (1:99)/100; par(mfcol=c(1,2))
# 
# # 1st Accumulator
# plot(x,pdf1(dt=x,pa,1),type="l");abline(h=0)
# lines(density(samp1[1,samp1[1,]<max(x)]),col="red")
# plot(x,cdf1(dt=x,pa,1),type="l",ylim=c(0,1))
# lines(quantile(samp1[1,],probs=probs),probs,col="red")
# abline(h=0); abline(h=1)
# 
# # 2nd Accumulator
# plot(x,pdf1(dt=x,pa,2),type="l");abline(h=0)
# lines(density(samp1[2,samp1[2,]<max(x)]),col="red")
# plot(x,cdf1(dt=x,pa,2),type="l",ylim=c(0,1))
# lines(quantile(samp1[2,],probs=probs),probs,col="red")
# abline(h=0); abline(h=1)
