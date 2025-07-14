## computes covariance of S and M, generate ellipses for plotting
library(ellipse)

## assemble covariance from pp where

## upper / lower blocks of Cov
P3 <- function(pp){
  Sig <- diag(pp*(1-pp))
  Sig[1,2] <- -pp[1]*pp[2]
  Sig[1,3] <- -pp[1]*pp[3]
  Sig[2,3] <- -pp[2]*pp[3]
  Sig[lower.tri(Sig)] <- Sig[upper.tri(Sig)]  
  return(Sig)
}

## entire Cov
ConstructCov <- function(pp,nom){
  Zero <- matrix(0,nrow=3,ncol=3)
  return(rbind(cbind(P3(pp[1:3]),Zero),nom*cbind(Zero,P3(pp[4:6]))))
}

## point estimates of stage shift and mortality reduction
SMest <- function(ppC,ppS){
  Sest <- (sum(ppC[,2]) - sum(ppS[,2]))/sum(ppC[,2])
  Mest <- (sum(ppC[2,]) - sum(ppS[2,]))/sum(ppC[2,])
  return(c(Sest,Mest))
}


## Covariance of SM using analytic gradient
## pp = (control prob death AND early/no cancer, control prob alive AND late, control prob dead AND late,
##        screen prob death AND early/no cancer,  screen prob alive AND late,  screen prob dead AND late)
##  n = control sample size
##  m = screen sample size
ConstructCovSM <- function(pp,n,m){
  Cov <- ConstructCov(pp,nom=n/m)
  rc <- pp[2] + pp[3]
  rs <- pp[5] + pp[6]
  sc <- pp[1] + pp[3]
  ss <- pp[4] + pp[6]
  grad1 <- (rs/rc)*c(0,1/rc,1/rc,0,-1/rs,-1/rs)
  grad2 <- (ss/sc)*c(1/sc,0,1/sc,-1/ss,0,-1/ss)
  Grad <- cbind(grad1,grad2)
  gCov <- (t(Grad)%*%Cov%*%Grad)/n
  return(gCov)
}


## Covariance of SM using marginal distribution
##.  useful for sensitivity analyses
##
## pp = (rc,rs,sc,ss,rho)
##
##  rc = control prob late
##  rs = screen prob late
##. sc = control prob death
##. ss = screen prob death
##. rho = correlation between S and M
##. n = control size
##  m = screen size
ConstructCovSMMargin <- function(pp,n,m){
  sm11 <- (pp[2]^2/pp[1]^3)*(1-pp[1]) + (n/m)*(1/pp[1]^2)*pp[2]*(1-pp[2])
  sm22 <- (pp[4]^2/pp[3]^3)*(1-pp[3]) + (n/m)*(1/pp[3]^2)*pp[4]*(1-pp[4])
  sm12 <- pp[5]*sqrt(sm11)*sqrt(sm22)
  gCov <- matrix(c(sm11,sm12,sm12,sm22),nrow=2)/n
  return(gCov)
}

## returns a set of x,y points denoting elliptical confidence region for S,M
## arguments
##        pp = parameters or estimates, see ConstructCovSM, ConstructCovSMMargin for format
##      cent = center of ellipse
##         n = control sample size
##         m = screen sample size
##  marginal = should Confidence region be constructed only using marginal probs of death and late and corr
CreateEllipse <- function(pp,cent,n,m,alpha=0.05,marginal=FALSE){
  if(!marginal){
    gCov <- ConstructCovSM(pp,n,m)
  } else {
    gCov <- ConstructCovSMMargin(pp,n,m)
  }
  el <- ellipse(gCov,centre=cent,level=1-alpha)
  el <- as.data.frame(el)
  colnames(el) <- c("x","y")
  return(el)
}


## format joint probability table for printing
FormatTable <- function(tab){
  colnames(tab) <- c("No Cancer","Early","Late")
  rownames(tab) <- c("No Death","Death")
  return(tab)
}








## Code modified from FiellerRatio in package twopartm_0.1.0
## the original code can be accessed by installing package and typing
## getMethod("FiellerRatio","numeric")
## returns a list with elements
##
##   res = results (i.e. interval and point estimate)
## latex = latex formatted results
##  type = either (0= [lb,ub]), (1=[-inf,bound1] U [bound2,inf]), or (2=[-inf,inf])
FiellerCI <- function (xest, yest, V, alpha = 0.05,round_digits=1) {
  if (!is.numeric(xest) | !is.numeric(yest) | length(xest) != 
      1 | length(yest) != 1) {
    stop("xest or yest should be the estimates of two variables.")
  }
  if (is.vector(V) | is.null(V)) {
    stop("Covariance matrix of two estimates should be given.")
  }
  Q <- xest/yest
  varx <- V[1, 1]
  vary <- V[2, 2]
  covxy <- V[1, 2]
  z <- qnorm(1 - alpha/2)
  if (yest^2/vary > z^2) {
    m = xest * yest - z^2 * covxy
    minB <- (m - sqrt(m^2 - (yest^2 - z^2 * vary) * (xest^2 - 
                                                       z^2 * varx)))/(yest^2 - z^2 * vary)
    maxB <- (m + sqrt(m^2 - (yest^2 - z^2 * vary) * (xest^2 - 
                                                       z^2 * varx)))/(yest^2 - z^2 * vary)
    return(list(res=c(ratio = Q, min = minB, max = maxB),
                latex=paste0("$(",round(minB,round_digits),",",round(maxB,round_digits),")$"),
                type=0))
  }
  else {
    z2_unb = yest^2/vary + (xest * vary - yest * covxy)^2/(vary *
                                                             (varx * vary - covxy^2))
    if (z2_unb > z^2) {
      m = xest * yest - z^2 * covxy
      B1 <- (m - sqrt(m^2 - (yest^2 - z^2 * vary) * (xest^2 -
                                                       z^2 * varx)))/(yest^2 - z^2 * vary)
      B2 <- (m + sqrt(m^2 - (yest^2 - z^2 * vary) * (xest^2 -
                                                       z^2 * varx)))/(yest^2 - z^2 * vary)
      minB <- min(B1, B2)
      maxB <- max(B1, B2)
      mess <- "Notice: The denominator is not significantly different from zero, and the confidence set of the ratio is combining [-inf,bound1] and [bound2,inf]."
      return(list(res=c(ratio = Q, bound1 = minB, bound2 = maxB),
                  latex=paste0("($-\\infty,",round(minB,round_digits),") \\cup (",round(maxB,round_digits),",\\infty$)"),
                  message=mess,
                  type=1))
    }
    else {
      return(list(res=c(ratio=Q,min=-Inf,max=Inf),latex="$(-\\infty,\\infty)$",type=2))
    }
  }
}
