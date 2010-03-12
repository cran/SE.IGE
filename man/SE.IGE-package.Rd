\name{SE.IGE-package}
\alias{SE.IGE-package}
\alias{SE.IGE}
\docType{package}
\title{Package for standard errors of estimated genetic parameters}
\description{Package for standard errors of estimated genetic parameters
  predicts standard (SE) errors of estimated direct, indirect and total 
  genetic variance based on one-way ANOVA's. Also for optimum family size}
\details{
  \tabular{ll}{
  Package: \tab SE.IGE\cr
  Type: \tab Package\cr
  Version: \tab 1.0\cr
  Date: \tab 2010-03-12\cr
  License: \tab What license is it under?\cr
  LazyLoad: \tab yes\cr
 }
 Most important functions:
 SE_random: SE for schemes with groups composed at random
 SE_2families: SE for schemes with groups composed at random
 nopt_random: optimum family sizes for schemes with groups composed at random
 nopt_2families:  optimum family sizes for schemes with groups composed of 2 families
}
\author{Piter Bijma; piter.bijma@wur.nl}
\references{Bijma, 2010. Estimating social genetic effects: precision of 
   estimates and optimum designs. Submitted to Genetics March 2010}
\keyword{ package }
\examples{
  N <- 100
  n <- 20
  vpd <- 100
  vps <- 10
  h2d <- 0.3
  h2s <- 0.3
  rgds <- 0.1
  reds <- 0.1
  nw <- 4
  r <- 0.5
  se_r <- SE_random(N,n,vpd,vps,h2d,h2s,rgds,reds,nw,r)
  se_2f <- SE_2families(N,n,vpd,vps,h2d,h2s,rgds,reds,nw,r)
  T <- N*n
  optima_r <- nopt_random(T,vpd,vps,h2d,h2s,rgds,reds,nw,r)
  optima_2f <- nopt_2families(T,vpd,vps,h2d,h2s,rgds,reds,nw,r)
}
