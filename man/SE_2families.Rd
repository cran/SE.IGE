\name{SE_2families}
\alias{SE_2families}
\title{standard error for estimated genetic parameters for 2-family groups}
\description{predicts the standard error of estimated genetic parameters
  for direct and indirect effects when groups of interacting individuals
  are composed of two families, each contributing half the individuals}
\usage{SE_2families(N,n,vpd,vps,h2d,h2s,rgds,reds,nw,r)}
\arguments{
  \item{N}{number of families, numeric}
  \item{n}{family size, numeric}
  \item{vpd}{phenotypic (i.e. full) variance of direct effect, numeric}
  \item{vps}{phenotypic variance of indirect ("associative) effect, numeric}
  \item{h2d}{heritability of direct effect, Var(DGE) = h2d*vpd, numeric}
  \item{h2s}{heritability of indirect effect, Var(IGE) = h2s*vps, numeric}
  \item{rgds}{genetic correlation direct-indirect effect}
  \item{reds}{residual correlation direct-indirect effect}
  \item{nw}{group size, i.e. total number of individuals in a group, numeric}
  \item{r}{additive genetic relatedness between family members, 
           e.g. r=0.25 for half sibs, numeric}
}
\value{
   \item{SEvad}{standard error of estimated direct genetic variance}
   \item{SEvas}{standard error of estimated indirect genetic variance}
   \item{SEvtbv}{standard error of estimated total genetic variance}
   \item{SET2}{standard error of estimated ratio total genetic 
                over phenotypic variance}
   \item{SEcads}{standard error of estimated direct-indirect 
                  genetic covariance}
   \item{SErgds}{standard error of estimated direct-indirect
                  genetic correlation}
}
\references{Bijma, 2010. Estimating social genetic effects: precision of 
   estimates and optimum designs. Submitted to Genetics March 2010}
\author{P.Bijma; piter.bijma@wur.nl}
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
  se_2f <- SE_2families(N,n,vpd,vps,h2d,h2s,rgds,reds,nw,r)
}
\keyword{utilities}

