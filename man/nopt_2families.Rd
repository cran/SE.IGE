\name{nopt_2families}
\alias{nopt_2families}
\title{function for optimum family sizes and standard errors of
       genetic parameters}
\description{function calculating optimum family sizes and corresponding 
             standard errors of genetic parameters for direct, indirect,
             total genetic variance, and for T2. For groups composed of two
             families}
\usage{nopt_2families(T,vpd,vps,h2d,h2s,rgds,reds,nw,r)}
\arguments{
  \item{T}{total number of individuals T, #families*family size, numeric}
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
\value{returns a matrix, 4 rows x 3 columns
      row=1: direct genetic effects
      row=2: indirect genetic effects
      row=3: total genetic effects
      row=4: ratio of total genetic variance over phenotypic variance, T2
      column=1: optimum number of families for this row
      column=2: optimum family size for this row
      column=3: standard error for this row
}
\seealso{
  \code{\link{vvarb}}
}
\references{Bijma, 2010. Estimating social genetic effects: precision of 
   estimates and optimum designs. Submitted to Genetics March 2010}
\author{P.Bijma; piter.bijma@wur.nl}
\examples{
  T <- 2000
  vpd <- 100
  vps <- 10
  h2d <- 0.3
  h2s <- 0.3
  rgds <- 0.1
  reds <- 0.1
  nw <- 4
  r <- 0.5
  optima <- nopt_2families(T,vpd,vps,h2d,h2s,rgds,reds,nw,r)
}
\keyword{utilities}

