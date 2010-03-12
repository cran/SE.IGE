\name{vvarb}
\alias{vvarb}
\title{Sampling variance of estimated between-family variance}
\description{
  Calculates the sampling variance of an estimated
  between-family variance in a one-way ANOVA.}
\usage{vvarb(var,varb,N,n)}
\arguments{
  \item{var}{full (or phenotypic) variance, numeric}
  \item{varb}{between-family variance, numeric}
  \item{N}{number of families, numeric}
  \item{n}{family size, numeric}
}
\value{sampling variance of between family variance, numeric}
\references{e.g. Lynch and Walsh, 1998, Sinauer, p. 561}
\author{P. Bijma; piter.bijma@wur.nl}
\examples{
var <- 100
varb <- 10
N <- 100
n <- 20
var_varb <- vvarb(var,varb,N,n)
}
\keyword{utilities}

