\name{SE_corr}
\alias{SE_corr}
\title{standard error of estimated correlation}
\description{predicts the standard error of an estimated genetic correlation
  making use of the sampling variance of the (co)variances of both parameters}
\usage{SE_corr(v1,v2,c12,SEv1,SEv2,SEc12)}
\arguments{
  \item{v1}{true variance of first variable, numeric} 
  \item{v2}{true variance of second variable, numeric} 
  \item{c12}{true covariance between first and second variable}
  \item{SEv1}{standard error of estimated variance of first variable,
    numeric}
  \item{SEv2}{standard error of estimated variance of second variable,
    numeric}
  \item{SEc12}{standard error of estimated covariance between
    first and second variable, numeric}
}
\value{standard error of the estimated genetic correlation, numeric}
\references{Bijma, 2010. Estimating social genetic effects: precision of 
   estimates and optimum designs. Submitted to Genetics March 2010}
\author{P.Bijma; piter.bijma@wur.nl}
\examples{
  v1 <- 100
  v2 <- 400
  c12 <- 130
  SEv1 <- 35
  SEv2 <- 80
  SEc12 <- 50
  Serror_corr <- SE_corr(v1,v2,c12,SEv1,SEv2,SEc12)
}
\keyword{utilities}


