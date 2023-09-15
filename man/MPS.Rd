\name{MPS}
\alias{MPS}
\title{Compute the MPS using Nelder-Mead algorithm of Nakagami distribution under first-failure type I progressive censoring data}
\description{
Compute the MPS using Nelder-Mead algorithm of Nakagami distribution under first-failure type I progressive censoring data. Information matrix using observed information principle is also computed.
}
\usage{
MPS(para, X, R, k)
}
\arguments{
  \item{para}{(shape,scale) parameters of Nakagami distribution}
  \item{X}{Failure times}
  \item{R}{The vector of removals scheme}
  \item{k}{size of group}
}
\value{
 returns a list which includes
  \item{par}{estimated values}
  \item{value}{value of log-likelihood function at the estimated value}
  \item{inform}{infomation matrix using observed information principle}
}
\author{
Riyadh Al-Mosawi
}
\examples{
MPS(para, X, R, k) 
}
