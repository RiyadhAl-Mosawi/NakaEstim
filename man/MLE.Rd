\name{MLE}
\alias{MLE}
\title{Compute the MLE using Nelder-Mead algorithm of Nakagami distribution under first-failure type I progressive censoring data}
\description{
Compute the MLE using Nelder-Mead algorithm of Nakagami distribution under first-failure type I progressive censoring data. Information matrix using observed information principle is also computed.
}
\usage{
MLE(para, X, R, k)
}
%- maybe also 'usage' for other objects documented here.
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
MLE(para, X, R, k) 
{
    .Call(`_NDLibraryFunction_MLE`, para, X, R, k)
  }
}
