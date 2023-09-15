\name{BS_MPS}
\alias{BS_MPS}
\title{Bootsrap samples based on product of spacing function of Nakagami distribution under first-failure type I progressive censoring data}
\description{
Bootsrap samples based on product of spacing function of Nakagami distribution under first-failure type I progressive censoring data}
\usage{
BS_MPS(para, R, lw, up, l, u, t, Boot, k, verbose=0, display_progress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{para}{(shape,scale) parameters of Nakagami distribution}
  \item{R}{The vector of removals scheme}
  \item{lw}{The vector of lower limits of the parameters}
  \item{up}{The vector of upper limits of the parameters}
  \item{l}{: The lower specification limit}
  \item{u}{: The upper specification limit}
  \item{t}{: The target value}
  \item{Boot}{The size of bootstrap sample}
  \item{k}{Size of each group}
  \item{verbose}{Non-zero: to display the results of each iteration}
  \item{display_progress}{TRUE: to display the progress bar of the process}
}
\details{Obtaing a parametric bootstrap samples for the first-failure Nakagami distribution based on product of spacing function under first-failure type I progressive censoring data}
\value{
A matrix with 6 columns and Boots rows. The first three columms represnt the estiamtes of the shape, scale and Spmk paratemers, respectively. The next three columns represnt the estimated standrad errors using observed information matrix of the shape, scale and Spmk paratemers, respectively.
}
\author{Riyadh Al-Mosawi}
\examples{BS_MPS(c(1,1), c(0,0,0,0,10), c(0.5,0),c(5,5), 1, 10, 1, 100, 2, 0, TRUE) }
