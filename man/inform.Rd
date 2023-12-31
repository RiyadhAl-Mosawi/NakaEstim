\name{inform}
\alias{inform}
\title{Information matrix using missing information principle of Liouse}
\description{
}
\usage{
inform(para, X, R, n, k, lower, upper)
}
\arguments{
  \item{para}{(shape,scale) parameters of Nakagami distribution}
  \item{X}{Failure times}
  \item{R}{The vector of removals scheme}
  \item{n}{Sample size i.e the multiplication of number of groups(m) by size of each group (k)}
  \item{k}{size of group}
  \item{lower}{lower limit of integral}
  \item{upper}{upper limit of integral}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (para, X, R, n, k, lower, upper) 
{
    .Call(`_NDLibraryFunction_inform`, para, X, R, n, k, lower, 
        upper)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
