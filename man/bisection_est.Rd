\name{bisection_est}
\alias{bisection_est}
\title{Bisection Method}
\description{
  Performs bisection method to solve f(x)=0.
}
\usage{
bisection_est(f, a, b, n = 1000, tol = 1e-7)
}
\arguments{
  \item{f}{The iteration function f(x).}
  \item{a, b}{Initial guess, where f(a)f(b)<0}
  \item{n}{Number of iteration}
}
\examples{
bisection_est(function(x) cos(x)-x, a=0, b=1)
bisection_est(function(x) sqrt(3+x)-x, a=0, b=3)
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
function (x)
{
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
