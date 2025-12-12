\name{ProgressiveSample}
\alias{ProgressiveSample}
\title{
Generation of Progressively Type-II Censored Samples
}
\description{
Generates a progressively Type-II censored sample from a continuous lifetime
distribution using a quantile-based simulation approach.
}
\usage{
ProgressiveSample(qfunc, alpha, beta, N, M, R)
}
\arguments{
  \item{qfunc}{
A quantile function corresponding to the lifetime distribution.
It should take arguments \code{u}, \code{alpha}, and \code{beta}.
}
  \item{alpha}{
The shape parameter of the lifetime distribution.
}
  \item{beta}{
The scale parameter of the lifetime distribution.
}
  \item{N}{
Total number of units placed on test.
}
  \item{M}{
Number of observed failures.
}
  \item{R}{
A progressive censoring scheme vector of length \code{M}, where \code{R[i]}
denotes the number of surviving units removed immediately after the \code{i}-th failure.
}
}
\details{
In a progressive Type-II censoring experiment, a total of \code{N} units are
placed on test and only \code{M} failures are observed. After each observed
failure, a predetermined number of surviving units are removed according to
the censoring scheme \code{R}. The function implements a general simulation
algorithm based on the inverse distribution (quantile) method.
}
\value{
A numeric vector of length \code{M} containing the observed failure times under
the specified progressive censoring scheme.
}
\references{
Balakrishnan, N. and Aggarwala, R. (2000).
\emph{Progressive Censoring: Theory, Methods, and Applications}.
Birkh\"auser, Boston.
}
\author{
Kundan Singh
}
\note{
The censoring scheme must satisfy the condition \code{sum(R) + M = N}.
}
\seealso{
\code{\link{runif}}, \code{\link{qweibull}}
}
\examples{
# Weibull quantile function
q_weibull <- function(u, alpha, beta) {
  ((-log(1 - u)) / beta)^(1 / alpha)
}

# Total and effective sample sizes
N <- 50
M <- 35

# True parameter values
alpha <- 1.7
beta  <- 1.5

# Progressive censoring scheme
R <- c(rep(0, M - 1), N - M)

# Generate progressively censored sample
xsamp <- ProgressiveSample(
  qfunc = q_weibull,
  alpha = alpha,
  beta  = beta,
  N     = N,
  M     = M,
  R     = R
)

xsamp
}

