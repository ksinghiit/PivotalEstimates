\name{PivotalEstimates}
\alias{PivotalEstimates}
\title{
Pivotal-Based Estimation under Block Progressive Censoring
}
\description{
Computes pivotal-based estimates of the common shape parameter and
block-specific and pooled scale parameters for progressively censored
samples under a block progressive censoring scheme.
}
\usage{
PivotalEstimates(aa, bb, rr, gt, mm, nn)
}
\arguments{
  \item{aa}{
Lower bound of the search interval for the shape parameter \eqn{\alpha} of two-parameter Weibull distribution
used in the bisection method.
}
  \item{bb}{
Upper bound of the search interval for the shape parameter \eqn{\alpha}
used in the bisection method.
}
  \item{rr}{
A list of progressive censoring schemes. Each element is a vector
\eqn{R_i} corresponding to the censoring scheme of the \eqn{i}-th block.
}
  \item{gt}{
A list of observed progressively censored samples. Each element
corresponds to the observed failure times from a block under Weibull model.
}
  \item{mm}{
A numeric vector containing the effective sample sizes
\eqn{(s_1, s_2, \ldots, s_k)} for each block.
}
  \item{nn}{
A numeric vector containing the initial sample sizes
\eqn{(n_1, n_2, \ldots, n_k)} for each block.
}
}
\details{
The function implements a pivotal-based inferential procedure for
block progressively censored data. A pivotal quantity for the common
shape parameter \eqn{\alpha} is constructed and solved using the
bisection method. Conditional on the estimated \eqn{\alpha}, block-wise
pivotal quantities are used to obtain estimates of the scale parameters
\eqn{\beta_i}. A pooled estimate of the scale parameter is obtained by
inverse-variance weighting across blocks.

Monte Carlo sampling is used to approximate the pivotal distributions,
with an initial burn-in period discarded to stabilize the estimates.
}
\value{
A list with the following components:
\item{alpha}{Pivotal estimate of the common shape parameter \eqn{\alpha}.}
\item{beta_i}{Numeric vector of pivotal estimates of the block-specific
scale parameters \eqn{\beta_1, \beta_2, \ldots, \beta_k}.}
\item{beta}{Pivotal estimate of the pooled scale parameter \eqn{\beta}.}
}
\references{
Balakrishnan, N. and Aggarwala, R. (2000).
\emph{Progressive Censoring: Theory, Methods, and Applications}.
Birkh\"auser, Boston.

Balakrishnan, N. (2007).
Progressive censoring methodology: An appraisal.
\emph{Test}, 16, 211--296.
}
\author{
Kundan Singh
}
\note{
The function uses Monte Carlo simulation with burn-in. The number of
iterations and burn-in length are controlled internally by the global
constants \code{NN} and \code{N0}. The censoring schemes must satisfy
\eqn{\sum R_i + s_i = n_i} for each block.
}
\seealso{
\code{\link{ProgressiveSample}}, \code{\link{bisection_est}}, \code{\link{rchisq}}
}
\examples{
# Example assumes progressively censored block data already generated

n <- 100 #sample slize
N <- c(40, 33, 27) # n1, n2, n3 ...each block sample size
M <- c(35, 30, 25) # s1, s2, s3 ...each block effective sample size
k <- length(M) #Block size

#Initial Value of Parameters
alp <- 1.7; bet <- 1.5;


NN=2000; #Total generate the pivotal estimates
N0=1000; #Consider the burn in obtained initial pivotal estimates

#Defined Censoring Schemes
R <- list()
for(r in 1:k){
  R[[r]] <- c(rep(0,(M[r]-1)),(N[r]-M[r])) # Censoring Schemes
}

# -------------------------
# Weibull quantile function
# F^{-1}(u) = [ -log(1-u) / beta ]^{1/alpha}
# -------------------------
q_weibull <- function(u, alpha, beta) {
  ((-log(1 - u)) / beta)^(1 / alpha)
}

# -------------------------
# Progressive censored sample generator
# qfunc : quantile function F^{-1}(u, alpha, beta)
# alpha, beta : parameters
# N : total sample size (not directly used but included for clarity)
# M : number of observed failures
# R : censoring scheme vector of length M
# -------------------------

# -------------------------
# Generate block progressively censored sample
# -------------------------


smpl_list <- list()
for(sm in 1:k){
  smpl_list[[sm]] <- ProgressiveSample(
    qfunc = q_weibull,
    alpha = alp,
    beta  = bet,
    N     = N[sm],
    M     = M[sm],
    R     = R[[sm]]
  ) + (rnorm(M[sm], 0, 0.001))
}
x <- smpl_list; #x;


gtt <- list()
for(i in 1:k){
  gtt[[i]] <- x[[i]];#Weibull distribution
}
#gtt;

gx <- function(xt){
  return(xt)#weibull
}

# Obtain pivotal estimates
pvt_val <- PivotalEstimates(
  aa = 1.2,
  bb = 2.3,
  rr = R,
  gt = x,
  mm = M,
  nn = N
)

cat('\n pivotal Estimates of alpha...\n', pvt_val[[1]])

cat('\n pivotal Estimates of beta_i...\n', pvt_val[[2]])

cat('\n pivotal Estimates of beta...\n', pvt_val[[3]])

}

