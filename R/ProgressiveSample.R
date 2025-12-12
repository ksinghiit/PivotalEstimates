# -------------------------
# Progressive censored sample generator
# qfunc : quantile function F^{-1}(u, alpha, beta)
# alpha, beta : parameters
# N : total sample size (not directly used but included for clarity)
# M : number of observed failures
# R : censoring scheme vector of length M
# -------------------------
ProgressiveSample <- function(qfunc, alpha, beta, N, M, R) {

  if (length(R) != M) {
    stop("Length of R must be equal to M.")
  }

  # Step 1: generate uniform variables
  w  <- runif(M)
  v  <- numeric(M)
  u  <- numeric(M)

  # Step 2: construct v_i as in the progressive censoring algorithm
  for (i in 1:M) {
    denom   <- i + sum(R[(M - i + 1):M])
    v[i]    <- w[i]^(1 / denom)
  }

  # Step 3: construct u_i from v_i
  for (i in 1:M) {
    u[i] <- 1 - prod(v[(M - i + 1):M])
  }

  # Step 4: transform uniforms to Weibull via quantile function
  x <- qfunc(u, alpha, beta)
  return(x)
}
