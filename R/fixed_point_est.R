# ------------------------------
# Fixed point iteration
# Solve: x = g(x)
# ------------------------------
fixed_point_est <- function(gfun, init_val, tol = 0.001, max_iter = 1000, ...) {

  x <- numeric()           # store iterates like your code
  x[1] <- init_val
  itr <- 1

  while (itr <= max_iter) {

    # Next iteration
    x[itr + 1] <- gfun(x[itr])

    # Convergence check
    if (abs(x[itr + 1] - x[itr]) < tol) {
      return(x[itr + 1])   # final fixed point estimate
    }

    itr <- itr + 1
  }

  # If max_iter reached
  return(x[itr])
}
