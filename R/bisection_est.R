# ------------------------------
# Bisection iteration Method
# Solve: f(x) = 0
# ------------------------------
##............. Bisection Method
bisection_est <- function(fun, a, b, tol = 1e-6, max_iter = 10000, ...) {

  fa <- fun(a, ...)
  fb <- fun(b, ...)

  # check for sign change
  if (fa * fb > 0) {
    stop("Bisection requires f(a) and f(b) to have opposite signs.")
  }

  iter <- 0

  while ((b - a) / 2 > tol && iter < max_iter) {
    c <- (a + b) / 2         # midpoint
    fc <- fun(c, ...)

    # root found exactly
    if (fc == 0)
      return(c)

    # choose the subinterval with sign change
    if (fa * fc < 0) {
      b  <- c
      fb <- fc
    } else {
      a  <- c
      fa <- fc
    }

    iter <- iter + 1
  }

  return((a + b) / 2)        # approximate root
}

