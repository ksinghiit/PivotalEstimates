````md
# PivotalEstimates
An R package for **pivotal-based parameter estimation** under
**Block Progressive Type-II Censoring Schemes (BPCS)**, with a focus on
the **Weibull distribution**.

---

## Overview

In block progressive censoring, experimental units are divided into
several homogeneous blocks (testing groups). Each block follows its own
progressive censoring scheme, while the lifetime distribution shares a
**common shape parameter** across blocks and allows **block-specific
scale parameters**.

This package implements a **Monte Carlo pivotal inference procedure**
to estimate:

- a common shape parameter \( \alpha \),
- block-specific scale parameters \( \beta_1, \ldots, \beta_k \),
- a pooled scale parameter \( \beta \).

---

## Statistical Model

The current implementation focuses on the **Weibull distribution**.
For block \( i \), lifetimes follow

\[
X_i \sim \text{Weibull}(\alpha, \beta_i),
\]

with survival function

\[
S(x) = \exp\left(-\beta_i x^{\alpha}\right).
\]

- \( \alpha \) is common across all blocks,
- \( \beta_i \) captures block-specific testing effects.

---

## Methods

The estimation procedure is based on **pivotal quantities**:

1. A pivotal equation for the common shape parameter \( \alpha \) is
   constructed from progressively censored data and solved using the
   **bisection method**.
2. Conditional on the estimated \( \alpha \), pivotal quantities are
   used to obtain estimates of each \( \beta_i \).
3. A pooled estimate of \( \beta \) is computed using
   **inverse-variance weighting**.
4. Monte Carlo simulation with burn-in is used to stabilize estimates.

---

## Package Dependencies

This package integrates functionality from:

- **ProgressiveSample** – for generating progressively censored samples
- **FixedPoint** – for numerical root finding (bisection method)

Dependencies are handled internally using proper namespace imports.

---

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ksinghiit/PivotalEstimates")
````

---

## Main Function

### `PivotalEstimates()`

```r
PivotalEstimates(aa, bb, rr, gt, mm, nn)
```

**Arguments**

* `aa`, `bb` : lower and upper bounds for solving the pivotal equation for ( \alpha )
* `rr` : list of progressive censoring schemes (one per block)
* `gt` : list of observed progressively censored samples (X) for Weibull distribution but it can be change for other Shape-Scale Family of distribution
* `mm` : vector of effective sample sizes
* `nn` : vector of initial sample sizes

**Returns**

A list containing:

* `alpha`  : pivotal estimate of the common shape parameter
* `beta_i` : pivotal estimates of block-specific scale parameters
* `beta`   : pooled pivotal estimate of the scale parameter

---

## Example Workflow (Weibull Model)

```r
library(PivotalEstimates)

# Block structure
N <- c(40, 33, 27)   # initial sample sizes
M <- c(35, 30, 25)   # effective sample sizes

# Progressive censoring schemes
R <- list(
  c(rep(0, M[1]-1), N[1]-M[1]),
  c(rep(0, M[2]-1), N[2]-M[2]),
  c(rep(0, M[3]-1), N[3]-M[3])
)

# Weibull quantile function
q_weibull <- function(u, alpha, beta) {
  ((-log(1 - u)) / beta)^(1 / alpha)
}

# Generate block-wise censored samples
x <- lapply(1:3, function(i)
  ProgressiveSample(q_weibull, alpha = 1.7, beta = 1.5,
                    N = N[i], M = M[i], R = R[[i]])
)

# Pivotal estimation
pvt <- PivotalEstimates(
  aa = 1.2,
  bb = 3,
  rr = R,
  gt = x,
  mm = M,
  nn = N
)

pvt
```

---

## Applications

* Reliability analysis
* Life-testing experiments
* Progressive and block censoring studies
* Monte Carlo inference under censoring

---

## Author

**Kundan Singh**

```

