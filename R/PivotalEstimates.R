#' @importFrom ProgressiveSample ProgressiveSample
#' @importFrom FixedPoint fixed_point_est
#' @export fixed_point_est

#' @importFrom FixedPoint bisection_est
#' @export bisection_est
#' devtools::document()

#' Pivotal-based parameter estimation for block progressively censored data
#'
#' This function computes pivotal-based estimates of the common shape parameter
#' \eqn{\alpha} and the block-specific scale parameters \eqn{\beta_i} for
#' block progressively censored samples (BPCS). A pivotal equation for
#' \eqn{\alpha} is solved using the bisection method, and conditional on the
#' resulting \eqn{\alpha}, pivotal quantities are used to obtain estimates of
#' \eqn{\beta_i}. A pooled estimate of \eqn{\beta} is also computed via
#' inverse-variance weighting across blocks. Monte Carlo simulation is used,
#' with a burn-in period discarded to stabilize the estimates.
#'
#' @param aa Lower bound of the search interval for \eqn{\alpha}.
#' @param bb Upper bound of the search interval for \eqn{\alpha}.
#' @param rr A list of progressive censoring schemes, one vector per block.
#' @param gt A list of observed progressively censored samples, one vector per block.
#' @param mm A numeric vector of effective sample sizes \eqn{(s_1,\ldots,s_k)}.
#' @param nn A numeric vector of initial sample sizes \eqn{(n_1,\ldots,n_k)}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{alpha}}{Pivotal estimate of the common shape parameter \eqn{\alpha}.}
#'   \item{\code{beta_i}}{Vector of pivotal estimates of block-specific scale parameters \eqn{\beta_1,\ldots,\beta_k}.}
#'   \item{\code{beta}}{Pivotal estimate of the pooled scale parameter \eqn{\beta}.}
#' }
#'
#' @details The censoring schemes must satisfy \eqn{\sum_{j=1}^{s_i} R_{ij} + s_i = n_i}
#' for each block \eqn{i}. The number of Monte Carlo iterations and burn-in length
#' are controlled by the objects \code{NN} and \code{N0} used in the function environment.
#'
#' @seealso \code{\link{bisection_est}}, \code{\link{ProgressiveSample}}
#'
#' @export


##...pivotal based inference

PivotalEstimates <- function(aa, bb, rr, gt, mm, nn){

  pivalp <- function(ap){ # pivotal function for alpha
    apstr <- mean(rchisq(5,(2*sum(mm+1))));
    smi <- 0;
    for(i in 1:length(mm)){
      sm21 <- sm22 <- 0;
      for(l in 1:(mm[i]-1)){
        sm21 <- sm21 + ((rr[[i]][l])+1)*((gt[[i]][l])^ap);
        sm22 <- sm22 + ((rr[[i]][l])+1)
      }
      demit <- sm21 + ((nn[i]-sm22)*((gt[[i]][(mm[i])])^ap));#W_imi
      smj <- 0;
      for(j in 1:(mm[i]-1)){
        sm11 <- sm12 <- 0;
        for(l in 1:(j)){
          sm11 <- sm11 + ((rr[[i]][l])+1)*((gt[[i]][l])^ap);
          sm12 <- sm12 +	((rr[[i]][l])+1)
        }
        numrt <- sm11 + ((nn[i]-sm12)*((gt[[i]][j])^ap)); #W_ij, j=2,3,...,(mm[i]-1)
        smj <- smj + log(numrt/demit);
      }
      smi <- smi + smj;
    }
    piv_alp <- -2*(smi);
    pp <- piv_alp-apstr;
    return(pp);
  }

  #..pivotal_bet_i
  piv_bti <- function(ap, rr, gt, mm, nn){
    pivotbti <-  numeric();
    for(i in 1:length(mm)){
      btstr <- rchisq(1,(2*mm[i]));
      sm21 <- sm22 <- 0;
      for(l in 1:(mm[i]-1)){
        sm21 <- sm21 + ((rr[[i]][l])+1)*((gt[[i]][l])^ap);
        sm22 <- sm22 + ((rr[[i]][l])+1)
      }
      GTi <- 2*(sm21 + ((nn[i]-sm22)*((gt[[i]][mm[i]])^ap)));
      pivotbti[i] <- btstr/GTi;
    }
    return(pivotbti);
  }
  #iteration
  pvt_ap <- pvt_bt <- pvt_sf <- pvt_hf <- pvt_mdtf <- numeric();
  pvt_bti <- matrix(0, nrow=NN, ncol=length(M));
  for(ii in 1:NN){
    pvt_ap_val <- bisection_est(pivalp, a=aa, b=bb);
    pvt_ap[ii] <- pvt_ap_val; #pivotal alpha
    pvt_bti_val <- piv_bti(pvt_ap[ii], rr, gt, mm, nn);
    pvt_bti[ii,] <- pvt_bti_val;#pivotal beta_i
    varbti <- etinv <- numeric(); etbti <- 0;
    for(i in 1:length(M)){
      varbti[i] <- var(pvt_bti[,i])
      etinv[i] <- 1/varbti[i] ;
      etbti <- etbti + (etinv[i]*(pvt_bti[ii,i]));
    }
    pvt_bt[ii] <- sum(etbti)/sum(etinv); #pivotal beta
  }
  #pivotal estimates
  pvt_ap_est <- mean(pvt_ap[(N0+1):NN]); #pivt alp
  pvt_bti_est <- numeric()
  for(i in 1:k){
    pvt_bti_est[i] <- mean(pvt_bti[,i][(N0+1):NN]);#pivt beti
  }
  pvt_bt_est <- mean(pvt_bt[(N0+1):NN]); #pivt bet

  pivotal_est <- list()
  pivotal_est[[1]] <- pvt_ap_est;
  pivotal_est[[2]] <- pvt_bti_est;
  pivotal_est[[3]] <- pvt_bt_est;

  return(pivotal_est)
}


