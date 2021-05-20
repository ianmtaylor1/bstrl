# This file contains debug functions: temporary functions that either expose
# some inner workings of the package or are copied and pasted segments of code
# that can be exported and run separately for debugging.

#' @export
Z_fc_smcmc_weights <- function(j, Z, Z2, m, u, cmpdata, aBM, bBM, directratio=TRUE) {
  n1 <- cmpdata[[1]][[1]]$n1

  weights <- rep(-Inf, n1+1)
  for (i in seq_len(n1+1)) {
    if (i > n1) {
      rec.i <- n1+j
    } else {
      rec.i <- i
    }
    Z.prop <- Z
    Z.prop[j] <- rec.i
    # First, the priors of Z and Z2
    weights[i] <- (
      calc.log.Z2prior(n1, Z2, Z.prop, aBM, bBM)
      + calc.log.Z2prior(n1, Z.prop, c(), aBM, bBM)
    )
    # If the prior has nonzero probability, calculate likelihoods
    # Since anything + -Inf == -Inf in R, this just saves time.
    if (weights[i] > -Inf) {
      weights[i] <- weights[i] + if (directratio) {
        smcmc.log.lkl.ratio.Z(cmpdata, m, u, Z.prop, Z2, j)
      } else {
        calc.log.lkl.tracing(cmpdata[[1]], m, u, c(), Z.prop) + calc.log.lkl.tracing(cmpdata[[2]], m, u, Z.prop, Z2)
      }
    }
  }
  # Semi-normalize, just to avoid too much over/underflow in exp
  weights <- weights - max(weights)

  return(weights)
}

#' @export
Z2_fc_smcmc_weights <- function(j, Z, Z2, m, u, cmpdata, aBM, bBM, directratio=TRUE) {
  # This function is copied and pasted from r_Z2_fc_smcmc
  n1 <- cmpdata[[1]][[1]]$n1
  n2 <- cmpdata[[1]][[1]]$n2

  weights <- rep(-Inf, n1+n2+1)
  # All n1 + n2 linked possibilities
  for (i in seq_len(n1+n2+1)) {
    if (i > n1 + n2) {
      rec.i <- n1+n2+j
    } else {
      rec.i <- i
    }
    Z2.prop <- Z2
    Z2.prop[j] <- rec.i
    # First, the priors of Z and Z2
    weights[i] <- calc.log.Z2prior(n1, Z2.prop, Z, aBM, bBM)
    # If the prior has nonzero probability, calculate likelihoods
    # Since anything + -Inf == -Inf in R, this just saves time.
    if (weights[i] > -Inf) {
      weights[i] <- weights[i] + if (directratio) {
        smcmc.log.lkl.ratio.Z2(cmpdata, m, u, Z, Z2.prop, j)
      } else {
        calc.log.lkl.tracing(cmpdata[[2]], m, u, Z, Z2.prop)
      }
    }
  }
  # Semi-normalize, just to avoid too much over/underflow in exp
  weights <- weights - max(weights)

  return(weights)
}
