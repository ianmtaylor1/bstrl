# This file contains componentwise full conditional samplers for Z and Z2
# 3/23/2021 - Currently used with the SMCMC transition and approximate
# jumping kernels. To do: make flexible enough to deal with more than 3 files


# Z, Z2, m, u - current parameter values
# cmpdata - a list-of-lists format of comparison data
# aBM, bBM - Z/Z2 prior parameters.
#
# Returns: the new value of Z
r_Z_fc_smcmc <- function(Z, Z2, m, u, cmpdata, aBM, bBM) {
  n1 <- cmpdata[[1]][[1]]$n1

  for (j in 1:length(Z)) {
    weights <- rep(0, n1+1)
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
        weights[i] <- weights[i] + (
          calc.log.lkl.tracing(cmpdata[[1]], m, u, c(), Z.prop)
          + calc.log.lkl.tracing(cmpdata[[2]], m, u, Z.prop, Z2)
        )
      }
    }
    # Semi-normalize, just to avoid too much over/underflow in exp
    weights <- weights - max(weights)
    # Select new value, put it directly in Z
    Z[j] <- sample(c(seq_len(n1), n1+j), 1, prob=exp(weights))
  }

  return(Z)
}

# Z, Z2, m, u - current parameter values
# cmpdata - a list-of-lists format of comparison data
# aBM, bBM - Z/Z2 prior parameters.
#
# Returns: the new value of Z2
r_Z2_fc_smcmc <- function(Z, Z2, m, u, cmpdata, aBM, bBM) {
  n1 <- cmpdata[[1]][[1]]$n1
  n2 <- cmpdata[[1]][[1]]$n2

  for (j in 1:length(Z2)) {
    weights <- rep(0, n1+n2+1)
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
        weights[i] <- weights[i] + calc.log.lkl.tracing(cmpdata[[2]], m, u, Z, Z2.prop)
      }
    }
    # Semi-normalize, just to avoid too much over/underflow in exp
    weights <- weights - max(weights)
    # Select new value, put it directly in Z
    Z2[j] <- sample(c(seq_len(n1+n2), n1+n2+j), 1, prob=exp(weights))
  }

  return(Z2)
}
