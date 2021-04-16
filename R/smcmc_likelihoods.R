# This file contains functions to calculate optimized and simplified likelihoods
# for use with SMCMC sampling.

# future note: the goal in generalizing these functions to the streaming case is
# to have a single function that can work for the specified component *of the
# specified Z vector*.
#
# But that is TBD


# Compute the ratio of two likelihoods for nearly identical states of Z.
# The difference is in the single component j, 1 <= j <= n2. We compare the
# ratio for the state of Z as it is supplied, vs the state of Z where the
# component j is unlinked. If component j is already unlinked in the value of Z
# that is passed in, then we return a value of 0 (log of a ratio of 1).
smcmc.log.lkl.ratio.Z <- function(cmpdata, m, u, Z, Z2, j) {
  n1 <- cmpdata[[1]][[1]]$n1
  n2 <- cmpdata[[2]][[1]]$n2
  n3 <- cmpdata[[2]][[2]]$n2
  # If Z[j] is unlinked don't need to do anything
  if (Z[j] == n1 + j) {
    return(0)
  } else {
    totals <- rep(0, length(m))
    # First add the comparison to file 1.
    totals <- totals + comparison(cmpdata, 1, 2, Z[j], j)
    # If anything in file 3 is linked here, add that follow-on link
    k <- match(n1 + j, Z2)
    if (!is.na(k)) {
      totals <- totals + comparison(cmpdata, 1, 3, Z[j], k)
    }
    # we want the ratio of m to u for every level of disagreement recorded
    return(sum(log(m / u) * totals))
  }
}

# Compute the ratio of two likelihoods for nearly identical states of Z2.
# The difference is in the single component j, 1 <= j <= n3. We compare the
# ratio for the state of Z2 as it is supplied, vs the state of Z2 where the
# component j is unlinked. If component j is already unlinked in the value of Z2
# that is passed in, then we return a value of 0 (log of a ratio of 1).
smcmc.log.lkl.ratio.Z2 <- function(cmpdata, m, u, Z, Z2, j) {
  n1 <- cmpdata[[1]][[1]]$n1
  n2 <- cmpdata[[2]][[1]]$n2
  n3 <- cmpdata[[2]][[2]]$n2
  # If Z2[j] is unlinked we don't even need to proceed
  if (Z2[j] == n1+n2+j) {
    return(0)
  } else {
    totals <- rep(0, length(m))
    #Where does Z2[j] point?
    if (Z2[j] <= n1) {
      # File 1: only need to add that comparison
      total <- totals + comparison(cmpdata, 1, 3, Z2[j], j)
    } else {
      # File 2: need to add that comparison and any follow-on comparisons in file 1
      totals <- totals + comparison(cmpdata, 2, 3, Z2[j] - n1, j)
      if (Z[Z2[j] - n1] <= n1) {
        totals <- totals + comparison(cmpdata, 1, 3, Z[Z2[j] - n1], j)
      }
    }
    # we want the ratio of m to u for every level of disagreement recorded
    return(sum(log(m / u) * totals))
  }
}
