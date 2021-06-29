# Function to calculate the log of the prior of Z2 given Z and the hyperparameters
# aBM and bBM. If the value of Z2 is impossible given the value of Z (i.e.
# due to linking to a non-candidate) this returns -Inf := log(0)
# In R, negative infinity plus any finite value equals negative infinity, and all
# finite values compare greater than negative infinity. So this value can be used
# in acceptance ratio computations without worrying.
# Parameters:
#       n1 - Number of records in file 1
#       Z,Z2 - parameter values. Z2 is the links originating from file k, and Z
#              is the links originating from files 2 through k-1, concatenated
#              together (file 1 is excluded because it's unnecessary under the
#              assumption of no duplicates within files.)
#       aBM, bBM - prior parameters for beta link prior
# Notes:
#   Assumes length(Z2) <= n1 + length(Z) !!!
calc.log.Z2prior <- function(n1, Z2, Z, aBM, bBM) {
  # Calculate the log of the prior and return
  if (valid.link.state(n1, Z, Z2)) {
    nprev <- length(Z) # Total records in previous files except file 1
    nlast <- length(Z2) # Records in the most recent file
    # Candidates are any unlinked entries in file 1, plus all entries in file 2
    ncand <- n1 + nprev - sum(Z < n1 + seq_len(nprev))
    # How many links are there currently between file 3 and the candidate set?
    nlinked <- sum(Z2 <= n1 + nprev)
    # Work in log scale: replacing log(x!) with lgamma(x+1) and log(beta(a,b))
    # with lbeta(a,b)
    logprior <- (
      (lgamma(ncand - nlinked + 1) - lgamma(ncand + 1))
      + (lbeta(aBM + nlinked, bBM + nlast - nlinked) - lbeta(aBM, bBM))
    )
    return(logprior)
  } else {
    # Impossible -> probability zero -> log probability -Inf.
    return(-Inf)
  }
}


calc.log.Z2prior.flat <- function(n1, Z2, Z) {
  if (check.no.duplicates(Z2)) {
    return(0)
  } else {
    return(-Inf)
  }
}


calc.log.Z2prior.noinvalid <- function(n1, Z2, Z, aBM, bBM) {
  # Minimal validity check: only check for duplicates within Z2, not compatibility
  # with Z1 and Z2 together
  if (check.no.duplicates(Z2)) {
    nprev <- length(Z) # Total records in previous files except file 1
    nlast <- length(Z2) # Records in the most recent file
    # Candidates in this implementation are all previous records: no non-candidates
    # are excluded, no invalid states are given zero probability in the prior
    ncand <- n1 + nprev
    # How many links are there currently between file 3 and the candidate set?
    nlinked <- sum(Z2 <= n1 + nprev)
    # Work in log scale: replacing log(x!) with lgamma(x+1) and log(beta(a,b))
    # with lbeta(a,b)
    logprior <- (
      (lgamma(ncand - nlinked + 1) - lgamma(ncand + 1))
      + (lbeta(aBM + nlinked, bBM + nlast - nlinked) - lbeta(aBM, bBM))
    )
    return(logprior)
  } else {
    return(-Inf)
  }
}
