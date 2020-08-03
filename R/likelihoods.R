# This file contains functions to calculate likelihoods and priors during
# tripartite record linkage.

# Function to calculate the likelihood of new data from given parameter values
# cmp.1to3 and cmp.2to3 are the comparison matrices from the output of
# BRL::compareRecords
calc.log.lkl <- function(cmp.1to3, cmp.2to3, n1, n2, n3, m, u, Z, Z2) {
  # 1. Which rows of cmp.1to3 and cmp.2to3 correspond to matches?
  matchrows.1to3 <- (which(Z2 <= n1) - 1) * n1 + Z2[Z2 <= n1]
  #                 ((what records in file 3 have links in file 1?) - 1) * n1
  #                 + (what records in file 1 are they linked to?)
  matchrows.2to3 <- (which((Z2 > n1) & (Z2 <= n1+n2)) - 1) * n2 + (Z2[(Z2 > n1) & (Z2 <= n1+n2)] - n1)
  #                 ((what records in file 3 have links in file 2?) - 1) * n2
  #                 + (what records in file 2 are they linked to?)
  # 2. Which rows of cmp.2to3 are not candidate links?
  # Any pairs (i,j) where i is a record in file 2 with a link in file 1
  noncand.1to3 <- c(outer((seq_len(n3) - 1) * n1, Z[Z <= n1], "+"))
  # Multiply each comparison matrix by log(m) and log(u), filter rows later
  lmprob.1to3 <- c(cmp.1to3 %*% log(m))
  luprob.1to3 <- c(cmp.1to3 %*% log(u))
  lmprob.2to3 <- c(cmp.2to3 %*% log(m))
  luprob.2to3 <- c(cmp.2to3 %*% log(u))
  # Filter to matches and nonmatches, appropriately, add the log probabilities of each
  nonmatchrows.1to3 <- setdiff(seq_len(nrow(cmp.1to3)), c(matchrows.1to3, noncand.1to3))
  nonmatchrows.2to3 <- setdiff(seq_len(nrow(cmp.2to3)), matchrows.2to3)
  loglkl <- (
    sum(lmprob.1to3[matchrows.1to3])
    + sum(luprob.1to3[nonmatchrows.1to3])
    + sum(lmprob.2to3[matchrows.2to3])
    + sum(luprob.2to3[nonmatchrows.2to3])
  )
  # Exponentiate the log likelihood and return
  return(loglkl)
}

# Function to calculate the likelihood of new data from given parameter values
# This likelihood function will trace links transitively instead of excluding
# non-candidates.
calc.log.lkl.tracing <- function(cmp.1to3, cmp.2to3, n1, n2, n3, m, u, Z, Z2) {
  # 1. Which rows of cmp.2to3 correspond to matches?
  matchrows.2to3 <- (which((Z2 > n1) & (Z2 <= n1+n2)) - 1) * n2 + (Z2[(Z2 > n1) & (Z2 <= n1+n2)] - n1)
  #                 ((what records in file 3 have links in file 2?) - 1) * n2
  #                 + (what records in file 2 are they linked to?)
  # 2. Which rows of cmp.1to3 correspond to matches (traced or direct?) matches?
  # First, create a copy of Z2 except following traced links where applicable.
  # I.e., any links to file 2 should be replaced with a record in file 1, if the linked
  # records in file 2 are further linked to file 1. Take advantage of the fact that
  # unlinked records essentially point to themselves.
  Z2.traced <- Z2
  Z2.traced[which((Z2 > n1) & (Z2 <= n1+n2))] <- Z[Z2[which((Z2 > n1) & (Z2 <= n1+n2))] - n1]
  # Now, use this new Z2 copy to find match rows
  matchrows.1to3 <- (which(Z2.traced <= n1) - 1) * n1 + Z2.traced[Z2.traced <= n1]
  #                 ((what records in file 3 have links in file 1?) - 1) * n1
  #                 + (what records in file 1 are they linked to?)
  # Multiply each comparison matrix by log(m) and log(u), filter rows later
  lmprob.1to3 <- c(cmp.1to3 %*% log(m))
  luprob.1to3 <- c(cmp.1to3 %*% log(u))
  lmprob.2to3 <- c(cmp.2to3 %*% log(m))
  luprob.2to3 <- c(cmp.2to3 %*% log(u))
  # Filter to matches and nonmatches, appropriately, add the log probabilities of each
  nonmatchrows.1to3 <- setdiff(seq_len(nrow(cmp.1to3)), matchrows.1to3)
  nonmatchrows.2to3 <- setdiff(seq_len(nrow(cmp.2to3)), matchrows.2to3)
  loglkl <- (
    sum(lmprob.1to3[matchrows.1to3])
    + sum(luprob.1to3[nonmatchrows.1to3])
    + sum(lmprob.2to3[matchrows.2to3])
    + sum(luprob.2to3[nonmatchrows.2to3])
  )
  # Exponentiate the log likelihood and return
  return(loglkl)
}

# Function to calculate the log of the prior of Z2 given Z and the hyperparameters
# aBM and bBM. If the value of Z2 is impossible given the value of Z (i.e.
# due to linking to a non-candidate) this returns -Inf := log(0)
# In R, negative infinity plus any finite value equals negative infinity, and all
# finite values compare greater than negative infinity. So this value can be used
# in acceptance ratio computations without worrying.
calc.log.Z2prior <- function(n1, n2, n3, Z2, Z, aBM, bBM) {
  # Candidates are any unlinked entries in file 1, plus all entries in file 2
  cand <- c(setdiff(seq_len(n1), Z), n1 + seq_len(n2))
  noncand <- Z[Z <= n1]
  # How many links are there currently between file 3 and the candidate set?
  nlinked <- sum(Z2 <= n1 + n2)
  # Calculate the log of the prior and return
  if (any(Z2 %in% noncand)) {
    # Impossible -> probability zero -> log probability -Inf.
    return(-Inf)
  } else {
    # Work in log scale: replacing log(x!) with lgamma(x+1) and log(beta(a,b))
    # with lbeta(a,b)
    logprior <- (
      (lgamma(length(cand) - nlinked + 1) - lgamma(length(cand) + 1))
      + (lbeta(aBM + nlinked, bBM + n3 - nlinked) - lbeta(aBM, bBM))
    )
    return(logprior)
  }
}