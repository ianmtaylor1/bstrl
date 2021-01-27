# This file contains functions to calculate likelihoods and priors during
# tripartite record linkage.

# Function to calculate the likelihood of new data from given parameter values.
# Operates on comparisons between a file k and all k-1 previous files (k-1 total
# comparison objects)
# Parameters:
#   cmpdata - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#   m,u,Z,Z2 - parameter values. Z2 is the links originating from file k, and Z
#              is the links originating from files 2 through k-1, concatenated
#              together (file 1 is excluded because it's unnecessary under the
#              assumption of no duplicates within files.)
#   do.trace - Whether to perform link tracing
# Notes:
#   Eventually, I want to deprecate calc.log.lkl.tracing in favor of just using
#   this function with the flag. But it's fine for now.
calc.log.lkl <- function(cmpdata, m, u, Z, Z2, do.trace=FALSE) {
  # Calculate the disagreement level counts
  counts <- disag.counts(cmpdata, Z, Z2, do.trace)
  # The log likelihood is these counts combined with log m and log u
  return(sum(counts$match * log(m)) + sum(counts$nonmatch * log(u)))
}

# Function to calculate the likelihood of new data from given parameter values
# This likelihood function will trace links transitively instead of excluding
# non-candidates.
# Parameters:
#   cmpdata - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#   m,u,Z,Z2 - parameter values. Z2 is the links originating from file k, and Z
#              is the links originating from files 2 through k-1, concatenated
#              together (file 1 is excluded because it's unnecessary under the
#              assumption of no duplicates within files.)
calc.log.lkl.tracing <- function(cmpdata, m, u, Z, Z2) {
  return(calc.log.lkl(cmpdata, m, u, Z, Z2, do.trace=TRUE))
}

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

# Function to calculate the log prior for Z2 when using the empirical prior
# where the probability of a pair being linked is proportional to the number
# of fields which compare equal (plus a base weight). If the value of Z2 is
# impossible given the value of Z (i.e. due to linking to a non-candidate) this
# returns -Inf := log(0)
# In R, negative infinity plus any finite value equals negative infinity, and all
# finite values compare greater than negative infinity. So this value can be used
# in acceptance ratio computations without worrying.
# Parameters:
#       cmpdata - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#       Z,Z2 - parameter values. Z2 is the links originating from file k, and Z
#              is the links originating from files 2 through k-1, concatenated
#              together (file 1 is excluded because it's unnecessary under the
#              assumption of no duplicates within files.)
#       aBM, bBM - prior parameters for beta link prior
#       baseweight - base amount added to the number of equal fields when
#              computing link weight
# Notes:
#   Assumes length(Z2) <= n1 + length(Z) !!!
calc.log.Z2prior.empirical <- function(cmpdata, Z, Z2, aBM, bBM, baseweight=0.01) {
  # Calculates in the order of:
  # 1. Probability of the observed "linked" indicators (beta-binomial-ish)
  # 2. Probability of each record's state, in order, conditioned on the
  #    indicators and the previous records' states. (Links weighted empirically)

  # Required file sizes
  k <- length(cmpdata) + 1 # Number of total files
  ns <- rep(0, k) # Vector of length k of file sizes
  for (file in 1:(k-1)) {
    ns[file] <- cmpdata[[file]]$n1
  }
  ns[k] <- cmpdata[[1]]$n2
  n1 <- ns[1]
  nlast <- ns[k]
  nprev <- sum(ns[1:(k-1)]) # Or nprev <- n1 + length(Z)

  if (valid.link.state(n1, Z, Z2)) {
    # Probaility of these indicators
    nlinked <- sum(Z2 <= nprev)
    logprob <- lbeta(aBM + nlinked, bBM + nlast - nlinked) - lbeta(aBM, bBM)
    # Add on conditional log probabilities of each observed pair
    # (conditional probability of unlinked is always 1)
    which.linked <- which(Z2 <= nprev)
    # Store candidate records. Will change as links are created
    cand <- setdiff(seq_len(nprev), Z[Z < n1 + seq_len(length(Z))])
    for (recj in which.linked) {
      # Vector to store the weights based on field agreement
      weights <- rep(0, length(cand))
      for (file in seq_len(k-1)) {
        # Find the candidates that are in this file
        filestart <- sum(ns[seq_len(file-1)]) + 1
        fileend <- filestart + ns[file] - 1
        in.file <- (cand >= filestart) & (cand <= fileend)
        # Pull out the number of matching fields for those candidates
        pairindices <- (recj - 1) * ns[file] + (cand[in.file] - filestart + 1) # (j-1)*n1 + i
        weights[in.file] <- attr(cmpdata[[file]]$comparisons, "numfieldsequal")[pairindices]
      }
      weights <- (weights + baseweight) / sum(weights + baseweight)
      logprob <- logprob + log(weights[match(Z2[recj], cand)])
      cand <- setdiff(cand, Z2[recj])
    }
    return(logprob)
  } else {
    return(-Inf)
  }
}

################################################################################
# Helper Functions for Likelihood/Prior functions ##############################
################################################################################

# Determines whether the given values of Z, Z2 are viable according to the
# candidate set rules. I.e. returns true if no records are linked to the same
# record in file 1. Else returns false.
# Parameters:
#   Z2 = links involving the most recent file
#   Z = links between all previous files
#   offset = An offset of the starting position of the indices of Z. Any records
#            with indices less than or equal to offset are assumed to be the end
#            of their chain. This can be used to exclude file 1, for example, if
#            its links are excluded from Z to save space
# Notes: assumption is that if a record j is unlinked then Z[j - offset] == j
valid.link.state <- function(offset, Z, Z2) {
  noncand <- Z[Z < offset + seq_len(length(Z))] # Non-candidates: records with links in later files
  return(!any(Z2 %in% noncand))
}

