# This file contains functions to calculate likelihoods and priors during
# tripartite record linkage.

# Function to calculate the likelihood of new data from given parameter values
# cmp.1to3 and cmp.2to3 are the comparison objects from the output of
# BRL::compareRecords
calc.log.lkl <- function(cmp.1to3, cmp.2to3, m, u, Z, Z2) {
  # 1. Which rows of cmp.1to3 and cmp.2to3 correspond to matches?
  matchrows.1to3 <- matchrows(cmp.1to3, Z2, offset=0)
  matchrows.2to3 <- matchrows(cmp.2to3, Z2, offset=cmp.1to3$n1)
  # 2. Which rows of cmp.2to3 are not candidate links?
  noncand.1to3 <- noncandrows(cmp.1to3, Z, offset=0)
  # 3. Compute filtered column sums: for match and nonmatch pairs, how many times
  #    does a given level of a given field appear?
  total.count <- attr(cmp.1to3$comparisons, "totals") + attr(cmp.2to3$comparisons, "totals")
  match.count <- colSums(cmp.1to3$comparisons[matchrows.1to3,,drop=FALSE]) + colSums(cmp.2to3$comparisons[matchrows.2to3,,drop=FALSE])
  nonmatch.count <- total.count - match.count - colSums(cmp.1to3$comparisons[noncand.1to3,,drop=FALSE])
  # 4. Calculate and return log likelihood based on count of field levels for
  #    matches and non-matches
  loglkl <- (
    sum(match.count * log(m))
    + sum(nonmatch.count * log(u))
  )
  return(loglkl)
}

# Function to calculate the likelihood of new data from given parameter values
# This likelihood function will trace links transitively instead of excluding
# non-candidates.
calc.log.lkl.tracing <- function(cmp.1to3, cmp.2to3, m, u, Z, Z2) {
  # 1. Which rows of cmp.2to3 correspond to matches?
  matchrows.2to3 <- matchrows(cmp.2to3, Z2, offset=cmp.1to3$n1)
  # 2. Which rows of cmp.1to3 correspond to matches (traced or direct?) matches?
  matchrows.1to3 <- matchrows(cmp.1to3, traced(Z2, Z, steps=1, offset=cmp.1to3$n1), offset=0)
  # 3. Compute filtered column sums: for match and nonmatch pairs, how many times
  #    does a given level of a given field appear?
  total.count <- attr(cmp.1to3$comparisons, "totals") + attr(cmp.2to3$comparisons, "totals")
  match.count <- colSums(cmp.1to3$comparisons[matchrows.1to3,,drop=FALSE]) + colSums(cmp.2to3$comparisons[matchrows.2to3,,drop=FALSE])
  nonmatch.count <- total.count - match.count # (no need to remove non-candidiates when link tracing)
  # 4. Calculate and return log likelihood based on count of field levels for
  #    matches and non-matches
  loglkl <- (
    sum(match.count * log(m))
    + sum(nonmatch.count * log(u))
  )
  return(loglkl)
}

# Function to calculate the log of the prior of Z2 given Z and the hyperparameters
# aBM and bBM. If the value of Z2 is impossible given the value of Z (i.e.
# due to linking to a non-candidate) this returns -Inf := log(0)
# In R, negative infinity plus any finite value equals negative infinity, and all
# finite values compare greater than negative infinity. So this value can be used
# in acceptance ratio computations without worrying.
calc.log.Z2prior <- function(n1, n2, n3, Z2, Z, aBM, bBM) {
  # Calculate the log of the prior and return
  if (valid.link.state(n1, Z, Z2)) {
    # Candidates are any unlinked entries in file 1, plus all entries in file 2
    cand <- c(setdiff(seq_len(n1), Z), n1 + seq_len(n2))
    # How many links are there currently between file 3 and the candidate set?
    nlinked <- sum(Z2 <= n1 + n2)
    # Work in log scale: replacing log(x!) with lgamma(x+1) and log(beta(a,b))
    # with lbeta(a,b)
    logprior <- (
      (lgamma(length(cand) - nlinked + 1) - lgamma(length(cand) + 1))
      + (lbeta(aBM + nlinked, bBM + n3 - nlinked) - lbeta(aBM, bBM))
    )
    return(logprior)
  } else {
    # Impossible -> probability zero -> log probability -Inf.
    return(-Inf)
  }
}

################################################################################
# Helper Functions for Likelihood/Prior functions ##############################
################################################################################

# Determines whether the given values of Z, Z2 are viable according to the
# candidate set rules. I.e. returns true if no records are linked to the same
# record in file 1. Else returns false.
valid.link.state <- function(n1, Z, Z2) {
  noncand <- Z[Z <= n1] # Non-candidates: records in file 1 with links
  return(!any(Z2 %in% noncand))
}

# Determines which rows in comparison data correspond to matches, according to
# the match vector supplied
# Parameters:
#   cmpdata = comparison data in the format returned by BRL::compareRecords
#   Z2 = a vector of matches the length of "file 2" in the comparison data. Each
#        entry is for an element of "file 2"
#   offset = an integer equal to the total number of records in files "earlier"
#            than "file 1" in the comparison data. This offset is essentially
#            subtracted from the values of Z2 to rule out matches with any of
#            these earlier files.
# Example usage:
#   In a streaming setting, cmpdata compares file 3 to file 5. The file sizes are
#   n1, n2, n3, n4, and n5 respectively. Z2 is a vector of length n5 with values
#   from 1 to n1+n2+n3+n4+n5. Values between n1+n2+1 and n1+n2+n3 represent the
#   links with file 3. Values larger than n1+n2+n3+n4 represent unlinked records
#   You would call
#     matchrows(cmpdata, Z2, offset=n1+n2)
#   to return the rows in cmpdata the correspond to matches between file 3 and
#   file 5.
# Notes:
#   Link tracing is not done. Link tracing must be pre-incorporated into Z2.
matchrows <- function(cmpdata, Z2, offset=0) {
  n1 <- cmpdata$n1
  n2 <- cmpdata$n2
  # Check that we have a correct Z2
  if (length(Z2) != n2) {
    stop("matchrows() - length of Z2 does not match cmpdata's file 2")
  }
  # Boolean index of records in "file 2" which are linked to "file 1"
  linked <- (Z2 > offset) & (Z2 <= offset+n1)
  # For links, return j*n1 + i (subtracting offset from i)
  return((which(linked) - 1) * n1 + (Z2[linked] - offset))
}

# Determines which rows in comparison data correspond to noncandidates, according to
# the match vector supplied
# Parameters:
#   cmpdata = comparison data in the format returned by BRL::compareRecords. The
#             "file 1" in this structure corresponds to a file number m, and
#             "file 2" corresponds to a file number k > m
#   Z = a vector of matches. Its length is unimportant, but it must obey this
#       constraint: An index i appears in Z, corresponding to a record in
#       cmpdata's "file 1", IF AND ONLY IF it is linked to by a record in some
#       file numbered l, m < l < k. I.e. if an index i appears in Z then that
#       record is not a candidate.
#   offset = an integer equal to the total number of records in files "earlier"
#            than "file 1" in the comparison data.
# Example usage:
#   In a streaming setting, cmpdata compares file 2 to file 5. The file sizes are
#   n1, n2, n3, n4, and n5 respectively. Z is a vector of length n3+n4 with values
#   from 1 to n1+n2+n3+n4. Values between n1+1 and n1+n2 represent the
#   links with file 2.
#   You would call
#     noncandrows(cmpdata, Z, offset=n1)
#   to return the rows in cmpdata that correspond to matches involving a record
#   in file 2 which are not candidates.
noncandrows <- function(cmpdata, Z=c(), offset=0) {
  n1 <- cmpdata$n1
  n2 <- cmpdata$n2
  # Indices of records in "file 1" which are linked to by the vector Z
  noncand.records <- Z[(Z > offset) & (Z <= offset+n1)]
  # Non-candidate pairs are any pair involving any of these records
  return(c(outer((seq_len(n2) - 1) * n1, noncand.records, "+")))
}

# Function to return a link-traced version of the link vector Z2, following down
# the vector Z for the given number of steps
# Parameters:
#   Z2 = The vector to be traced
#   Z = The vector containing the next steps to be followed
#   steps = The number of additional steps to take
#   offset = An offset of the starting position of the indices of Z. Any records
#            with indices less than or equal to offset are assumed to be the end
#            of their chain. This can be used to exclude file 1, for example, if
#            its links are excluded from Z to save space
# Returns:
#   Z2, but with links traces 'steps' further steps according to the links
#   defined in Z.
trace <- function(Z2, Z=Z2, steps=0, offset=0) {
  # Vector which will be updated with traces
  Z2.traced <- Z2
  n.prev.records <- length(Z) + offset
  for (step in 1:steps) {
    # At this step, which records are linked to previous records which can be
    # further traced?
    traceable <- (Z2.traced > offset) & (Z2.traced <= n.prev.records)
    # Follow one more step with Z
    Z2.traced[traceable] <- Z[Z2.traced[traceable] - offset]
  }
  Z2.traced
}
