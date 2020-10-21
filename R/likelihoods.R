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
calc.log.lkl <- function(cmpdata, m, u, Z, Z2) {
  # How many files are we dealing with?
  nfiles <- length(cmpdata) + 1
  # Which records in files 1 through nfiles-1 are not candidates? i.e. are matched
  # to by a later file already according to Z
  n1 <- cmpdata[[1]]$n1
  noncand <- Z[Z <= n1 + seq_len(length(Z))]
  # Loop through each previous file and accumulate log likelihood
  loglkl <- 0
  totalprevrecords <- 0
  for (file in 1:(nfiles-1)) {
    # Which rows in this comparison object correspond to matches?
    filematchrows <- matchrows(cmpdata[[file]], Z2, offset=totalprevrecords)
    match.count <- colSums(cmpdata[[file]]$comparisons[filematchrows,,drop=FALSE])
    # Which rows in this comparison object correspond to non-candidates?
    filenoncandrows <- noncandrows(cmpdata[[file]], Z=noncand, offset=totalprevrecords)
    noncand.count <- colSums(cmpdata[[file]]$comparisons[filenoncandrows,,drop=FALSE])
    # Which rows does that leave for nonmatches?
    nonmatch.count <- attr(cmpdata[[file]]$comparisons, "totals") - match.count - noncand.count
    # Update the log likelihood
    loglkl <- loglkl + sum(match.count * log(m)) + sum(nonmatch.count + log(u))
    # Update the number of records in previous files
    totalprevrecords <- totalprevrecords + cmpdata[[file]]$n1
  }
  # After going through all files, return the total log likelihood
  return(loglkl)
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
  # How many files are we dealing with?
  nfiles <- length(cmpdata) + 1
  # How many records are in the first file? (needed later for tracing)
  n1 <- cmpdata[[1]]$n1
  # Loop through each previous file and accumulate log likelihood
  loglkl <- 0
  totalprevrecords <- 0
  for (file in 1:(nfiles-1)) {
    # Which rows in this comparison object correspond to matches?
    # Need to trace one less than the file difference between the two compared files
    Z2.traced <- trace(Z2, Z, steps=nfiles - file - 1, offset=n1)
    filematchrows <- matchrows(cmpdata[[file]], Z2.traced, offset=totalprevrecords)
    match.count <- colSums(cmpdata[[file]]$comparisons[filematchrows,,drop=FALSE])
    # Which rows does that leave for nonmatches?
    nonmatch.count <- attr(cmpdata[[file]]$comparisons, "totals") - match.count
    # Update the log likelihood
    loglkl <- loglkl + sum(match.count * log(m)) + sum(nonmatch.count + log(u))
    # Update the number of records in previous files
    totalprevrecords <- totalprevrecords + cmpdata[[file]]$n1
  }
  # After going through all files, return the total log likelihood
  return(loglkl)
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
    ncand <- n1 + nprev - sum(Z <= n1 + seq_len(nprev))
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
  noncand <- Z[Z <= offset + seq_len(length(Z))] # Non-candidates: records with links in later files
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
