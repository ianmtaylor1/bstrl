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

