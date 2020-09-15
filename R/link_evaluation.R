# This file contains functions for evaluating the links output by the sampler

# Given an estimated Z,Z2, and a true Z,Z2: determine the precision of the
# estimates.
#' @export
precision <- function(n1, Z, Z2, true.Z, true.Z2) {
  n2 <- length(Z)
  n3 <- length(Z2)
  # First, create a copy of Z2 (and true.Z2) except following traced links where applicable.
  # I.e., any links to file 2 should be replaced with a record in file 1, if the linked
  # records in file 2 are further linked to file 1. Take advantage of the fact that
  # unlinked records essentially point to themselves.
  Z2.traced <- Z2
  true.Z2.traced <- true.Z2
  Z2.traced[which((Z2 > n1) & (Z2 <= n1+n2))] <- Z[Z2[which((Z2 > n1) & (Z2 <= n1+n2))] - n1]
  true.Z2.traced[which((true.Z2 > n1) & (true.Z2 <= n1+n2))] <- true.Z[true.Z2[which((true.Z2 > n1) & (true.Z2 <= n1+n2))] - n1]
  # Next, calculate the number of true positive, and total estimated links
  tot.links <- sum(Z2 <= n1+n2) + sum(Z <= n1) + sum((Z2.traced != Z2)&(Z2.traced <= n1))
  truepos.links <- (
    # True positives between files 1 and 2
    sum((Z <= n1)&(Z == true.Z)) +
    # True positives between files 2 and 3
    sum((Z2 > n1)&(Z2 <= n1+n2)&(Z2 == true.Z2)) +
    # True positives between files 1 and 3
    sum((Z2.traced <= n1)&(Z2.traced == true.Z2.traced))
  )
  # Return precision
  return(truepos.links/tot.links)
}

# Given an estimated Z,Z2, and a true Z,Z2: determine the recall of the
# estimates.
#' @export
recall <- function(n1, Z, Z2, true.Z, true.Z2) {
  n2 <- length(Z)
  n3 <- length(Z2)
  # First, create a copy of Z2 (and true.Z2) except following traced links where applicable.
  # I.e., any links to file 2 should be replaced with a record in file 1, if the linked
  # records in file 2 are further linked to file 1. Take advantage of the fact that
  # unlinked records essentially point to themselves.
  Z2.traced <- Z2
  true.Z2.traced <- true.Z2
  Z2.traced[which((Z2 > n1) & (Z2 <= n1+n2))] <- Z[Z2[which((Z2 > n1) & (Z2 <= n1+n2))] - n1]
  true.Z2.traced[which((true.Z2 > n1) & (true.Z2 <= n1+n2))] <- true.Z[true.Z2[which((true.Z2 > n1) & (true.Z2 <= n1+n2))] - n1]
  # Next, calculate the number of true positive, and total estimated links
  true.tot.links <- sum(true.Z2 <= n1+n2) + sum(true.Z <= n1) + sum((true.Z2.traced != true.Z2)&(true.Z2.traced <= n1))
  truepos.links <- (
    # True positives between files 1 and 2
    sum((Z <= n1)&(Z == true.Z)) +
    # True positives between files 2 and 3
    sum((Z2 > n1)&(Z2 <= n1+n2)&(Z2 == true.Z2)) +
    # True positives between files 1 and 3
    sum((Z2.traced <= n1)&(Z2.traced == true.Z2.traced))
  )
  # Return recall
  return(truepos.links/true.tot.links)
}

