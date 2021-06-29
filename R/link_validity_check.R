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
  Z2linked <- Z2[Z2 <= offset + length(Z)]
  if (any(Z2linked %in% noncand)) {
    # Invalid state #1: Z2 links to records with links in Z1
    return(FALSE)
  } else if (!check.no.duplicates(noncand)) {
    # Invalid state #2: Z1 contains duplicate links
    return(FALSE)
  } else if (!check.no.duplicates(Z2linked)) {
    # Invalid state #3: Z2 contains duplicate links
    return(FALSE)
  }
  return(TRUE)
}

# Check that a vetor Z contains no duplicate values. This is a minimum validity
# check, that checks whether a link parameter vector Z (Z1, Z2, etc) doesn't
# have any duplicates in itself.
check.no.duplicates <- function(Z) {
  if (length(unique(Z)) != length(Z)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}
