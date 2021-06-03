# This file contains functions useful for handling comparison data objects (or
# lists of comparison data objects). Mainly for extracting disagreement level
# counts for matches and nonmatches based on current values of Z2 and Z.
# These functions are referenced in likelihood evaluation, m/u full conditional
# sampling, and postprocessing of samples for m/u full conditional comparison
# summaries.


# Given a list of comparison data objects and a link state (Z and Z2), returns
# the total counts of disagreement levels for matches and nonmatches
# Parameters:
#   cmpdata - "list of lists" comparison data format
#       Z,Z2 - parameter values. Z2 is the links originating from file k, and Z
#              is the links originating from files 2 through k-1, concatenated
#              together (file 1 is excluded because it's unnecessary under the
#              assumption of no duplicates within files.)
# Returns:
#   A list with two components, match and nonmatch, which are the total counts
#   of disagreement levels for matches and nonmatches respectfully. Each vector
#   is the length of the number of columns of a comparison matrix, which is the
#   total number of disagreement levels of all fields that are compared.
disag.counts.smcmc.fast <- function(cmpdata, Z, Z2) {
  # How many files are we dealing with?
  nfiles <- length(cmpdata) + 1
  # Size of files: file 1 is special, and let's keep the rest in a vector
  n1 <- cmpdata[[1]][[1]]$n1
  ns <- rep(0, nfiles)
  for (i in seq_len(nfiles - 1)) {
    ns[i] <- cmpdata[[nfiles - 1]][[i]]$n1
  }
  ns[nfiles] <- cmpdata[[nfiles - 1]][[nfiles - 1]]$n2
  # Total disagreement levels
  nlevels <- sum(cmpdata[[1]][[1]]$nDisagLevs)

  # Concatenate Z and Z2 into a single vector we can repeatedly trace
  Z.all <- c(Z, Z2)

  # Initial match counts:
  match.count <- rep(0, nlevels)
  nonmatch.count <- rep(0, nlevels)

  # Match counts: one pair at a time, tracing

  # The initial "trace" is just from no links to the base form of Z/Z2
  Z.prev <- n1 + seq_len(sum(ns) - n1)
  Z.traced <- Z.all
  # Take the maximum number of steps
  for (steps in seq_len(nfiles - 1)) {
    # Which links did tracing change?
    different <- which(Z.traced != Z.prev)

    if (length(different) == 0) break # If nothing has changed we don't have to keep looping

    # Convert each end of the link to file,idx pair format
    fileidx.from <- infile(n1 + different, ns)
    fileidx.to <- infile(Z.traced[different], ns)

    # For each new link, tally its disagreement level counts
    for (i in seq_len(different)) {
      match.count <- match.count +
        comparison(cmpdata, fileidx.to$file[i], fileidx.from$file[i],
                   fileidx.to$idx[i], fileidx.from$idx[i])
    }

    # Trace 1 more step
    Z.prev <- Z.traced
    Z.traced <- trace(Z.traced, Z.all, steps=1)
  }

  # Nonmatch counts: totals minus match counts
  for (f1 in seq_len(nfiles - 1)) {
    for (f2 in seq_len(f1)) {
      nonmatch.count <- nonmatch.count + attr(cmpdata[[f1]][[f2]]$comparisons, "totals")
    }
  }
  nonmatch.count <- nonmatch.count - match.count

  # return both
  return(list(match=match.count, nonmatch=nonmatch.count))
}

# Given a vector of record indices (globalized) and a vector of file sizes,
# tell which file each record belongs to
# Returns a list of two vectors: file and idx, each the same length as records.
# file has indices of the file a record belongs to (1 to length(ns))
# idx has the corresponding index within the file
# e.g.
# records = c(5, 12, 27)
# ns = c(6, 3, 10, 10)
# Return:
#   list(file=c(1, 3, 4), idx=c(5, 3, 8))
infile <- function(records, ns) {
  # Thresholds of the cumulative number of records through file k
  thresh <- cumsum(ns)

  # Find the file number of each record
  fileno <- findInterval(records, thresh, left.open=T) + 1

  # Find the index within each file
  idxno <- records - (thresh - ns)[fileno]

  return(list(file=fileno, idx=idxno))
}


# Given a list of comparison data objects and a link state (Z and Z2), returns
# the total counts of disagreement levels for matches and nonmatches
# Parameters:
#   cmpdata - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#       Z,Z2 - parameter values. Z2 is the links originating from file k, and Z
#              is the links originating from files 2 through k-1, concatenated
#              together (file 1 is excluded because it's unnecessary under the
#              assumption of no duplicates within files.)
#   do.trace - Logical, whether to do link tracing.
# Returns:
#   A list with two components, match and nonmatch, which are the total counts
#   of disagreement levels for matches and nonmatches respectfully. Each vector
#   is the length of the number of columns of a comparison matrix, which is the
#   total number of disagreement levels of all fields that are compared.
disag.counts <- function(cmpdata, Z, Z2, do.trace) {
  # How many files are we dealing with?
  nfiles <- length(cmpdata) + 1
  # Which records in files 1 through nfiles-1 are not candidates? i.e. are matched
  # to by a later file already according to Z
  n1 <- cmpdata[[1]]$n1
  noncand <- Z[Z < n1 + seq_len(length(Z))]
  # Loop through each previous file and accumulate log likelihood
  tot.match.count <- tot.nonmatch.count <- rep(0, sum(cmpdata[[1]]$nDisagLevs))
  totalprevrecords <- 0
  for (file in seq_len(nfiles-1)) {
    # Matches and non-candidates change depending on whether tracing is done
    if (do.trace) {
      # Which rows in this comparison object correspond to matches?
      # Need to trace one less than the file difference between the two compared files
      Z2.traced <- trace(Z2, Z, steps=nfiles - file - 1, offset=n1)
      filematchrows <- matchrows(cmpdata[[file]], Z2.traced, offset=totalprevrecords)
      file.match.count <- colSums(cmpdata[[file]]$comparisons[filematchrows,,drop=FALSE])
      # There are no noncandidates with tracing, therefore sums are zero
      file.noncand.count <- rep(0, sum(cmpdata[[file]]$nDisagLevs))
    } else {
      # Which rows in this comparison object correspond to matches?
      filematchrows <- matchrows(cmpdata[[file]], Z2, offset=totalprevrecords)
      file.match.count <- colSums(cmpdata[[file]]$comparisons[filematchrows,,drop=FALSE])
      # Which rows in this comparison object correspond to non-candidates?
      filenoncandrows <- noncandrows(cmpdata[[file]], Z=noncand, offset=totalprevrecords)
      file.noncand.count <- colSums(cmpdata[[file]]$comparisons[filenoncandrows,,drop=FALSE])
    }
    # Which rows does that leave for nonmatches?
    file.nonmatch.count <- attr(cmpdata[[file]]$comparisons, "totals") - file.match.count - file.noncand.count
    # Update the total counts so far
    tot.match.count <- tot.match.count + file.match.count
    tot.nonmatch.count <- tot.nonmatch.count + file.nonmatch.count
    # Update the number of records in previous files
    totalprevrecords <- totalprevrecords + cmpdata[[file]]$n1
  }
  # Package the counts into one list and return
  return(list(match=tot.match.count, nonmatch=tot.nonmatch.count))
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
  for (step in seq_len(steps)) {
    # At this step, which records are linked to previous records which can be
    # further traced?
    traceable <- (Z2.traced > offset) & (Z2.traced <= n.prev.records)
    # Follow one more step with Z
    Z2.traced[traceable] <- Z[Z2.traced[traceable] - offset]
  }
  Z2.traced
}


# Function to return a single comparison between any two records in any two files
# Parameters:
#   cmpdata = a "triangular" list of lists. first list contains comparisons with file 2
#             second contains comparisons with file 3 in order, etc.
#   file1 = index of the first file
#   file2 = index of the second file. file2 > file1
#   rec1 = index of the first record. 1 <= rec1 <= n_{file1}
#   rec2 = index of the second record. 1 <= rec2 <= n_{file2}
comparison <- function(cmpdata, file1, file2, rec1, rec2) {
  stopifnot(file1 < file2)

  filepair <- cmpdata[[file2 - 1]][[file1]]
  nf1 <- filepair$n1
  nf2 <- filepair$n2

  stopifnot(rec1 >= 1, rec1 <= nf1, rec2 >= 1, rec2 <= nf2)
  stopifnot(nrow(filepair$comparisons) == nf1 * nf2)

  filepair$comparisons[(rec2 - 1) * nf1 + rec1, ]
}
