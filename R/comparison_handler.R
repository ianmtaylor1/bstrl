

# Tally up disagreement counts for comparisons between a file k and all k-1
# previous files, WHERE file k is the LAST FILE linked by sl.
# This function is really identical to disag.counts.onefile but having it
# separate makes logic easier to understand.
# Parameters:
#   lastfilecmps - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#   sl = streaminglinks object defining current link state Must represent links
#        between exactly k files
# Returns:
#   the total counts of disagreement levels for matches and nonmatches as a list
#   with elements match and nonmatch
disag.counts.lastfile <- function(lastfilecmps, sl) {

  # Enforce that these comparisons must be for the last file in sl. This is the
  # only real contribution of having this as a separate function, but it makes
  # logic elsewhere in the code more legible.
  stopifnot(length(lastfilecmps) == nfiles(sl) - 1)

  disag.counts.onefile(lastfilecmps, sl)
}

# Tally up disagreement counts for comparisons between a file k and all k-1
# previous files.
# Parameters:
#   filecmps - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#   sl = streaminglinks object defining current link state Must represent links
#        between at least k files
# Returns:
#   the total counts of disagreement levels for matches and nonmatches as a list
#   with elements match and nonmatch
disag.counts.onefile <- function(filecmps, sl) {

  # First, extract list of all linked pairs of records in local indexing to
  # separate them by file.
  links <- alllinks(sl, idx="local")

  # Find the number of the latest file to compare
  nfiles <- length(filecmps) + 1

  tot.match.count <- tot.nonmatch.count <- rep(0, sum(filecmps[[1]]$nDisagLevs))

  # Go through each comparison data object file by file
  for (f in seq_len(nfiles - 1)) {
    pairidx <- (links$file1 == f) & (links$file2 == nfiles)
    n1 <- filecmps[[f]]$n1

    if (sum(pairidx) > 0) {
      leftrecords <- links$record1[pairidx]
      rightrecords <- links$record2[pairidx]
      rows <- (rightrecords - 1) * n1 + leftrecords
      match.count <- colSums(filecmps[[f]]$comparisons[rows,,drop=F])
    } else {
      match.count <- rep(0, sum(filecmps[[f]]$nDisagLevs))
    }
    nonmatch.count <- attr(filecmps[[f]]$comparisons, "totals") - match.count

    # Add this file's contribution to total
    tot.match.count <- tot.match.count + match.count
    tot.nonmatch.count <- tot.nonmatch.count + nonmatch.count
  }

  return(list(match=tot.match.count, nonmatch=tot.nonmatch.count))
}


# Tally up disagreement counts for comparisons involving a file k or any later
# files. i.e. first find counts for file k and previous files, then file k+1 and
# previous files, and so on until the last file and previous files.
# Parameters:
#   cmpdata - a list of lists of comparison data objects. (Standard triangular
#             comparison data format.)
#   sl = streaminglinks object defining current link state
# Returns:
#   the total counts of disagreement levels for matches and nonmatches as a list
#   with elements match and nonmatch
disag.counts.tailfiles <- function(cmpdata, sl, startfile) {
  stopifnot(length(cmpdata) == nfiles(sl) - 1)
  stopifnot(startfile >= 2, startfile <= nfiles(sl))

  tot.match.count <- tot.nonmatch.count <- 0

  for (f in seq(startfile, nfiles(sl))) {
    tmp <- disag.counts.onefile(cmpdata[[f - 1]], sl)
    tot.match.count <- tot.match.count + tmp$match
    tot.nonmatch.count <- tot.nonmatch.count + tmp$nonmatch
  }

  return(list(match=tot.match.count, nonmatch=tot.nonmatch.count))
}


disag.counts.allfiles <- function(cmpdata, sl) {
  disag.counts.tailfiles(cmpdata, sl, 2)
}
