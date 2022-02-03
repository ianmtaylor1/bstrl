

# Calculate the log likelihood of comparisons only involving the last file
# Parameters:
#   lastfilecmps = list of k-1 comparison data objects, comparing each of the first
#     k-1 files with the kth file. This is the last element in the triangular
#     list-of-lists format for all comparison data.
#   m,u = current parameter values
#   sl = streaminglinks object defining current link state
# Returns:
#   a single number, the log likelihood of Gamma^{(k)}, all comparisons
#   involving the most recent file k.
calc.log.lkl.lastfile <- function(lastfilecmps, m, u, sl) {
  # Calculate the disagreement level counts
  counts <- disag.counts.lastfile(lastfilecmps, sl)
  # The log likelihood is these counts combined with log m and log u
  return(sum(counts$match * log(m)) + sum(counts$nonmatch * log(u)))
}

calc.log.lkl.tailfiles <- function(cmpdata, m, u, sl, startfile) {
  # Calculate the disagreement level counts
  counts <- disag.counts.tailfiles(cmpdata, sl, startfile)
  # The log likelihood is these counts combined with log m and log u
  return(sum(counts$match * log(m)) + sum(counts$nonmatch * log(u)))
}

calc.log.lkl.allfiles <- function(cmpdata, m, u, sl) {
  # Calculate the disagreement level counts
  counts <- disag.counts.allfiles(cmpdata, sl)
  # The log likelihood is these counts combined with log m and log u
  return(sum(counts$match * log(m)) + sum(counts$nonmatch * log(u)))
}
