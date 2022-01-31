

# Calculate the log likelihood of comparisons only involving the last file
# Parameters:
#   cmpdata = list of k-1 comparison data objects, comparing each of the first
#     k-1 files with the kth file. This is the last element in the triangular
#     list-of-lists format for all comparison data.
#   m,u = current parameter values
#   sl = streaminglinks object defining current link state
# Returns:
#   a single number, the log likelihood of Gamma^{(k)}, all comparisons
#   involving the most recent file k.
calc.log.lkl.lastfile <- function(cmpdata, m, u, sl) {
  # Calculate the disagreement level counts
  counts <- disag.counts.lastfile(cmpdata, sl)
  # The log likelihood is these counts combined with log m and log u
  return(sum(counts$match * log(m)) + sum(counts$nonmatch * log(u)))
}
