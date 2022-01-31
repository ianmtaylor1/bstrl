

# Calculate the log prior of links
# sl = streaming link object
# aBM, bBM = prior parameters
# vec = whether to calculate the prior for all Z vectors or just the last one,
#       log(p(Z1)p(Z2)...p(Zk)) or log(p(Zk))
log.Zprior <- function(sl, aBM, bBM, vec = c("all", "last")) {
  vec <- match.arg(vec)
  # Global indexes for the start and end of each file
  fileend <- cumsum(sl$ns)
  filestart <- fileend - sl$ns + 1
  # Which file to start and end with
  logp <- 0
  first <- if (vec == "all") { 2 } else { length(sl$ns) }
  last <- length(sl$ns)
  # Calculate the log prior
  for (f in seq(first, last)) {
    Zk <- sl$Z[seq(filestart[f], fileend[f])]
    N <- sl$ns[f]
    nprev <- fileend[f-1]
    nlinked <- sum(Zk < seq(filestart[f], fileend[f]))
    logp <- logp + (
      (lgamma(nprev - nlinked + 1) - lgamma(nprev + 1))
      + (lbeta(aBM + nlinked, bBM + N - nlinked) - lbeta(aBM, bBM))
    )
  }

  return(logp)
}
