# Functions in this perform bipartite linkage, either by wrapping the BRL package
# or by doing its own sampling


# Function to do bipartite bayesian RL - all but the point estimate part.
# Returns direct Gibbs samples without post-processing or finding point
# estimates (e.g. for PPRB).
# Copy/paste from the BRL source code, just removing the call to linkRecords()
# Therefore, arguments are the same (minus the loss weight arguments for the
# point estimates)
# Returns a list with five items:
#   Z - a matrix with n2 rows and nIter-burn columns
#   m - a matrix with \sum_{f=1}^F L_f rows and nIter-burn columns
#   u - a matrix with \sum_{f=1}^F L_f rows and nIter-burn columns
#   m.fc.pars - Used for further sampling. Contains the values that, when added
#     to the prior parameters for m, allow you to sample from the conjugate full
#     conditional of m given Z and the comparisons.
#   u.fc.pars - Same as ^ but for u.
#' @export
bipartiteRL.precmp <- function(cmpdata, nIter=1000, burn=round(nIter*.1), a=1, b=1, aBM=1, bBM=1, seed=0) {

  # 'burn' is not a parameter of the other functions, so we control it here
  if( !is.numeric(burn) | (burn<0) | (burn>=nIter) )
    stop("burn should be an integer that satisfies 0 <= burn < nIter")

  # Gibbs sampling from posterior of bipartite matchings
  chain <- BRL::bipartiteGibbs(cmpdata, nIter, a, b, aBM, bBM, seed)

  # Calculate comparison summaries for each value of Z
  total.counts <- colSums(cmpdata$comparisons)
  m.fc.pars <- matrix(0, nrow=nrow(chain$m), ncol=nIter)
  u.fc.pars <- matrix(0, nrow=nrow(chain$m), ncol=nIter)
  for (i in 1:nIter) {
    match.idx <- matchrows(cmpdata, Z)
    match.counts <- colSums(cmpdata$comparisons[match.idx,,drop=FALSE])
    m.fc.pars[,i] <- match.counts
    u.fc.pars[,i] <- total.counts - match.counts
  }

  # Filter the burn-in iterations
  iterfilter <- setdiff(1:nIter,seq_len(burn))
  list(Z=chain$Z[,iterfilter,drop=FALSE],
       m=chain$m[,iterfilter,drop=FALSE],
       u=chain$u[,iterfilter,drop=FALSE],
       m.fc.pars=m.fc.pars[,iterfilter,drop=FALSE],
       u.fc.pars=u.fc.pars[,iterfilter,drop=FALSE])
}

#' @export
bipartiteRL <- function(df1, df2, flds=NULL, flds1=NULL, flds2=NULL, types=NULL, breaks=c(0,.25,.5),
                      nIter=1000, burn=round(nIter*.1), a=1, b=1, aBM=1, bBM=1, seed=0){

  # create comparison data
  myCompData <- BRL::compareRecords(df1, df2, flds, flds1, flds2, types, breaks)

  # Call the function for precompared records
  BRL.gibbs.precmp(myCompData, nIter, burn, a, b, aBM, bBM, seed)
}
