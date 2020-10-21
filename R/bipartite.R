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
bipartiteRL.precmp <- function(cmpdata, nIter=1000, burn=round(nIter*.1), a=1, b=1, aBM=1, bBM=1,
                               seed=0, method="BRL", blocksize=NULL) {

  # 'burn' is not a parameter of the other functions, so we control it here
  if( !is.numeric(burn) | (burn<0) | (burn>=nIter) )
    stop("burn should be an integer that satisfies 0 <= burn < nIter")

  # Gibbs sampling from posterior of bipartite matchings
  if (method == "BRL") {
    chain <- BRL::bipartiteGibbs(cmpdata, nIter, a, b, aBM, bBM, seed)
  } else if (method == "LB") {
    # Draw initial values for m and u from the prior
    if (length(a) == 1) a <- rep(a, sum(cmpdata$nDisagLevs))
    if (length(b) == 1) b <- rep(b, sum(cmpdata$nDisagLevs))
    m.curr <- rdirichlet.multi(a, cmpdata$nDisagLevs)
    u.curr <- rdirichlet.multi(b, cmpdata$nDisagLevs)
    # Initial Z is all unlinked
    Z.curr <- cmpdata$n1 + seq_len(cmpdata$n2)
    # Initialize arrays to hold samples
    m.samples <- u.samples <- matrix(0, nrow=sum(cmpdata$nDisagLevs), ncol=nIter)
    Z.samples <- matrix(0, nrow=cmpdata$n2, ncol=nIter)
    # Iterate through sampler
    for (i in 1:nIter) {
      # 1. Full conditional update of Z conditioned on m and u.
      tmp <- draw.Z2.informed(list(cmpdata), c(), Z.curr, m, u, aBM, bBM, blocksize=blocksize)
      Z.prop <- tmp$Z2
      mhmod <- tmp$mod
      log.alpha <- (
        calc.log.lkl(list(cmpdata), m.curr, u.curr, c(), Z.prop)
        - calc.log.lkl(list(cmpdata), m.curr, u.curr, c(), Z.curr)
        + calc.log.Z2prior(cmpdata$n1, Z.prop, c(), aBM, bBM)
        - calc.log.Z2prior(cmpdata$n1, Z.curr, c(), aBM, bBM)
        + log(mhmod)
      )
      if (-rexp(1) < log.alpha) {
        Z.curr <- Z.prop
      }
      Z.samples[,i] <- Z.curr
      # 2. Full conditional update of m and u conditioned on Z.
      tmp <- r_m_u_fc(list(cmpdata), c(), Z.curr, a, b, 0, 0)
      m.curr <- m.samples[,i] <- tmp$m
      u.curr <- u.samples[,i] <- tmp$u
    }
    # Put into same format as BRL
    chain <- list(m=m.samples, u=u.samples, Z=Z.samples)
  } else {
    stop(paste0("Method '", method, "' unknown."))
  }

  # Calculate comparison summaries for each value of Z
  # At bipartite stage, there is only one cmpdata object and tracing is irrelevant
  total.counts <- colSums(cmpdata$comparisons)
  m.fc.pars <- matrix(0, nrow=nrow(chain$m), ncol=nIter)
  u.fc.pars <- matrix(0, nrow=nrow(chain$m), ncol=nIter)
  for (i in 1:nIter) {
    tmp <- disag.counts(list(cmpdata), Z=c(), Z2=chain$Z[,i], do.trace=FALSE)
    m.fc.pars[,i] <- tmp$match
    u.fc.pars[,i] <- tmp$nonmatch
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
