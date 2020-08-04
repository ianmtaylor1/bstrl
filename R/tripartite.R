# This file contains the main functions to do tripartite (or streaming) linkage

# Tripartite RL function built on the BRL package. Arguments the same as
# BRL.gibbs, but with the additional df3 and flds3 arguments.
# DESIGN CHOICES
#  - Discard burn in BRL.gibbs initially, or pass burn=0 and discard burn
#    at end of this function?
#' @export
tripartiteRL.precmp <- function(cmpdata.1to2, cmpdata.1to3, cmpdata.2to3, trace=FALSE,
                                nIter=1000, burn=round(nIter*.1), a=1, b=1, aBM=1, bBM=1, seed=0) {
  # 1. Size of files
  n1 <- cmpdata.1to3$n1
  n2 <- cmpdata.2to3$n1
  n3 <- cmpdata.1to3$n2
  # 2. Do bipartite RL between df1 and df2 using BRL.gibbs()
  #    Do not burn anything now! Keep it all, and burn later
  bipartite.samp <- BRL.gibbs.precmp(cmpdata.1to2, nIter, 0, a, b, aBM, bBM, seed)
  # 3. Create precomputed data structures
  comparisons.1to3 <- preproc.cmpdata(cmpdata.1to3)
  comparisons.2to3 <- preproc.cmpdata(cmpdata.2to3)
  # 4. Set up empty arrays of the appropriate size that will eventually be returned
  m.samples  <- matrix(0, nrow=nrow(bipartite.samp$m), ncol=nIter)
  u.samples  <- matrix(0, nrow=nrow(bipartite.samp$u), ncol=nIter)
  Z.samples  <- matrix(0, nrow=nrow(bipartite.samp$Z), ncol=nIter)
  Z2.samples <- matrix(0, nrow=n3,                     ncol=nIter)
  # 5. Initialize the first sample (m,u,Z: first sample from bipartite. Z2: unlinked)
  m.samples[,1] <- bipartite.samp$m[,1]
  u.samples[,1] <- bipartite.samp$u[,1]
  Z.samples[,1] <- bipartite.samp$Z[,1]
  Z2.samples[,1] <- n1 + n2 + seq_len(n3)
  # Create array to hold record of acceptances: first row = m,u,Z. second row = Z2
  # Contains all 1's by default, zeros written on rejection
  accepted <- matrix(1, nrow=2, ncol=nIter)
  # Choose the likelihood function based on tracing flag
  if (trace) {
    ell <- calc.log.lkl.tracing
  } else {
    ell <- calc.log.lkl
  }
  # 6. Perform M-H for each iteration
  for (i in seq(2, nIter)) {
    # What are the current values of the parameters?
    m.curr <- m.samples[,i-1]
    u.curr <- u.samples[,i-1]
    Z.curr <- Z.samples[,i-1]
    Z2.curr <- Z2.samples[,i-1]
    # 6.1: Propose new m,u,Z collectively. Take sequential samples from bipartite
    # TODO Future: resample with replacement from bipartite posterior samples?
    m.prop <- bipartite.samp$m[,i]
    u.prop <- bipartite.samp$u[,i]
    Z.prop <- bipartite.samp$Z[,i]
    # Decide whether to accept new values
    log.alpha1 <- (
      ell(comparisons.1to3, comparisons.2to3, n1, n2, n3, m.prop, u.prop, Z.prop, Z2.curr)
      - ell(comparisons.1to3, comparisons.2to3, n1, n2, n3, m.curr, u.curr, Z.curr, Z2.curr)
      + calc.log.Z2prior(n1, n2, n3, Z2.curr, Z.prop, aBM, bBM)
      - calc.log.Z2prior(n1, n2, n3, Z2.curr, Z.curr, aBM, bBM)
    )
    if (-rexp(1) < log.alpha1) {
      # accept
      m.samples[,i] <- m.prop
      u.samples[,i] <- u.prop
      Z.samples[,i] <- Z.prop
    } else {
      # reject
      m.samples[,i] <- m.curr
      u.samples[,i] <- u.curr
      Z.samples[,i] <- Z.curr
      accepted[1,i] <- 0
    }
    # Reset "current" values of parameters
    m.curr <- m.samples[,i]
    u.curr <- u.samples[,i]
    Z.curr <- Z.samples[,i]
    # 6.2: Propose new Z2
    tmp <- draw.Z2.informed(n1, n2, n3, Z.curr, Z2.curr, m.curr, u.curr,
                            comparisons.1to3, comparisons.2to3,
                            aBM, bBM, trace=trace)
    Z2.prop <- tmp$Z2
    # Decide whether to accept new values
    log.alpha2 <- (
      ell(comparisons.1to3, comparisons.2to3, n1, n2, n3, m.curr, u.curr, Z.curr, Z2.prop)
      - ell(comparisons.1to3, comparisons.2to3, n1, n2, n3, m.curr, u.curr, Z.curr, Z2.curr)
      + calc.log.Z2prior(n1, n2, n3, Z2.prop, Z.curr, aBM, bBM)
      - calc.log.Z2prior(n1, n2, n3, Z2.curr, Z.curr, aBM, bBM)
      + log(tmp$mod) # Modifier for proposal probability from the informed proposal
    )
    if (-rexp(1) < log.alpha2) {
      # accept
      Z2.samples[,i] <- Z2.prop
    } else {
      # reject
      Z2.samples[,i] <- Z2.curr
      accepted[2,i] <- 0
    }
    # Print updates
    if (i  %% max(nIter %/% 20, 1) == 0) {
      cat("Iteration ", i, "\n")
    }
  }
  # 7. Burn the initial samples from both the bipartite and tripartite links, and...
  keptsamples <- setdiff(1:nIter,seq_len(burn))
  # 8. Return new list comprised of only accepted samples and new matchings, Z2
  return(list(Z1=Z.samples[,keptsamples,drop=FALSE],
              Z2=Z2.samples[,keptsamples,drop=FALSE],
              m=m.samples[,keptsamples,drop=FALSE],
              u=u.samples[,keptsamples,drop=FALSE],
              accepted=accepted[,keptsamples,drop=FALSE]))
}

#' @export
tripartiteRL <- function(df1, df2, df3,
                         flds=NULL, flds1=NULL, flds2=NULL, flds3=NULL,
                         types=NULL, breaks=c(0,.25,.5),
                         trace=FALSE,
                         nIter=1000, burn=round(nIter*.1), a=1, b=1, aBM=1, bBM=1, seed=0) {
  # Create comparison data using BRL::compareRecords
  cmpdata.1to2 <- compareRecords(df1, df2, flds, flds1=flds1, flds2=flds2,
                                 types=types, breaks=breaks)
  cmpdata.1to3 <- compareRecords(df1, df3, flds, flds1=flds1, flds2=flds3,
                                 types=types, breaks=breaks)
  cmpdata.2to3 <- compareRecords(df2, df3, flds, flds1=flds2, flds2=flds3,
                                 types=types, breaks=breaks)
  # Call function with precompared files
  tripartiteRL.precmp(cmpdata.1to2, cmpdata.1to3, cmpdata.2to3, trace,
                      nIter, burn, a, b, aBM, bBM, seed)
}


