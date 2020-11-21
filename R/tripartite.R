# This file contains the main functions to do tripartite (or streaming) linkage

# Tripartite RL function built on the BRL package. Arguments the same as
# BRL.gibbs, but with the additional df3 and flds3 arguments.
# DESIGN CHOICES
#  - Discard burn in BRL.gibbs initially, or pass burn=0 and discard burn
#    at end of this function?
#' @export
tripartiteRL.precmp <- function(cmpdata.1to2, cmpdata.1to3, cmpdata.2to3, trace=FALSE,
                                bipartite.samp=NULL,
                                nIter.bi=1200, burn.bi=round(nIter.bi*.1), bipartite.method="BRL",
                                nIter.tri=nIter.bi-burn.bi, burn.tri=round(nIter.tri*.1),
                                pprb.method="ordered", Z2blocksize=NULL,
                                a=1, b=1, aBM=1, bBM=1, seed=0) {
  # Parameter checking
  if (! is.element(pprb.method, c("ordered","permuted","resampled"))) {
    stop("pprb.method must be one of 'ordered', 'permuted', or 'resampled'")
  }
  if (is.element(pprb.method, c("ordered","permuted")) && (nIter.tri != nIter.bi - burn.bi)) {
    warntext <- paste0("If pprb.method is '", pprb.method, "', ",
                       "we require nIter.tri == nIter.bi - burn.bi\n",
                       "Setting parameter values:\n",
                       "nIter.tri <- ", nIter.bi - burn.bi, "\n",
                       "burn.tri  <- ", burn.tri)
    warning(warntext)
    nIter.tri <- nIter.bi - burn.bi
  }

  # 1. Size of files
  n1 <- cmpdata.1to3$n1
  n2 <- cmpdata.2to3$n1
  n3 <- cmpdata.1to3$n2

  # 2. Do bipartite RL between df1 and df2 using bipartiteRL.gibbs()
  #    Do not burn anything now! Keep it all, and burn later
  if (is.null(bipartite.samp)) {
    cat("Beginning bipartite sampling\n")
    bipartite.samp <- bipartiteRL.precmp(cmpdata.1to2, nIter.bi, burn.bi, a, b, aBM, bBM, seed, method=bipartite.method, blocksize=Z2blocksize)
  } else {
    cat("Bipartite samples passed in, using those. Other bipartite sampling",
        "options will be ignored.\n")
    # Reset nIter.bi and burn.bi based on the passed in samples
    burn.bi <- 0
    nIter.bi <- ncol(bipartite.samp$Z)
    # Would normally depend on bipartite sampling function to set seed. Do that
    # here instead
    set.seed(seed)
  }


  # 3. Create precomputed data structures
  comparisons.1to3 <- preproc.cmpdata(cmpdata.1to3)
  comparisons.2to3 <- preproc.cmpdata(cmpdata.2to3)
  cmpdata.list <- list(comparisons.1to3, comparisons.2to3)
  nDisagLevs <- comparisons.1to3$nDisagLevs

  # 4. Set up empty arrays of the appropriate size that will eventually be returned
  m.samples  <- matrix(0, nrow=nrow(bipartite.samp$m), ncol=nIter.tri)
  u.samples  <- matrix(0, nrow=nrow(bipartite.samp$u), ncol=nIter.tri)
  Z.samples  <- matrix(0, nrow=nrow(bipartite.samp$Z), ncol=nIter.tri)
  Z2.samples <- matrix(0, nrow=n3,                     ncol=nIter.tri)

  # 5. If PPRB is permuted, produce a list of the proposal indexes by permuting
  if (pprb.method == "permuted") {
    pprb.index <- sample(nIter.tri) # Note: nIter.tri == nIter.bi - burn.bi necessarily
  }

  # 6. Initialize the first sample (m,u,Z: first sample from bipartite. Z2: unlinked)
  if (pprb.method == "resampled") {
    startidx <- sample(nIter.bi - burn.bi, size=1)
  } else if (pprb.method == "permuted") {
    startidx <- pprb.index[1]
  } else { # Ordered
    startidx <- 1
  }
  m.samples[,1] <- bipartite.samp$m[,startidx]
  u.samples[,1] <- bipartite.samp$u[,startidx]
  Z.samples[,1] <- bipartite.samp$Z[,startidx]
  Z.curr.idx <- startidx # Keep track of index of current Z in previous samples
  Z2.samples[,1] <- n1 + n2 + seq_len(n3)

  # Create array to hold record of acceptances: first row = m,u,Z. second row = Z2
  # Contains all 1's by default, zeros written on rejection
  accepted <- matrix(1, nrow=2, ncol=nIter.tri)
  rownames(accepted) <- c("PPRB Z", "Locally Balanced Z2")
  # Creates array to track impossible proposals: proposals of m,u,Z from pprb
  # which have probability zero
  pprb.impossible <- rep(0, nIter.tri)
  # Track how different Z1 proposals are from the current value
  Zprop.diffs <- rep(0, nIter.tri)
  # Track the PPRB acceptance ratio
  pprb.log.ratio <- rep(0, nIter.tri)

  # 7. Perform M-H for each iteration
  cat("Beginning tripartite sampling\n")
  for (i in seq(2, nIter.tri)) {
    # What are the current values of the parameters?
    m.curr <- m.samples[,i-1]
    u.curr <- u.samples[,i-1]
    Z.curr <- Z.samples[,i-1]
    Z2.curr <- Z2.samples[,i-1]

    # 7.1: Propose new Z from PPRB
    if (pprb.method == "resampled") {
      # Go through all bipartite samples until you find one that is a valid
      # link state, then propose that one. There is guaranteed to always be at
      # least one because the current state is valid. All invalid link states
      # have probability zero in the target so we can avoid proposing them.
      Z.prop.idx <- Z.curr.idx # Default: propose current state
      for (j in sample(ncol(bipartite.samp$Z))) {
        if (j == Z.curr.idx) { # Skip current state
          next
        }
        if (valid.link.state(n1, bipartite.samp$Z[,j], Z2.curr)) {
          Z.prop.idx <- j
          break
        }
      }
      # If there were literally no possible states, we are now proposing the
      # current state
    } else {
      # What value should be proposed?
      if (pprb.method == "permuted") {
        Z.prop.idx <- pprb.index[i]
      } else {
        # Ordered
        Z.prop.idx <- i
      }
    }
    Z.prop <- bipartite.samp$Z[,Z.prop.idx]
    # How different is z.prop from z.curr?
    Zprop.diffs[i] <- sum(Z.prop != Z.curr)
    # What is the log of the MH acceptance ratio?
    log.alpha1 <- (
      calc.log.lkl(cmpdata.list, m.curr, u.curr, Z.prop, Z2.curr, do.trace=trace)
      - calc.log.lkl(cmpdata.list, m.curr, u.curr, Z.curr, Z2.curr, do.trace=trace)
      + calc.log.Z2prior(n1, Z2.curr, Z.prop, aBM, bBM)
      - calc.log.Z2prior(n1, Z2.curr, Z.curr, aBM, bBM)
      + ddirichlet.multi(m.curr, bipartite.samp$m.fc.pars[,Z.prop.idx] + a, nDisagLevs, log=TRUE)
      + ddirichlet.multi(u.curr, bipartite.samp$u.fc.pars[,Z.prop.idx] + b, nDisagLevs, log=TRUE)
      - ddirichlet.multi(m.curr, bipartite.samp$m.fc.pars[,Z.curr.idx] + a, nDisagLevs, log=TRUE)
      - ddirichlet.multi(u.curr, bipartite.samp$u.fc.pars[,Z.curr.idx] + b, nDisagLevs, log=TRUE)
    )
    # Store the ratio
    pprb.log.ratio[i] <- log.alpha1
    # Decide whether to accept new values
    if (-rexp(1) < log.alpha1) {
      # accept
      Z.samples[,i] <- Z.prop
      Z.curr.idx <- Z.prop.idx
    } else {
      # reject
      Z.samples[,i] <- Z.curr
      accepted[1,i] <- 0
      if (log.alpha1 == -Inf) {
        pprb.impossible[i] <- 1
      }
    }
    # Reset "current" values of parameters
    Z.curr <- Z.samples[,i]

    # 7.2: Draw new m, u from full conditional distributions
    tmp <- r_m_u_fc(cmpdata.list, Z.curr, Z2.curr, a, b,
                    bipartite.samp$m.fc.pars[,Z.curr.idx],
                    bipartite.samp$u.fc.pars[,Z.curr.idx])
    m.samples[,i] <- m.curr <- tmp$m
    u.samples[,i] <- u.curr <- tmp$u

    # 7.3: Propose new Z2
    tmp <- draw.Z2.informed(cmpdata.list, Z.curr, Z2.curr, m.curr, u.curr,
                            aBM, bBM, trace=trace, blocksize=Z2blocksize)
    Z2.prop <- tmp$Z2
    # Decide whether to accept new values
    log.alpha2 <- (
      calc.log.lkl(cmpdata.list, m.curr, u.curr, Z.curr, Z2.prop, do.trace=trace)
      - calc.log.lkl(cmpdata.list, m.curr, u.curr, Z.curr, Z2.curr, do.trace=trace)
      + calc.log.Z2prior(n1, Z2.prop, Z.curr, aBM, bBM)
      - calc.log.Z2prior(n1, Z2.curr, Z.curr, aBM, bBM)
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
    if (i  %% max(nIter.tri %/% 20, 1) == 0) {
      cat("Iteration ", i, "\n")
    }
  }
  # 8. Burn the initial samples from both the bipartite and tripartite links, and...
  keptsamples <- setdiff(1:nIter.tri,seq_len(burn.tri))
  # 9. Return new list comprised of only accepted samples and new matchings, Z2
  return(list(Z1=Z.samples[,keptsamples,drop=FALSE],
              Z2=Z2.samples[,keptsamples,drop=FALSE],
              m=m.samples[,keptsamples,drop=FALSE],
              u=u.samples[,keptsamples,drop=FALSE],
              accepted=accepted[,keptsamples,drop=FALSE],
              pprb.impossible=pprb.impossible[keptsamples],
              Zprop.diffs=Zprop.diffs[keptsamples],
              pprb.log.ratio=pprb.log.ratio[keptsamples]))
}

#' @export
tripartiteRL <- function(df1, df2, df3,
                         flds=NULL, flds1=NULL, flds2=NULL, flds3=NULL,
                         types=NULL, breaks=c(0,.25,.5),
                         ...) {
  # Create comparison data using BRL::compareRecords
  cmpdata.1to2 <- BRL::compareRecords(df1, df2, flds, flds1=flds1, flds2=flds2,
                                 types=types, breaks=breaks)
  cmpdata.1to3 <- BRL::compareRecords(df1, df3, flds, flds1=flds1, flds2=flds3,
                                 types=types, breaks=breaks)
  cmpdata.2to3 <- BRL::compareRecords(df2, df3, flds, flds1=flds2, flds2=flds3,
                                 types=types, breaks=breaks)
  # Call function with precompared files
  tripartiteRL.precmp(cmpdata.1to2, cmpdata.1to3, cmpdata.2to3, ...)
}


################################################################################
# Helper Functions for Tripartite Linkage Functions ############################
################################################################################

