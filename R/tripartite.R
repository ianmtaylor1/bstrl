# This file contains the main functions to do tripartite (or streaming) linkage

# Tripartite RL function built on the BRL package. Arguments the same as
# BRL.gibbs, but with the additional df3 and flds3 arguments.
#' @export
tripartiteRL.precmp <- function(cmpdata.1to2, cmpdata.1to3, cmpdata.2to3, trace=FALSE,
                                nIter.bi=1200, burn.bi=round(nIter.bi*.1),
                                nIter.tri=nIter.bi-burn.bi, burn.tri=round(nIter.tri*.1),
                                pprb.method="ordered", pprb.adapt=FALSE, pprb.adapt.goal=0.23,
                                pprb.adapt.iter=burn.tri, Z2blocksize=NULL,
                                a=1, b=1, aBM=1, bBM=1, seed=0,
                                print=TRUE) {
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
  if ((pprb.method != "resampled") && pprb.adapt) {
    warntext <- paste0("If pprb.method is '", pprb.method, "', ",
                       "we require pprb.adapt == FALSE\n",
                       "Setting parameter values:\n",
                       "pprb.adapt <- FALSE")
    warning(warntext)
    pprb.adapt <- FALSE
  }
  if (pprb.adapt && (pprb.adapt.iter > burn.tri)) {
    warntext <- paste0("If pprb.adapt is true but burn.tri < pprb.adapt.iter\n",
                       "Setting parameter values:\n",
                       "pprb.adapt.iter <- ", burn.tri)
    warning(warntext)
    pprb.adapt.iter <- burn.tri
  }

  # 1. Size of files
  n1 <- cmpdata.1to3$n1
  n2 <- cmpdata.2to3$n1
  n3 <- cmpdata.1to3$n2

  # 2. Do bipartite RL between df1 and df2 using BRL.gibbs()
  if (print) {
    cat("Beginning bipartite linkage...\n")
  }
  bipartite.samp <- BRL.gibbs.precmp(cmpdata.1to2, nIter.bi, burn.bi, a, b, aBM, bBM, seed)

  # 3. Create precomputed data structures
  comparisons.1to3 <- preproc.cmpdata(cmpdata.1to3)
  comparisons.2to3 <- preproc.cmpdata(cmpdata.2to3)

  # 4. Set up empty arrays of the appropriate size that will eventually be returned
  m.samples  <- matrix(0, nrow=nrow(bipartite.samp$m), ncol=nIter.tri)
  u.samples  <- matrix(0, nrow=nrow(bipartite.samp$u), ncol=nIter.tri)
  Z.samples  <- matrix(0, nrow=nrow(bipartite.samp$Z), ncol=nIter.tri)
  Z2.samples <- matrix(0, nrow=n3,                     ncol=nIter.tri)

  # 5. Specific variables for different pprb.methods
  # If PPRB is permuted, produce a list of the proposal indexes by permuting
  if (pprb.method == "permuted") {
    pprb.index <- sample(nIter.tri) # Note: nIter.tri == nIter.bi - burn.bi necessarily
  }
  # How many times should the PPRB step be repeated? Will be potentially
  # increased if pprb.method is 'resampled' and pprb.adapt is true
  pprb.repeat <- 1

  # 6. Initialize the first sample (m,u,Z: first sample from bipartite. Z2: unlinked)
  if (pprb.method == "resampled") {
    startidx <- sample(nIter.bi - burn.bi, size=1)
  } else if (pprb.method == "permuted") {
    startidx <- pprb.index[1]
  } else { # Ordered
    startidx <- 1
  }
  m.samples[,1] <- m.curr <- bipartite.samp$m[,startidx]
  u.samples[,1] <- u.curr <- bipartite.samp$u[,startidx]
  Z.samples[,1] <- Z.curr <- bipartite.samp$Z[,startidx]
  Z2.samples[,1] <- Z2.curr <- n1 + n2 + seq_len(n3)

  # Create array to hold record of acceptances: first row = m,u,Z. second row = Z2
  # Contains all 0's by default, oness written on accept
  accepted <- matrix(0, nrow=2, ncol=nIter.tri)
  accepted[,1] <- 1
  # Creates array to track impossible proposals: proposals of m,u,Z from pprb
  # which have probability zero
  pprb.impossible <- rep(0, nIter.tri)

  # Choose the likelihood function based on tracing flag
  if (trace) {
    ell <- calc.log.lkl.tracing
  } else {
    ell <- calc.log.lkl
  }


  # 7. Perform M-H for each iteration
  if (print) {
    cat("Beginning tripartite linkage...\n")
    cat(paste0("Options:\nPPRB method = ",pprb.method,"\n"))
    if(pprb.adapt) {
      cat("Adaptation: on, ")
      cat(paste0("Goal = ", round(100*pprb.adapt.goal, 0), "% after ",
                 pprb.adapt.iter, " iterations\n"))
    } else {
      cat("Adaptation: off\n")
    }
  }
  for (i in seq(2, nIter.tri)) {

    # 7.1: Propose new m,u,Z collectively from PPRB
    for (pprb.step in 1:pprb.repeat) {
      # Which values should be proposed?
      if (pprb.method == "resampled") {
        # Go through all bipartite samples until you find one that is a valid
        # link state, then propose that one. There is guaranteed to always be at
        # least one because the current state is valid. All invalid link states
        # have probability zero in the target so we can avoid proposing them.
        for (j in sample(ncol(bipartite.samp$Z))) {
          if (valid.link.state(n1, bipartite.samp$Z[,j], Z2.curr)) {
            break # Leaves j at its current value outside the loop
          }
        }
        m.prop <- bipartite.samp$m[,j]
        u.prop <- bipartite.samp$u[,j]
        Z.prop <- bipartite.samp$Z[,j]
      } else {
        # What value should be proposed?
        if (pprb.method == "permuted") {
          m.prop <- bipartite.samp$m[,pprb.index[i]]
          u.prop <- bipartite.samp$u[,pprb.index[i]]
          Z.prop <- bipartite.samp$Z[,pprb.index[i]]
        } else {
          # Ordered
          m.prop <- bipartite.samp$m[,i]
          u.prop <- bipartite.samp$u[,i]
          Z.prop <- bipartite.samp$Z[,i]
        }
      }
      # What is the log of the MH acceptance ratio?
      log.alpha1 <- (
        ell(comparisons.1to3, comparisons.2to3, n1, n2, n3, m.prop, u.prop, Z.prop, Z2.curr)
        - ell(comparisons.1to3, comparisons.2to3, n1, n2, n3, m.curr, u.curr, Z.curr, Z2.curr)
        + calc.log.Z2prior(n1, n2, n3, Z2.curr, Z.prop, aBM, bBM)
        - calc.log.Z2prior(n1, n2, n3, Z2.curr, Z.curr, aBM, bBM)
      )
      # Decide whether to accept new values
      if (-rexp(1) < log.alpha1) {
        # accept
        m.samples[,i] <- m.prop
        u.samples[,i] <- u.prop
        Z.samples[,i] <- Z.prop
        accepted[1,i] <- 1
      } else {
        # reject
        m.samples[,i] <- m.curr
        u.samples[,i] <- u.curr
        Z.samples[,i] <- Z.curr
        if (log.alpha1 == -Inf) {
          pprb.impossible[i] <- 1
        }
      }
      # Reset "current" values of parameters
      m.curr <- m.samples[,i]
      u.curr <- u.samples[,i]
      Z.curr <- Z.samples[,i]
    }

    # 7.2: Propose new Z2
    tmp <- draw.Z2.informed(n1, n2, n3, Z.curr, Z2.curr, m.curr, u.curr,
                            comparisons.1to3, comparisons.2to3,
                            aBM, bBM, trace=trace, blocksize=Z2blocksize)
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
      accepted[2,i] <- 1
    } else {
      # reject
      Z2.samples[,i] <- Z2.curr
    }

    # Reset "current" value of Z2 parameter
    Z2.curr <- Z2.samples[,i]

    # 7.3 Periodic check-ins and updates

    # Print 5% periodic updates
    if (print && (i  %% max(nIter.tri %/% 20, 1) == 0)) {
      cat("Iteration ", i, "\n")
    }

    # Handle pprb adaptation period
    if (pprb.adapt && (i == pprb.adapt.iter)) {
      est.accept.rate <- (sum(accepted[1,1:pprb.adapt.iter]) + 1) / (pprb.adapt.iter + 2)
      pprb.repeat <- ceiling(log(1 - pprb.adapt.goal) / log(1 - est.accept.rate))
      if (print) {
        msgtext <- paste0("Adapting PPRB: Estimated acceptance rate = ",
                          round(100*est.accept.rate, 1), "%\n",
                          "Setting repeat = ", pprb.repeat, "\n")
        cat(msgtext)
      }
    }

    # Output informational message about burn
    if (print && (i == burn.tri)) {
      cat("Ending burn-in period.\n")
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
              pprb.impossible=pprb.impossible[keptsamples]))
}

#' @export
tripartiteRL <- function(df1, df2, df3,
                         flds=NULL, flds1=NULL, flds2=NULL, flds3=NULL,
                         types=NULL, breaks=c(0,.25,.5),
                         trace=FALSE,
                         nIter.bi=1200, burn.bi=round(nIter.bi*.1),
                         nIter.tri=nIter.bi-burn.bi, burn.tri=round(nIter.tri*.1),
                         pprb.method="ordered", Z2blocksize=NULL,
                         a=1, b=1, aBM=1, bBM=1, seed=0) {
  # Create comparison data using BRL::compareRecords
  cmpdata.1to2 <- compareRecords(df1, df2, flds, flds1=flds1, flds2=flds2,
                                 types=types, breaks=breaks)
  cmpdata.1to3 <- compareRecords(df1, df3, flds, flds1=flds1, flds2=flds3,
                                 types=types, breaks=breaks)
  cmpdata.2to3 <- compareRecords(df2, df3, flds, flds1=flds2, flds2=flds3,
                                 types=types, breaks=breaks)
  # Call function with precompared files
  tripartiteRL.precmp(cmpdata.1to2, cmpdata.1to3, cmpdata.2to3, trace,
                      nIter.bi, burn.bi, nIter.tri, burn.tri,
                      pprb.method, Z2blocksize,
                      a, b, aBM, bBM, seed)
}


################################################################################
# Helper Functions for Tripartite Linkage Functions ############################
################################################################################

