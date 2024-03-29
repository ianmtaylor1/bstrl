
#' Perform a PPRB update of record linkage with a new file
#'
#' @param state Existing record linkage state. Returned by either bipartiteRL,
#'   PPRBupdate, SMCMCupdate, or multifileRL.
#' @param newfile A data.frame representing the new file: one row per record
#' @param flds Names of fields in the new file to use for comparison. Only used
#'   if no global field names were specified in bipartiteRL initially.
#' @param nIter Number of iterations for which to run the PPRB sampler. By
#'   default, this is the same as the number of samples present in 'state'.
#' @param burn Number of initial iterations to discard. The total number of
#'   samples returned is nIter - burn.
#' @param blocksize Size of blocks to use for locally balanced proposals.
#'   Default performs unblocked locally balanced proposals.
#' @param threestep Whether to perform three Gibbs sampling steps per iteration,
#'   with past Z's updated with PPRB, m and u updated with full conditionals,
#'   and the current Z updated with locally balanced proposals. If false, a two
#'   step Gibbs sampler is used where past Z's, m and u are updated together
#'   using PPRB and the current Z is updated with locally balanced proposals
#' @param seed Random seed to set at the beginning of the MCMC run
#' @param refresh How often to output an update including the iteration number
#'   and percent complete. If refresh >= 1, taken as a number of iterations
#'   between messages (rounded). If 0 < refresh < 1, taken as the proportion of
#'   nIter. If refresh == 0, no messages are displayed.
#'
#' @return An object of class 'bstrlstate' containing posterior samples and
#'   necessary metadata for passing to future streaming updates.
#'
#' @examples
#' data(geco_small)
#' data(geco_small_result)
#'
#' # Add fifth file to previous four-file link result
#' file5.result <- PPRBupdate(geco_small_result, geco_small[[5]],
#'                            nIter=2, burn=1) # Very small run for example
#'
#' @export
PPRBupdate <- function(state, newfile, flds = NULL, nIter = NULL, burn = 0, blocksize = NULL,
                       threestep = TRUE, seed=0, refresh=0.1) {
  set.seed(seed)

  if (is.null(refresh)) {
    refresh <- 0
  } else if ((0 < refresh) && (refresh < 1)) {
    refresh <- ceiling(refresh * nIter)
  } else {
    refresh <- round(refresh)
  }

  # Create and append new comparison data, building up the triangular
  # list-of-lists format
  newcmps <- list()
  for (i in seq_along(state$files)) {
    f <- state$files[[i]]
    flds1 <- state$cmpdetails$flds[[i]]
    thiscmp <- compareRecords(
      df1 = f, df2 = newfile, flds = state$cmpdetails$fldsglobal,
      flds1 = flds1, flds2 = flds, types = state$cmpdetails$types,
      breaks = state$cmpdetails$breaks
    )
    newcmps <- c(newcmps, list(thiscmp))
  }
  cmpdata <- c(state$comparisons, list(newcmps))
  # Append the file to the list of files
  files <- c(state$files, list(newfile))
  # Append the new comparison field names to the comparison details
  cmpdetails <- state$cmpdetails
  cmpdetails$flds <- c(cmpdetails$flds, list(flds))

  # Default number of iterations is same as object passed in.
  if (is.null(nIter)) {
    nIter <- ncol(state$Z)
  }

  # Groupings for m and u
  nDisagLevs <- cmpdata[[1]][[1]]$nDisagLevs

  # Initialize streaming link object. Unlinked to new file, randomly chosen
  # value from existing state
  filesizes <- rep(0, length(files))
  for (i in seq_along(files)) {
    filesizes[i] <- nrow(files[[i]])
  }
  pprb.index.curr <- sample(ncol(state$Z), 1)
  slcurr <- streaminglinks(filesizes)
  Zpre <- state$Z[,pprb.index.curr]
  slcurr <- swapprefix(slcurr, Zpre, conflict = "error")
  # Initialize m and u
  mcurr <- state$m[,pprb.index.curr]
  ucurr <- state$u[,pprb.index.curr]

  # Initialize matrices to save new posterior samples and other quantities
  msave <- matrix(NA, nrow=nrow(state$m), ncol=nIter)
  usave <- matrix(NA, nrow=nrow(state$u), ncol=nIter)
  Zsave <- matrix(NA, nrow=length(savestate(slcurr)), ncol=nIter)
  pprb.index.save <- rep(NA, nIter)
  pprb.accepted <- rep(FALSE, nIter)

  samplingstarttime <- burnendtime <- Sys.time() # Start the clock to time sampling

  # Main PPRB process
  for (iter in seq_len(nIter)) {
    # New pprb proposal index
    pprb.index.prop <- sample(ncol(state$Z), 1)
    if (threestep) {
      # Sample m and u from full conditional
      tmp <- r_m_u_fc_pprb(newcmps, slcurr, state$priors$a, state$priors$b,
                           state$m.fc.pars[,pprb.index.curr], state$u.fc.pars[,pprb.index.curr])
      mcurr <- tmp$m
      ucurr <- tmp$u
      # Sample previous Z's with PPRB
      slprop <- swapprefix(slcurr, state$Z[,pprb.index.prop], conflict = "null")
      if (!is.null(slprop)) { # Check for impossible proposals
        log.alpha <- (
          calc.log.lkl.lastfile(newcmps, mcurr, ucurr, slprop) -
          calc.log.lkl.lastfile(newcmps, mcurr, ucurr, slcurr) +
          log.Zprior(slprop, state$priors$aBM, state$priors$bBM, vec="last") -
          log.Zprior(slcurr, state$priors$aBM, state$priors$bBM, vec="last") +
          ddirichlet.multi(mcurr, state$m.fc.pars[,pprb.index.prop] + state$priors$a, nDisagLevs, log=T) +
          ddirichlet.multi(ucurr, state$u.fc.pars[,pprb.index.prop] + state$priors$b, nDisagLevs, log=T) -
          ddirichlet.multi(mcurr, state$m.fc.pars[,pprb.index.curr] + state$priors$a, nDisagLevs, log=T) -
          ddirichlet.multi(ucurr, state$u.fc.pars[,pprb.index.curr] + state$priors$b, nDisagLevs, log=T)
        )

        if (log(stats::runif(1)) < log.alpha) {
          slcurr <- slprop
          pprb.index.curr <- pprb.index.prop
          pprb.accepted[iter] <- TRUE
        }
      }
    } else { # two step
      # Sample m, u, and previous Z's together with PPRB
      slprop <- swapprefix(slcurr, state$Z[,pprb.index.prop], conflict = "null")
      mprop <- state$m[,pprb.index.prop]
      uprop <- state$u[,pprb.index.prop]
      if (!is.null(slprop)) { # Check for impossible proposals
        log.alpha <- (
          calc.log.lkl.lastfile(newcmps, mprop, uprop, slprop) -
          calc.log.lkl.lastfile(newcmps, mcurr, ucurr, slcurr) +
          log.Zprior(slprop, state$priors$aBM, state$priors$bBM, vec="last") -
          log.Zprior(slcurr, state$priors$aBM, state$priors$bBM, vec="last")
        )

        if (log(stats::runif(1)) < log.alpha) {
          mcurr <- mprop
          ucurr <- uprop
          slcurr <- slprop
          pprb.index.curr <- pprb.index.prop
          pprb.accepted[iter] <- T
        }
      }
    }
    msave[,iter] <- mcurr
    usave[,iter] <- ucurr
    pprb.index.save[iter] <- pprb.index.curr

    # Sample Zk, the latest Z vector using LB proposals
    slcurr <- draw.Z.locbal.lastfile(newcmps, slcurr, mcurr, ucurr,
                                     state$priors$aBM, state$priors$bBM, blocksize=blocksize)
    Zsave[,iter] <- savestate(slcurr)

    # Update messages, if desired
    if ((refresh > 0) && (iter %% refresh == 0)) {
      message(
        iter, "/", nIter,
        " [", round(100*iter/nIter), "%]",
        if (iter <= burn) " (burn)" else ""
      )
    }

    # Log the end of the burn-in phase
    if (iter == burn) {
      burnendtime <- Sys.time()
    }
  }

  samplingendtime <- Sys.time() # Stop the clock

  # Store updated summary statistics of comparisons for use later
  m.fc.pars <- matrix(0, nrow=nrow(msave), ncol=nIter)
  u.fc.pars <- matrix(0, nrow=nrow(msave), ncol=nIter)
  for (i in 1:nIter) {
    tmp <- disag.counts.lastfile(newcmps, streaminglinks(filesizes, Zsave[,i]))
    m.fc.pars[,i] <- tmp$match + state$m.fc.pars[,pprb.index.save[i]]
    u.fc.pars[,i] <- tmp$nonmatch + state$u.fc.pars[,pprb.index.save[i]]
  }

  # Burn and Construct and return the new posterior state
  iterfilter <- setdiff(seq_len(nIter), seq_len(burn))
  structure(
    list(
      Z = Zsave[,iterfilter,drop=F],
      m = msave[,iterfilter,drop=F],
      u = usave[,iterfilter,drop=F],
      files = files,
      comparisons = cmpdata,
      priors = state$priors,
      cmpdetails = cmpdetails,
      m.fc.pars = m.fc.pars[,iterfilter,drop=F],
      u.fc.pars = u.fc.pars[,iterfilter,drop=F],
      diagnostics = list(
        burntime = as.double(burnendtime - samplingstarttime, units="secs"),
        samplingtime = as.double(samplingendtime - burnendtime, units="secs"),
        pprb.accepted = pprb.accepted[iterfilter]
      )
    ),
    class = "bstrlstate"
  )
}
