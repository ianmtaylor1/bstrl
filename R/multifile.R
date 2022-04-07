
#' Perform multifile record linkage via Gibbs sampling "from scratch"
#'
#' @param files A list of files
#' @param flds Vector of names of the fields on which to compare the records in
#'   each file
#' @param types Types of comparisons to use for each field
#' @param breaks Breaks to use for Levenshtein distance on string fields
#' @param nIter,burn MCMC run length parameters. The returned number of samples
#'   is nIter - burn.
#' @param a,b Prior parameters for m and u, respectively.
#' @param aBM,bBM Prior parameters for beta-linkage prior.
#' @param proposals Which kind of full conditional proposals to use for the link
#'   vectors.
#' @param blocksize What blocksize to use for locally balanced proposals. By
#'   default, LB proposals are not blocked
#' @param seed Random seed to set at beginning of MCMC run
#' @param refresh How often to output an update including the iteration number
#'   and percent complete. If refresh >= 1, taken as a number of iterations
#'   between messages (rounded). If 0 < refresh < 1, taken as the proportion of
#'   nIter. If refresh == 0, no messages are displayed.
#' @param maxtime Amount of time, in seconds, after which the sampler will
#'   terminate with however many samples it has produced up to that point. The
#'   sample matrix columns for any unproduced samples will be filled with NAs
#'
#' @return An object of class "bstrlstate"
#'
#' @export
multifileRL <- function(files, flds=NULL, types=NULL, breaks=c(0,.25,.5),
                        nIter=1000, burn=round(nIter*.1), a=1, b=1, aBM=1, bBM=1,
                        proposals=c("component", "LB"), blocksize=NULL, seed=0,
                        refresh=0.1, maxtime=Inf) {
  stopifnot(length(files) >= 2)
  proposals <- match.arg(proposals)
  set.seed(seed)

  if (is.null(refresh)) {
    refresh <- 0
  } else if ((0 < refresh) && (refresh < 1)) {
    refresh <- ceiling(refresh * nIter)
  } else {
    refresh <- round(refresh)
  }

  # Create comparison data from the list of files
  `%do%` <- foreach::`%do%`
  cmpdata <- foreach::foreach(j=seq(2, length(files))) %do% {
    foreach::foreach(i=seq_len(j-1)) %do% {
      compareRecords(files[[i]], files[[j]], flds=flds, types=types, breaks=breaks)
    }
  }

  filesizes <- rep(0, length(files))
  for (f in seq_along(files)) {
    filesizes[f] <- nrow(files[[f]])
  }

  # Metadata for storing in result
  priors <- list(a=a, b=b, aBM=aBM, bBM=bBM)
  cmpdetails <- list(fldsglobal=flds, types=types, breaks=breaks,
                     flds=rep(list(NULL), length(files)))

  # Initial state and save arrays
  slcurr <- streaminglinks(filesizes)

  msave <- usave <- matrix(NA, nrow=sum(cmpdata[[1]][[1]]$nDisagLevs), ncol=nIter)
  Zsave <- matrix(NA, nrow=length(savestate(slcurr)), ncol=nIter)

  samplingstarttime <- burnendtime <- Sys.time() # Start the clock to time sampling

  # Perform gibbs sampling
  completed <- T
  for (iter in seq_len(nIter)) {
    tmp <- r_m_u_fc_smcmc(cmpdata, slcurr, a, b)
    mcurr <- tmp$m
    ucurr <- tmp$u

    for (f in seq(2, length(files))) {
      if (proposals == "component") {
        slcurr <- draw.Z.componentwise(f, cmpdata, slcurr, mcurr, ucurr, aBM, bBM)
      } else if (proposals == "LB") {
        slcurr <- draw.Z.locbal(f, cmpdata, slcurr, mcurr, ucurr, aBM, bBM,
                                blocksize=blocksize)
      }
    }

    msave[,iter] <- mcurr
    usave[,iter] <- ucurr
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

    # Should we break out of sampling for exceeding time?
    if (as.double(Sys.time() - samplingstarttime, units="secs") > maxtime) {
      completed <- F
      break
    }
  }

  samplingendtime <- Sys.time() # Stop the clock
  if (!completed) {
    burnendtime <- samplingendtime
  }

  # Make note of how many iterations we actually reached
  if (!completed) {
    maxiter <- iter
  } else {
    maxiter <- nIter
  }

  # Post-process samples into summary statistics for m and u full conditionals
  m.fc.pars <- u.fc.pars <- matrix(NA, nrow=nrow(msave), ncol=nIter)
  for (s in seq_len(maxiter)) {
    tmp <- disag.counts.allfiles(cmpdata, streaminglinks(filesizes, Zsave[,s]))
    m.fc.pars[,s] <- tmp$match
    u.fc.pars[,s] <- tmp$nonmatch
  }

  # Create diagnostics and filters to burn
  diagnostics <- list(
    completed = completed,
    burntime = as.double(burnendtime - samplingstarttime, units="secs"),
    samplingtime = as.double(samplingendtime - burnendtime, units="secs")
  )
  if (completed) {
    # If we finished, burn as usual and include normal diagnostics
    iterfilter <- setdiff(seq_len(nIter), seq_len(burn))
    continue <- NULL
  } else {
    # If we didn't finish, don't burn yet and include info necessary for
    # resuming later
    iterfilter <- seq_len(nIter)
    continue <- list(
      maxiter = maxiter,
      needtoburn = burn,
      endingseed = .Random.seed,
      proposals = proposals,
      blocksize = blocksize
    )
  }
  # Burn and package results into final structure
  structure(
    list(
      m = msave[,iterfilter, drop=F],
      u = usave[,iterfilter, drop=F],
      Z = Zsave[,iterfilter, drop=F],
      files = files,
      comparisons = cmpdata,
      cmpdetails = cmpdetails,
      priors = priors,
      m.fc.pars = m.fc.pars[,iterfilter, drop=F],
      u.fc.pars = u.fc.pars[,iterfilter, drop=F],
      diagnostics = diagnostics,
      continue = continue
    ),
    class = "bstrlstate"
  )
}
