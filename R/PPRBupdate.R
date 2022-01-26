
#' @export
PPRBupdate <- function(state, newfile, flds = NULL, nIter = NULL, burn = 0, blocksize = NULL,
                       threestep = T, seed=0) {
  set.seed(seed)

  # Create and append new comparison data, building up the triangular
  # list-of-lists format
  newcmps <- list()
  for (i in seq_along(state$files)) {
    f <- state$files[[i]]
    flds1 <- state$cmpdetails$flds[[i]]
    thiscmp <- BRL::compareRecords(
      df1 = f, df2 = newfile, flds = state$cmpdetails$fldsglobal,
      flds1 = flds1, flds2 = flds2, types = state$cmpdetails$types,
      breaks = state$cmpdetails$breaks
    )
    newcmps <- c(newcmp, list(thiscmp))
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

  # Initialize streaming link object. Unlinked to new file, randomly chosen
  # value from existing state
  filesizes <- c()
  for (i in seq_along(files)) {
    filesizes[i] <- nrow(files[i])
  }
  slcurr <- streaminglinks(filesizes)
  Zpre <- state$Z[,sample(ncol(state$Z), 1)]
  slcurr <- swapprefix(sl, Zpre, conflict = "error")
  # Initialize m and u
  mcurr <- rep(NA, nrow(state$m))
  ucurr <- rep(NA, nrow(state$u))

  # Initialize matrices to save new posterior samples
  msave <- matrix(NA, nrow=nrow(state$m), ncol=nIter)
  usave <- matrix(NA, nrow=nrow(state$u), ncol=nIter)
  Zsave <- matrix(NA, nrow=length(savestate(sl)), ncol=nIter)

  # Main PPRB process
  for (iter in seq_len(nIter)) {
    if (threestep) {
      # Sample m and u from full conditional
      tmp <- list(m=mcurr, u=ucurr) # TODO
      mcurr <- tmp$m
      ucurr <- tmp$u
      # Sample previous Z's with PPRB
      pprbidx <- sample(ncol(state$Z), 1)
      slprop <- swapprefix(sl, state$Z[,pprbidx], conflict = "null")
      if (!is.null(slprop)) { # Check for impossible proposals
        log.alpha <- 0 # TODO

        if (log(runif(1)) < log.alpha) {
          slcurr <- slprop
        }
      }
    } else {
      # Sample m, u, and previous Z's together with PPRB
      pprbidx <- sample(ncol(state$Z), 1)
      slprop <- swapprefix(sl, state$Z[,pprbidx], conflict = "null")
      mprop <- state$m[,pprbidx]
      uprop <- state$u[,pprbidx]
      if (!is.null(slprop)) { # Check for impossible proposals
        log.alpha <- 0 # TODO

        if (log(runif(1)) < log.alpha) {
          mcurr <- mprop
          ucurr <- uprop
          slcurr <- slprop
        }
      }
    }
    msave[,iter] <- mcurr
    usave[,iter] <- ucurr

    # Sample Zk, the latest Z vector using LB proposals
    # TODO
    Zsave[,iter] <- savestate(slcurr)
  }

  # Construct and return the new posterior state
  structure(
    list(
      Z = Zsave,
      m = msave,
      u = uszve,
      files = files,
      comparisons = cmpdata,
      priors = state$priors,
      cmpdetails = cmpdetails,
    ),
    class = "bstrlstate"
  )
}
