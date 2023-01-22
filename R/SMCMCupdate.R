
#' Perform an SMCMC update of record linkage with a new file
#'
#' @param state Existing record linkage state. Returned by either bipartiteRL,
#'   PPRBupdate, SMCMCupdate, PRPBWTupdate, or multifileRL. This is used as the
#'   ensemble of samples in SMCMC update
#' @param newfile A data.frame representing the new file: one row per record
#' @param flds Names of fields in the new file to use for comparison. Only used
#'   if no global field names were specified in bipartiteRL initially.
#' @param nIter.jumping,nIter.transition Number of iterations to use in the
#'   jumping kernel and transition kernel, respectively, for each ensemble
#'   sample.
#' @param cores Number of cores to use for parallel execution. If cores == 1,
#'   update is run sequentially. A cluster is created using
#'   parallel::makeCluster().
#' @param proposals.jumping,proposals.transition Which kernel to use for Z
#'   updates in the jumping and transition kernels, respectively.
#' @param blocksize Size of blocks to use for locally balanced proposals.
#'   Default performs unblocked locally balanced proposals.
#' @param seed Random seed to set at the beginning of the MCMC run. This is
#'   ignored if cores > 1.
#'
#' @return An object of class 'bstrlstate' containing posterior samples and
#'   necessary metadata for passing to future streaming updates.
#'
#' @examples
#' data(geco_small)
#' data(geco_small_result)
#'
#' # Add fifth file to previous four-file link result
#' filtered <- thinsamples(geco_small_result, 2) # Filter ensemble to 2 - very small for example
#' file5.result <- SMCMCupdate(filtered, geco_small[[5]],
#'                             nIter.jumping=1, nIter.transition=1, # Very small run for example
#'                             proposals.jumping="LB", proposals.transition="LB",
#'                             blocksize=5)
#'
#' @export
SMCMCupdate <- function(state, newfile, flds=NULL, nIter.jumping=5, nIter.transition=10,
                        cores=1, proposals.jumping=c("component", "LB"),
                        proposals.transition=c("LB", "component"), blocksize=NULL,
                        seed=0) { # Future additions: directratio and fastmu?
  proposals.jumping <- match.arg(proposals.jumping)
  proposals.transition <- match.arg(proposals.transition)

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

  # Initialize streaming link object. Unlinked to new file, randomly chosen
  # value from existing state
  filesizes <- rep(0, length(files))
  for (i in seq_along(files)) {
    filesizes[i] <- nrow(files[[i]])
  }

  # How big is the ensemble we are working with?
  ensemblesize <- ncol(state$Z)

  # Pre-append unlinked z vector for new file
  newZ <- matrix(NA, nrow=nrow(state$Z)+filesizes[length(files)], ncol=ensemblesize)
  for (s in seq_len(ensemblesize)) {
    slcurr <- streaminglinks(filesizes)
    slcurr <- swapprefix(slcurr, state$Z[,s], conflict="error")
    newZ[,s] <- savestate(slcurr)
  }

  ensemble <- list(m=state$m, u=state$u, Z=newZ)

  updated <- coreSMCMCupdate(ensemble, state$priors, files, cmpdata,
                             nIter.jumping, nIter.transition, cores,
                             proposals.jumping, proposals.transition,
                             blocksize, seed)

  # Construct and return the new link state
  updated$files <- files
  updated$comparisons <- cmpdata
  updated$priors <- state$priors
  updated$cmpdetails <- cmpdetails

  return(updated)
}

# Function to perform the core SMCMC update, after file comparisons, and returning
# partial info
# The goal of this function is to be able to be used in SMCMCupdate as well as
# PPRBWTupdate. Therefore it doesn't do any file comparisons, but just performs
# the Jumping and Transition kernels.
#
# Most parameters are passed directly from the SMCMC function, where value
# validation is done. Note these attributes for specific parameters:
#  ensemble - a list. It only needs to contain the starting ensemble for m, u,
#    and Z, but can just be the state passed to SMCMCupdate
#
coreSMCMCupdate <- function(ensemble, priors, files, cmpdata,
                            nIter.jumping, nIter.transition,
                            cores, proposals.jumping, proposals.transition,
                            blocksize,
                            seed) {

  # Pull out comparisons with most recent file
  newcmps <- cmpdata[[length(cmpdata)]]

  # Groupings for m and u
  nDisagLevs <- cmpdata[[1]][[1]]$nDisagLevs

  # Initialize streaming link object. Unlinked to new file, randomly chosen
  # value from existing state
  filesizes <- rep(0, length(files))
  for (i in seq_along(files)) {
    filesizes[i] <- nrow(files[[i]])
  }

  # How big is the ensemble we are working with?
  ensemblesize <- ncol(ensemble$Z)

  # Unpack starting ensemble into list
  samplist <- list()
  for (i in seq_len(ensemblesize)) {
    samplist[[i]] <- list(
      m=ensemble$m[,i],
      u=ensemble$u[,i],
      sl=streaminglinks(filesizes, ensemble$Z[,i])
    )
  }

  # Create cluster if necessary and register for parallel execution
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    if (!is.null(seed)) {
      message("Random seed not set in parallel execution.")
    }
  } else {
    foreach::registerDoSEQ()
    if (!is.null(seed)) {
      set.seed(seed)
    }
  }
  `%dopar%` <- foreach::`%dopar%`

  samplingstart <- Sys.time() # Start the clock to time sampling

  # Now do the SMCMC update on each member of the ensemble
  if (nIter.jumping > 0) {
    samplist <- foreach::foreach(s=samplist, .inorder=TRUE) %dopar% {
      # Initial values
      mcurr <- s$m
      ucurr <- s$u
      slcurr <- s$sl

      # Jumping kernel for links to new file
      for (i in seq_len(nIter.jumping)) {
        if (proposals.jumping == "LB") {
          slcurr <- draw.Z.locbal.lastfile(newcmps, slcurr, mcurr, ucurr,
                                           priors$aBM, priors$bBM,
                                           blocksize=blocksize)
        } else if (proposals.jumping == "component") {
          slcurr <- draw.Z.componentwise(length(filesizes), cmpdata, slcurr, mcurr,
                                         ucurr, priors$aBM, priors$bBM)
        }
      }
    }
  }

  jumpingend <- Sys.time() # Record the end of jumping kernel

  # Transition kernel for all parameters. Outer loop is number of iterations,
  # inner (parallel) loop is each ensemble member point
  for (i in seq_len(nIter.transition)) {
    samplist <- foreach::foreach(s=samplist, .inorder=TRUE) %dopar% {
      # Initial values
      mcurr <- s$m
      ucurr <- s$u
      slcurr <- s$sl

      # m and u full conditional update
      tmp <- r_m_u_fc_smcmc(cmpdata, slcurr, priors$a, priors$b)
      mcurr <- tmp$m
      ucurr <- tmp$u

      # Links full conditional updates
      for (f in seq(2, length(files))) {
        if (proposals.transition == "LB") {
          slcurr <- draw.Z.locbal(f, cmpdata, slcurr, mcurr, ucurr,
                                  priors$aBM, priors$bBM,
                                  blocksize=blocksize)
        } else if (proposals.transition == "component") {
          slcurr <- draw.Z.componentwise(f, cmpdata, slcurr, mcurr, ucurr,
                                         priors$aBM, priors$bBM)
        }
      }

      ret <- list(m=mcurr, u=ucurr, sl=slcurr)

      # Only on last iteration: Calculate summary statistics of comparisons
      # for potential use later in PPRB (done here to take advantage of
      # parallelization)
      if (i == nIter.transition) {
        tmp <- disag.counts.allfiles(cmpdata, slcurr)

        ret$m.fc.pars <- tmp$match
        ret$u.fc.pars <- tmp$nonmatch
      }

      # Return a list of the ending values
      ret
    }
  }

  samplingend <- Sys.time() # Stop the clock

  if (cores > 1) {
    parallel::stopCluster(cl)
  }

  # Pack parallel results into arrays
  msamples <- matrix(NA, nrow=nrow(ensemble$m), ncol=ensemblesize)
  usamples <- matrix(NA, nrow=nrow(ensemble$u), ncol=ensemblesize)
  Zsamples <- matrix(NA, nrow=length(savestate(samplist[[1]]$sl)), ncol=ensemblesize)
  m.fc.pars <- matrix(0, nrow=nrow(ensemble$m), ncol=ensemblesize)
  u.fc.pars <- matrix(0, nrow=nrow(ensemble$u), ncol=ensemblesize)
  for (s in seq_len(ensemblesize)) {
    msamples[,s] <- samplist[[s]]$m
    usamples[,s] <- samplist[[s]]$u
    Zsamples[,s] <- savestate(samplist[[s]]$sl)
    m.fc.pars[,s] <- samplist[[s]]$m.fc.pars
    u.fc.pars[,s] <- samplist[[s]]$u.fc.pars
  }

  # Construct and return the new link state
  structure(
    list(
      Z = Zsamples,
      m = msamples,
      u = usamples,
      m.fc.pars = m.fc.pars,
      u.fc.pars = u.fc.pars,
      diagnostics = list(
        samplingtime = as.double(samplingend - samplingstart, units="secs"),
        jumpingtime = as.double(jumpingend - samplingstart, units="secs"),
        transitiontime = as.double(samplingend - jumpingend, units="secs")
      )
    ),
    class = "bstrlstate"
  )
}
