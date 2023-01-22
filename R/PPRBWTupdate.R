# File for containing PPRB-with-transition update, which combines benefits of
# PPRB and SMCMC updates

#' Perform a PPRB-with-Transition update of record linkage with a new file
#' @param state Existing record linkage state. Returned by either bipartiteRL,
#'   PPRBupdate, SMCMCupdate, PRPBWTupdate, or multifileRL. This is used as the
#'   ensemble of samples for the update
#' @param newfile A data.frame representing the new file: one row per record
#' @param flds Names of fields in the new file to use for comparison. Only used
#'   if no global field names were specified in bipartiteRL initially.
#' @param nIter.PPRB Number of iterations for which to run the PPRB sampler. By
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
#' @param refresh How often to output an update including the iteration number
#'   and percent complete. Outputs only during the PPRB application. If
#'   refresh >= 1, taken as a number of iterations
#'   between messages (rounded). If 0 < refresh < 1, taken as the proportion of
#'   nIter. If refresh == 0, no messages are displayed.
#' @param nIter.transition Number of iterations to use in the transition kernel
#'   following the PPRB update for each ensemble sample.
#' @param proposals.transition Which kernel to use for Z
#'   updates in the jumping and transition kernels, respectively.
#' @param cores Number of cores to use for parallel execution. If cores == 1,
#'   update is run sequentially. A cluster is created using
#'   parallel::makeCluster().
#' @param seed Random seed to set at the beginning of the MCMC run
#'
#' @return An object of class 'bstrlstate' containing posterior samples and
#'   necessary metadata for passing to future streaming updates.
#' @examples
#' data(geco_small)
#' data(geco_small_result)
#'
#' # Add fifth file to previous four-file link result
#' filtered <- thinsamples(geco_small_result, 2) # Filter ensemble to 2 - very small for example
#' file5.result <- PPRBWTupdate(filtered, geco_small[[5]],
#'                             nIter.PPRB=2, burn=1, nIter.transition=1, # Very small run for example
#'                             proposals.transition="LB",
#'                             blocksize=5)
#' @export
PPRBWTupdate <- function(state, newfile, flds = NULL,
                         nIter.PPRB = NULL, burn = 0, blocksize = NULL,
                         threestep = TRUE, refresh=0.1,
                         nIter.transition=10, proposals.transition=c("LB", "component"),
                         cores=1,
                         seed=0) {

  proposals.transition <- match.arg(proposals.transition)

  # First, we can outsource the PPRB update to the existing PPRB update function
  postpprb <- PPRBupdate(state, newfile, flds, nIter=nIter.PPRB, burn=burn, blocksize=blocksize,
                         threestep=threestep, seed=seed, refresh=refresh)

  # SMCMC transition kernel application starts here
  ensemble <- list(m=postpprb$m, u=postpprb$u, Z=postpprb$Z)

  updated <- coreSMCMCupdate(ensemble, postpprb$priors, postpprb$files, postpprb$comparisons,
                             nIter.jumping=0, nIter.transition=nIter.transition,
                             cores=cores,
                             proposals.jumping="LB", proposals.transition=proposals.transition,
                             blocksize=blocksize, seed=NULL)

  # Assemble final object and return
  updated$files <- postpprb$files
  updated$comparisons <- postpprb$comparisons
  updated$priors <- postpprb$priors
  updated$cmpdetails <- postpprb$cmpdetails
  # Combine timing/diagnostics
  updated$diagnostics <- list(
    pprb.burntime = postpprb$diagnostics$burntime,
    pprb.samplingtime = postpprb$diagnostics$samplingtime,
    pprb.accepted = postpprb$diagnostics$pprb.accepted,
    transitiontime = updated$diagnostics$transitiontime,
    itertimes = updated$diagnostics$itertimes
  )

  return(updated)

}
