# This file contains functions for post-processing objects of class bstrlstate
# (the output of sampling functions). These functions can thin samples
# to a more reasonable number, calculcate posterior linkage probabilities or
# other posterior quantities, etc.

#' Thin a bstrlstate object
#'
#' Thin a bstrlstate object by keeping only one sample every n, until the
#' desired count remains.
#'
#' @details This is useful to do before an SMCMC update. SMCMC produces
#' independent samples, so fewer are required to get the same quality estimates.
#'
#' @param state Object of class bstrlstate, output by bipartiteRL, SMCMCupdate,
#'   PPRBupdate, or multifileRL
#' @param count The number of desired samples after filtering
#'
#' @return An object of class bstrlstate, containing count samples.
#'
#' @export
thinsamples <- function(state, count) {
  N <- ncol(state$Z)
  if (count > 1) {
    every <- (N - 1) %/% (count - 1)
    iterfilter <- seq(N - (count - 1) * every, N, by=every)
  } else {
    iterfilter <- c(N)
  }

  structure(
    list(
      m = state$m[,iterfilter,drop=F],
      u = state$u[,iterfilter,drop=F],
      Z = state$Z[,iterfilter,drop=F],
      files = state$files,
      comparisons = state$comparisons,
      cmpdetails = state$cmpdetails,
      priors = state$priors,
      m.fc.pars = state$m.fc.pars[,iterfilter,drop=F],
      u.fc.pars = state$u.fc.pars[,iterfilter,drop=F]
    ),
    class = "bstrlstate"
  )
}
