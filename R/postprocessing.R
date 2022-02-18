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
  if ((count > 1) && (count <= N)) {
    every <- (N - 1) %/% (count - 1)
    iterfilter <- seq(N - (count - 1) * every, N, by=every)
  } else if (count > N) {
    warning("Cannot filter to greater number of samples")
    iterfilter <- seq_len(N)
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

#' Extract the links from a bstrlstate object into a list of streaminglinks
#' objects.
#'
#' @param state A bstrlstate object returned by one of the main RL functions.
#'
#' @return A list of streaminglinks objects, one per posterior sample contained
#'   in 'state'.
#'
#' @export
extractlinks <- function(state) {
  `%do%` <- foreach::`%do%`
  filesizes <- rep(NA, length(state$files))
  for (f in seq_along(state$files)) {
    filesizes[f] <- nrow(state$files[[f]])
  }
  foreach::foreach(i = seq_len(ncol(state$Z))) %do% {
    streaminglinks(filesizes, state$Z[,i])
  }
}

#' Calculate the precision of estimated links relative to true links
#'
#' @param sl.est streaminglinks object representing link estimates
#' @param sl.true streaminglinks object representing true links
#'
#' @return The precision of the estimated links.
#'
#' @export
precision <- function(sl.est, sl.true) {
  estlinks <- alllinks(sl.est, idx="global")
  istrue <- mapply(function(x, y, sl) islinked.gl(sl, x, y),
                   estlinks$idx1, estlinks$idx2,
                   MoreArgs = list(sl.true),
                   SIMPLIFY = TRUE)
  mean(istrue)
}


#' Calculate the recall of estimated links relative to true links
#'
#' @param sl.est streaminglinks object representing link estimates
#' @param sl.true streaminglinks object representing true links
#'
#' @return The recall of the estimated links.
#'
#' @export
recall <- function(sl.est, sl.true) {
  truelinks <- alllinks(sl.true, idx="global")
  isest <- mapply(function(x, y, sl) islinked.gl(sl, x, y),
                   truelinks$idx1, truelinks$idx2,
                   MoreArgs = list(sl.est),
                   SIMPLIFY = TRUE)
  mean(isest)
}
