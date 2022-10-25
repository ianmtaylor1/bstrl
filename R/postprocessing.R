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
#' @details
#' The function does one of three things when passed an unfinished
#' sampling state, e.g. from multifileRL after time limit expired. For
#' 'ignore', the desired burn is performed and any remaining samples are
#' discarded, returning only the number of completed post-burn samples. If
#' the link state does not have any completed post-burn samples, an empty
#' list is returned. If 'warn' (default), the same action is performed as
#' 'ignore' but a warning is issued. If 'fail', any unfinished link state will
#' cause the function to fail.
#'
#' @param state A bstrlstate object returned by one of the main RL functions.
#' @param unfinished What to do if passed an unfinished sampling state, e.g.
#'   from multifileRL. See details.
#'
#' @return A list of streaminglinks objects, one per posterior sample contained
#'   in 'state'.
#'
#' @export
extractlinks <- function(state, unfinished=c("warn", "ignore", "fail")) {
  unfinished <- match.arg(unfinished)
  `%do%` <- foreach::`%do%`

  # Handle unfinished sampling
  if (("continue" %in% names(state)) && !is.null(state$continue)) {
    if (unfinished == "fail") {
      stop("Cannot extract links from unfinished sampling")
    }
    burn <- state$continue$needtoburn
    maxiter <- state$continue$maxiter
    if (unfinished == "warn") {
      totaliter <- ncol(state$Z)
      warning("Unfinished sampling, extracting ", max(maxiter-burn, 0), " of ", totaliter-burn, " samples.")
    }
  } else {
    burn <- 0
    maxiter <- ncol(state$Z)
  }
  if (maxiter <= burn) {
    return(list())
  }

  filesizes <- rep(NA, length(state$files))
  for (f in seq_along(state$files)) {
    filesizes[f] <- nrow(state$files[[f]])
  }

  foreach::foreach(i = seq(burn + 1, maxiter)) %do% {
    i <- get("i") # Terrible hack to suppress CRAN check notes
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
