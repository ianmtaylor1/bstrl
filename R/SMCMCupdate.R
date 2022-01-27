
#' @export
SMCMCupdate <- function(state, newfile, flds=NULL, nIter.jumping, nIter.transition,
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


  # TODO: the actual SMCMC update


  # Construct and return the new link state
  structure(
    list(
      Z = NULL,
      m = NULL,
      u = NULL,
      files = files,
      comparisons = cmpdata,
      priors = state$priors,
      cmpdetails = cmpdetails,
    ),
    class = "bstrlstate"
  )
}
