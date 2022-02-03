
# Perform a full conditional update of a Z vector by drawing a new value for
# each of its components componentwise.
# Parameters:
#   file = The file number whose Z vector we wish to update.
#   cmpdata - a list of lists of comparison data objects. (Standard triangular
#             comparison data format.)
#   m,u = current parameter values
#   sl = streaminglinks object defining current link state
#   aBM,bBM = prior parameters for link vectors
# Returns:
#   a new state of streaming link object: potentially different in every single
# component for this file
draw.Z.componentwise <- function(file, cmpdata, sl, m, u, aBM, bBM) {
  stopifnot(file >= 2, file <= nfiles(sl))

  # Number of possibilities for each component to have
  possibilities <- sum(sl$ns[seq_len(file-1)])

  # Loop through doing the componentwise full conditional for all records in
  # the chosen file
  for (rec in seq_len(sl$ns[file])) {
    globalrec <- local.to.global(sl$ns, file, rec)
    slbase <- unlink.down.gl(sl, globalrec)

    logprobs <- rep(-Inf, possibilities + 1)

    # Fill in the probability value for the unlinked state
    logprobs[possibilities + 1] <- (
      calc.log.lkl.tailfiles(cmpdata, m, u, slbase, file) +
      log.Zprior(slbase, aBM, bBM, vec="all")
    )

    # Fill in the probability values for the linked states
    for (dest in seq_len(possibilities)) {
      slprop <- add.link.gl(slbase, dest, globalrec, conflict="null")
      if (!is.null(slprop)) {
        logprobs[dest] <- (
          calc.log.lkl.tailfiles(cmpdata, m, u, slprop, file) +
          log.Zprior(slprop, aBM, bBM, vec="all")
        )
      }
    }

    # Normalize, exponentiate probabilities and draw new state
    probs <- exp(logprobs - max(logprobs))
    dest <- sample(possibilities + 1, 1, prob=probs)
    if (dest > possibilities) {
      sl <- slbase
    } else {
      # If it had non-zero probability, it *should* be a valid link
      sl <- add.link.gl(slbase, dest, globalrec, conflict="error")
    }
  }

  sl
}
