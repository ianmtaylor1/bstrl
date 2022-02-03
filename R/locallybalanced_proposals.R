# This file contains functions to do locally balanced proposal updates for links

#################################
# Last file (i.e. most recent Z vector)
#################################


# Perform a locally balanced proposal step for the last file's links.
# Parameters:
#   lastfilecmps = list of k-1 comparison data objects, comparing each of the first
#     k-1 files with the kth file. This is the last element in the triangular
#     list-of-lists format for all comparison data.
#   m,u = current parameter values
#   sl = streaminglinks object defining current link state
#   aBM,bBM = prior parameters for link vectors
#   blocksize = blocksize for reduced complexity
# Returns:
#   a new state of streaming link object: different by a local step if proposal
#   was accepted, the same if proposal was rejected
draw.Z.locbal.lastfile <- function(lastfilecmps, sl, m, u, aBM, bBM, blocksize=NULL) {
  # Create blocks
  blocks <- createblocks(sl, file = nfiles(sl), blocksize = blocksize)
  iblock <- blocks$iblock
  jblock <- blocks$jblock

  # Calculate proposal probabilities
  probs <- calc.locbal.probs.lastfile(iblock, jblock, lastfilecmps, sl, m, u, aBM, bBM)

  # Pick proposal
  i.index <- sample(length(iblock), 1, prob=rowSums(probs))
  i <- iblock[i.index]
  j.index <- sample(length(jblock), 1, prob=probs[i.index,])
  j <- jblock[j.index]
  tmp <- performstep(sl, i, j)
  slprop <- tmp$state
  reverse <- tmp$reverse
  revi <- reverse[1]
  revj <- reverse[2]
  revi.index <- which(iblock == revi)
  revj.index <- which(jblock == revj)

  # Calculate back-step probabilities
  reverseprobs <- calc.locbal.probs.lastfile(iblock, jblock, lastfilecmps, slprop, m, u, aBM, bBM)

  # M-H acceptance step
  log.alpha <- (
    calc.log.lkl.lastfile(lastfilecmps, m, u, slprop) -
    calc.log.lkl.lastfile(lastfilecmps, m, u, sl) +
    log.Zprior(slprop, aBM, bBM, vec="last") -
    log.Zprior(sl, aBM, bBM, vec="last") +
    log(reverseprobs[revi.index, revj.index]) -
    log(probs[i.index, j.index])
  )
  if (log(runif(1)) < log.alpha) {
    return(slprop)
  } else {
    return(sl)
  }
}


# Calculate a matrix with proposal probabilities according to the locally
# balanced proposal weighting
calc.locbal.probs.lastfile <- function(iblock, jblock, lastfilecmps, sl, m, u, aBM, bBM) {
  logprobs <- matrix(0, nrow=length(iblock), ncol=length(jblock))

  # Calculate (value proportional to) log posterior for each step and normalize
  # by the log posterior at the current state
  for (i.index in seq_along(iblock)) {
    for (j.index in seq_along(jblock)) {
      i <- iblock[i.index]
      j <- jblock[j.index]
      slprop <- performstep(sl, i, j)$state
      logprobs[i.index, j.index] <- (
        calc.log.lkl.lastfile(lastfilecmps, m, u, slprop) +
        log.Zprior(slprop, aBM, bBM, vec="last")
      )
    }
  }
  logprobs <- logprobs - (
    calc.log.lkl.lastfile(lastfilecmps, m, u, sl) +
    log.Zprior(sl, aBM, bBM, vec="last")
  )

  # Bring out of log scale, use Barker weights function, normalize
  probs <- exp(logprobs)
  probs <- probs / (1 + probs)
  return(probs / sum(probs))
}


#######################################################################
# Arbitrary file
###################


# Perform a locally balanced proposal step for the last file's links.
# Parameters:
#   cmpdata - a list of lists of comparison data objects. (Standard triangular
#             comparison data format.)
#   m,u = current parameter values
#   sl = streaminglinks object defining current link state
#   aBM,bBM = prior parameters for link vectors
#   blocksize = blocksize for reduced complexity
# Returns:
#   a new state of streaming link object: different by a local step if proposal
#   was accepted, the same if proposal was rejected
draw.Z.locbal <- function(file, cmpdata, sl, m, u, aBM, bBM, blocksize=NULL) {
  # Create blocks
  blocks <- createblocks(sl, file = file, blocksize = blocksize)
  iblock <- blocks$iblock
  jblock <- blocks$jblock

  # Calculate proposal probabilities
  probs <- calc.locbal.probs(iblock, jblock, cmpdata, sl, m, u, aBM, bBM, file)

  # Pick proposal
  i.index <- sample(length(iblock), 1, prob=rowSums(probs))
  i <- iblock[i.index]
  j.index <- sample(length(jblock), 1, prob=probs[i.index,])
  j <- jblock[j.index]
  tmp <- performstep(sl, i, j)
  slprop <- tmp$state
  reverse <- tmp$reverse
  revi <- reverse[1]
  revj <- reverse[2]
  revi.index <- which(iblock == revi)
  revj.index <- which(jblock == revj)

  # Calculate back-step probabilities
  reverseprobs <- calc.locbal.probs(iblock, jblock, cmpdata, slprop, m, u, aBM, bBM, file)

  # M-H acceptance step
  log.alpha <- (
    calc.log.lkl.tailfiles(cmpdata, m, u, slprop, file) -
      calc.log.lkl.tailfiles(cmpdata, m, u, sl, file) +
      log.Zprior(slprop, aBM, bBM, vec="all") -
      log.Zprior(sl, aBM, bBM, vec="all") +
      log(reverseprobs[revi.index, revj.index]) -
      log(probs[i.index, j.index])
  )
  if (log(runif(1)) < log.alpha) {
    return(slprop)
  } else {
    return(sl)
  }
}


# Calculate a matrix with proposal probabilities according to the locally
# balanced proposal weighting
# This function is used with the locally balanced update for an arbitrary Z
# vector, not just the last Z vector.
calc.locbal.probs <- function(iblock, jblock, cmpdata, sl, m, u, aBM, bBM, startfile) {
  logprobs <- matrix(0, nrow=length(iblock), ncol=length(jblock))

  # Calculate (value proportional to) log posterior for each step and normalize
  # by the log posterior at the current state
  for (i.index in seq_along(iblock)) {
    for (j.index in seq_along(jblock)) {
      i <- iblock[i.index]
      j <- jblock[j.index]
      slprop <- performstep(sl, i, j)$state
      logprobs[i.index, j.index] <- (
        calc.log.lkl.tailfiles(cmpdata, m, u, slprop, startfile) +
          log.Zprior(slprop, aBM, bBM, vec="all")
      )
    }
  }
  logprobs <- logprobs - (
    calc.log.lkl.tailfiles(cmpdata, m, u, sl, startfile) +
      log.Zprior(sl, aBM, bBM, vec="all")
  )

  # Bring out of log scale, use Barker weights function, normalize
  probs <- exp(logprobs)
  probs <- probs / (1 + probs)
  return(probs / sum(probs))
}
