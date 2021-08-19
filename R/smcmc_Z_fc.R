# This file contains componentwise full conditional samplers for Z and Z2
# 3/23/2021 - Currently used with the SMCMC transition and approximate
# jumping kernels. To do: make flexible enough to deal with more than 3 files


# Z, Z2, m, u - current parameter values
# cmpdata - a list-of-lists format of comparison data
# aBM, bBM - Z/Z2 prior parameters.
# directratio - whether to calculate the likelihood portion of the full conditional
#     directly using likelihood ratios (i.e. using the unlinked state as reference
#     value). If FALSE, calculate the full likelihood for each possible state (i.e.
#     the "slow" evaluation.)
#
# Returns: the new value of Z
r_Z_fc_smcmc <- function(Z, Z2, m, u, cmpdata, aBM, bBM, directratio=TRUE, Z2prior=c("default", "flat", "noinv")) {
  n1 <- cmpdata[[1]][[1]]$n1

  # Parse Z2 prior choice
  Z2prior <- match.arg(Z2prior)

  for (j in seq_along(Z)) {
    weights <- rep(-Inf, n1+1)
    for (i in seq_len(n1+1)) {
      if (i > n1) {
        rec.i <- n1+j
      } else {
        rec.i <- i
      }
      Z.prop <- Z
      Z.prop[j] <- rec.i
      # First, the priors of Z and Z2
      weights[i] <- (
        calc.log.Z2prior(n1, Z.prop, c(), aBM, bBM) +  ## Z1 - same as before
        if (!valid.link.state(n1, Z.prop, Z2)) {       ## Z2 - multiple possibilities
          # Need to check for validity here because now not all priors do it
          -Inf
        } else if (Z2prior == "default") {
          calc.log.Z2prior(n1, Z2, Z.prop, aBM, bBM)
        } else if (Z2prior == "flat") {
          calc.log.Z2prior.flat(n1, Z2, Z.prop)
        } else if (Z2prior == "noinv") {
          calc.log.Z2prior.noinvalid(n1, Z2, Z.prop, aBM, bBM)
        } else {
          stop("Invalid Z2prior value. Should not be here, match.arg() must not have worked.")
        }
      )
      # If the prior has nonzero probability, calculate likelihoods
      # Since anything + -Inf == -Inf in R, this just saves time.
      if (weights[i] > -Inf) {
        weights[i] <- weights[i] + if (directratio) {
          smcmc.log.lkl.ratio.Z(cmpdata, m, u, Z.prop, Z2, j)
        } else {
          calc.log.lkl.tracing(cmpdata[[1]], m, u, c(), Z.prop) + calc.log.lkl.tracing(cmpdata[[2]], m, u, Z.prop, Z2)
        }
      }
    }
    # Semi-normalize, just to avoid too much over/underflow in exp
    weights <- weights - max(weights)
    # Select new value, put it directly in Z
    Z[j] <- sample(c(seq_len(n1), n1+j), 1, prob=exp(weights))
  }

  return(Z)
}

# Z, Z2, m, u - current parameter values
# cmpdata - a list-of-lists format of comparison data
# aBM, bBM - Z/Z2 prior parameters.
# directratio - whether to calculate the likelihood portion of the full conditional
#     directly using likelihood ratios (i.e. using the unlinked state as reference
#     value). If FALSE, calculate the full likelihood for each possible state (i.e.
#     the "slow" evaluation.)
#
# Returns: the new value of Z2
r_Z2_fc_smcmc <- function(Z, Z2, m, u, cmpdata, aBM, bBM, directratio=TRUE, Z2prior=c("default", "flat", "noinv")) {
  n1 <- cmpdata[[1]][[1]]$n1
  n2 <- cmpdata[[1]][[1]]$n2

  # Parse Z2 prior choice
  Z2prior <- match.arg(Z2prior)

  for (j in seq_along(Z2)) {
    weights <- rep(-Inf, n1+n2+1)
    # All n1 + n2 linked possibilities
    for (i in seq_len(n1+n2+1)) {
      if (i > n1 + n2) {
        rec.i <- n1+n2+j
      } else {
        rec.i <- i
      }
      Z2.prop <- Z2
      Z2.prop[j] <- rec.i
      # First, the priors of Z and Z2
      weights[i] <- if (!valid.link.state(n1, Z, Z2.prop)) {
        # Need to do validity checking here since not all priors do it anymore
        -Inf
      } else if (Z2prior == "default") {
        calc.log.Z2prior(n1, Z2.prop, Z, aBM, bBM)
      } else if (Z2prior == "flat") {
        calc.log.Z2prior.flat(n1, Z2.prop, Z)
      } else if (Z2prior == "noinv") {
        calc.log.Z2prior.noinvalid(n1, Z2.prop, Z, aBM, bBM)
      } else {
        stop("Invalid Z2prior value. Should not be here, match.arg() must not have worked.")
      }
      # If the prior has nonzero probability, calculate likelihoods
      # Since anything + -Inf == -Inf in R, this just saves time.
      if (weights[i] > -Inf) {
        weights[i] <- weights[i] + if (directratio) {
          smcmc.log.lkl.ratio.Z2(cmpdata, m, u, Z, Z2.prop, j)
        } else {
          calc.log.lkl.tracing(cmpdata[[2]], m, u, Z, Z2.prop)
        }
      }
    }
    # Semi-normalize, just to avoid too much over/underflow in exp
    weights <- weights - max(weights)
    # Select new value, put it directly in Z
    Z2[j] <- sample(c(seq_len(n1+n2), n1+n2+j), 1, prob=exp(weights))
  }

  return(Z2)
}



################################################################################
## Z and Z2 full conditionals based on Zanella (2020), but adapted for SMCMC
## or Gibbs sampling
################################################################################


# Propose, then accept or reject, a value of Z using a locally balanced proposal
# of the form of Zanella (2020). Takes the same arguments as the Sadinle (2017)
# full conditional functions except that it has no option to do directratio
# likelihood computation. Instead it takes a blocksize parameter.
r_Z_fc_smcmc_zanella <- function(Z, Z2, m, u, cmpdata, aBM, bBM, Z2prior=c("default", "flat", "noinv"), blocksize=NULL) {
  # Parse Z2 prior choice
  Z2prior <- match.arg(Z2prior)

}

r_Z2_fc_smcmc_zanella <- function(Z, Z2, m, u, cmpdata, aBM, bBM, Z2prior=c("default", "flat", "noinv"), blocksize=NULL) {
  # Parse Z2 prior choice
  Z2prior <- match.arg(Z2prior)

  # File sizes
  n1 <- cmpdata[[1]][[1]]$n1
  n2 <- cmpdata[[1]][[1]]$n2
  n3 <- cmpdata[[2]][[2]]$n2

  Z2.curr <- Z2 # To differentiate from proposed value

  # Create blocks of possibilities for movesShrink possibilities down to blocksize
  blocks <- create.blocks(blocksize, n1+n2, seq_len(n1+n2), Z2.curr)
  while ((length(blocks$iblock) == 0) || (length(blocks$jblock) == 0)) {
    blocks <- create.blocks(blocksize, n1+n2, seq_len(n1+n2), Z2.curr)
  }
  iblock <- blocks$iblock
  jblock <- blocks$jblock

  # What is the probability of making any given step?
  weights <- calc.Z2.stepmatrix.smcmc(iblock, jblock, cmpdata, m, u, Z, Z2.curr, aBM, bBM, Z2prior)

  # Sample i and j according to these
  i.draw <- iblock[sample(length(iblock), 1, prob=rowSums(weights))] # i marginally
  j.draw <- jblock[sample(length(jblock), 1, prob=weights[which(iblock == i.draw),])] # j conditionally
  tmp <- perform.Z2.step(nprev, Z2.curr, i.draw, j.draw) # This references
  Z2.prop <- tmp$Z2
  reverse.move <- tmp$rev

  # What are the probabilities of backwards steps
  rev.weights <- calc.Z2.stepmatrix.smcmc(iblock, jblock, cmpdata, m, u, Z, Z2.prop, aBM, bBM, Z2prior)

  # Return the proposed value and MH acceptance ratio component
  mhmod <- rev.weights[which(iblock == reverse.move[1]), which(jblock == reverse.move[2])] / weights[which(iblock == i.draw),which(jblock == j.draw)]

  # Metropolis-hastings step and return
  log.alpha <- (
    calc.log.lkl.tracing(cmpdata[[2]], m, u, Z, Z2.prop)
    - calc.log.lkl.tracing(cmpdata[[2]], m, u, Z, Z2.curr)
    + if (Z2prior == "default") {
        calc.log.Z2prior(n1, Z2.prop, Z, aBM, bBM)
      } else if (Z2prior == "flat") {
        calc.log.Z2prior.flat(n1, Z2.prop, Z)
      } else if (Z2prior == "noinv") {
        calc.log.Z2prior.noinvalid(n1, Z2.prop, Z, aBM, bBM)
      } else {
        stop("Invalid Z2prior value. Should not be here, match.arg() must not have worked.")
      }
    - if (Z2prior == "default") {
        calc.log.Z2prior(n1, Z2.curr, Z, aBM, bBM)
      } else if (Z2prior == "flat") {
        calc.log.Z2prior.flat(n1, Z2.curr, Z)
      } else if (Z2prior == "noinv") {
        calc.log.Z2prior.noinvalid(n1, Z2.curr, Z, aBM, bBM)
      } else {
        stop("Invalid Z2prior value. Should not be here, match.arg() must not have worked.")
      }
    + log(mhmod)
  )
  if (log(runif(1)) < log.alpha) {
    return(Z2.prop)
  } else {
    return(Z2.curr)
  }
}


# This function produces a matrix of probabilities for each step in an informed
# Z2 proposal. It is designed to work in the context of Zanella LB proposals within
# an SMCMC context.
# Returns: A matrix to hold, eventually, the proposal probability for a
# move i,j. Each row corresponds to an element of file 3. Each column
# corresponds to an element of the candidate set, in the order listed in `cand`.
# The value in M[i,j] is the probability of making the move i,j. Uses the
# Barker weights g(t) = t/(1+t)
# Parameters
# ivec - vector of indices of records in file 3 to consider for steps
# jvec - vector of indices of records in previous files (candidates) to consider for steps
# (everything else) - needed to calculate likelihoods for informed step probs
calc.Z2.stepmatrix.smcmc <- function(ivec, jvec,
                                     cmpdata, m, u, Z, Z2.curr, aBM, bBM,
                                     Z2prior) {
  # Required filesizes
  n1 <- cmpdata[[1]][[1]]$n1
  nprev <- 0
  for (file in seq_len(length(cmpdata))) {
    nprev <- nprev + cmpdata[[file]][[file]]$n1
  }
  # Initialize matrix
  weights <- matrix(0, nrow=length(ivec), ncol=length(jvec))
  # Fill in the posterior at each resulting position
  for (fromidx in 1:length(ivec)) {
    for (candidx in 1:length(jvec)) {
      i <- ivec[fromidx]
      j <- jvec[candidx]
      # What Z2 would be proposed from modifying the pair (i,j)?
      Z2.prop <- perform.Z2.step(nprev, Z2.curr, i, j)$Z2
      # Calculate posterior for new Z2
      weights[fromidx, candidx] <- if (!valid.link.state(n1, Z, Z2.prop)) {
        # Need to do validity checking here since not all priors do it anymore
        -Inf
      } else if (Z2prior == "default") {
        calc.log.Z2prior(n1, Z2.prop, Z, aBM, bBM)
      } else if (Z2prior == "flat") {
        calc.log.Z2prior.flat(n1, Z2.prop, Z)
      } else if (Z2prior == "noinv") {
        calc.log.Z2prior.noinvalid(n1, Z2.prop, Z, aBM, bBM)
      } else {
        stop("Invalid Z2prior value. Should not be here, match.arg() must not have worked.")
      }
      if (weights[fromidx, candidx] > -Inf) {
        weights[fromidx, candidx] <- weights[fromidx, candidx] + calc.log.lkl.tracing(cmpdata[[2]], m, u, Z, Z2.prop)
      }
    }
  }
  # Subtract the log posterior at the current state
  weights <- weights - (
    calc.log.lkl.tracing(cmpdata[[2]], m, u, Z, Z2.curr)
    + if (Z2prior == "default") {
      calc.log.Z2prior(n1, Z2.curr, Z, aBM, bBM)
    } else if (Z2prior == "flat") {
      calc.log.Z2prior.flat(n1, Z2.curr, Z)
    } else if (Z2prior == "noinv") {
      calc.log.Z2prior.noinvalid(n1, Z2.curr, Z, aBM, bBM)
    } else {
      stop("Invalid Z2prior value. Should not be here, match.arg() must not have worked.")
    }
  )
  # Bring out of log scale
  weights <- exp(weights)
  # Calculate g(t) for all elements of the matrix
  weights <- weights / (1 + weights)
  # Divide by the normalizing constant
  return(weights / sum(weights))
}
