# Functions in this file are for proposing values of Z2 during the tripartite
# MCMC. Each function does a different method: global (i.e. prior), local
# uninformed, and local informed.


# Function to draw Z2 from the predictive distribution, at a given iteration
draw.Z2.global <- function(n1, n2, n3, Z, aBM, bBM) {
  # Number of links comes from beta-binomial
  nlinks <- rbinom(n=1, size=n3, prob=rbeta(n=1, aBM, bBM))
  # Given the number of links, which records in file 3 will be linked?
  link.from <- sample.int(n3, nlinks)
  # Candidates are any unlinked entries in file 1, plus all entries in file 2
  cand <- c(setdiff(seq_len(n1), Z), n1 + seq_len(n2))
  # Randomly sample candidates to be linked
  link.to <- cand[sample(length(cand), nlinks)]
  # Construct Z2 initially as totally unlinked, then fill links as appropriate
  Z2 <- n1 + n2 + seq_len(n3)
  Z2[link.from] <- link.to
  Z2
}

# Draws Z2 as a local uninformed step based on an add/delete/swap move
# Results in a symmetric proposal distribution
draw.Z2.local <- function(n1, n2, n3, Z, Z2.curr) {
  # Candidates are any unlinked entries in file 1, plus all entries in file 2
  cand <- c(setdiff(seq_len(n1), Z), n1 + seq_len(n2))
  # Pick an i and a j for the proposed move
  i <- sample(n3, size=1)
  j <- cand[sample(length(cand), 1)]
  # Return the new Z2
  return(perform.Z2.step(n1, n2, Z2.curr, i, j)$Z2)
}

# Draws a proposal for Z2 using a locally-balanced pointwise-informed proposal
# distribution based on Zanella (2020). The uninformed kernel K is a proposal
# over add/delete/swap steps which is uniform over pairs of records i,j. This
# is the proposal implemented here as draw.Z2.local(). The multiplicative term
# uses the Barker weights g(t) = t/(1+t). Returns a list with two elements:
#   Z2 = proposed value of Z2
#   mod = proposal component of M-H acceptance ratio. So that
#         alpha = pi(Z2new, ...)/pi(Z2curr, ...) * mod
draw.Z2.informed <- function(n1, n2, n3, Z, Z2.curr,
                             m, u, cmp.1to3, cmp.2to3, aBM, bBM, trace=FALSE,
                             blocksize=NULL) {
  # Candidates are any unlinked entries in file 1, plus all entries in file 2
  cand <- c(setdiff(seq_len(n1), Z), n1 + seq_len(n2))
  # Shrink possibilities down to blocksize
  checkblocks <- FALSE # If true, will need to check blocks for outside links
  if (is.null(blocksize) || (blocksize >= n3)) {
    # No reduction necessary in file 3
    iblock <- seq_len(n3)
  } else {
    # Need to select subset of file 3
    iblock <- sample(n3, size=blocksize, replace=FALSE)
    checkblocks <- TRUE
  }
  if (is.null(blocksize) || (blocksize >= length(cand))) {
    # No reduction necessary in candidates
    jblock <- cand
  } else {
    # Need to choose subset of candidates
    jblock <- cand[sample(length(cand), size=blocksize, replace=FALSE)]
    checkblocks <- TRUE
  }
  if (checkblocks) {
    # Check for links outside the selected blocks
    # Which records in the iblock are linked, and linked to records NOT IN the jblock?
    iblock.linked <- iblock[Z2.curr[iblock] <= n1 + n2]
    iblock.remove <- iblock.linked[!(Z2.curr[iblock.linked] %in% jblock)]
    iblock <- setdiff(iblock, iblock.remove)
    # Which records in the jblock are linked, and linked to records NOT IN the iblock?
    iblock.not <- setdiff(seq_len(n3), iblock)
    jblock <- setdiff(jblock, Z2.curr[iblock.not])
  }
  # What is the probability of making any given step?
  weights <- calc.Z2.stepmatrix(iblock, jblock, n1, n2, n3, m, u, Z, Z2.curr, cmp.1to3, cmp.2to3, aBM, bBM, trace=trace)
  # Sample i and j according to these
  i.draw <- iblock[sample(length(iblock), 1, prob=rowSums(weights))] # i marginally
  j.draw <- jblock[sample(length(jblock), 1, prob=weights[which(iblock == i.draw),])] # j conditionally
  tmp <- perform.Z2.step(n1, n2, Z2.curr, i.draw, j.draw)
  Z2.prop <- tmp$Z2
  reverse.move <- tmp$rev
  # What are the probabilities of backwards steps
  rev.weights <- calc.Z2.stepmatrix(iblock, jblock, n1, n2, n3, m, u, Z, Z2.prop, cmp.1to3, cmp.2to3, aBM, bBM, trace=trace)
  # Return the proposed value and MH acceptance ratio component
  mhmod <- rev.weights[which(iblock == reverse.move[1]), which(jblock == reverse.move[2])] / weights[which(iblock == i.draw),which(jblock == j.draw)]
  return(list(Z2=Z2.prop, mod=mhmod))
}


################################################################################
## Helper functions for Z2 proposals ###########################################
################################################################################


# Given a current state of Z2, a choice i from 1 to n3 (i.e. an element of
# file 3), and a choice j between 1 and n1+n2 (i.e. an element of the
# candidate set), determine the appropriate add/delete/single swap/double swap.
# move to make and make it. Returns the new value of Z2 after the move has
# been made, and a pair i,j that can be used to make the reverse move from the
# new Z2 back to Z2.curr.
# Return format: list(Z2=<new state>, rev=c(i,j))
perform.Z2.step <- function(n1, n2, Z2.curr, i, j) {
  # Which nodes are i and j linked to, if any?
  i.linkedto <- Z2.curr[i] # An element of cand, or n1+n2+j
  j.linkedto <- which(Z2.curr == j) # An element of file 3, or 0
  if (length(j.linkedto) == 0) j.linkedto <- 0
  # Initialize Z2.prop as identical to current state
  Z2.prop <- Z2.curr
  # Initialize backwards move to be (i,j)
  reverse.move <- c(i,j)
  # Which move should be made? Make that move
  if (i.linkedto == j) { # delete
    Z2.prop[i] <- n1 + n2 + i
    reverse.move <- c(i,j) # Reverse = add
  } else if ((i.linkedto == n1+n2+i) && (j.linkedto == 0)) { # add
    Z2.prop[i] <- j
    reverse.move <- c(i,j) # Reverse = delete
  } else if ((i.linkedto == n1+n2+i) && (j.linkedto != 0)) { # single-swap 1
    Z2.prop[i] <- j
    Z2.prop[j.linkedto] <- n1 + n2 + j.linkedto
    reverse.move <- c(j.linkedto, j) # Reverse = single-swap back
  } else if ((i.linkedto <= n1+n2) && (j.linkedto == 0)) { # single-swap 2
    Z2.prop[i] <- j
    reverse.move <- c(i, i.linkedto) # Reverse = single-swap back
  } else if ((i.linkedto <= n1+n2) && (j.linkedto != 0)) { # double-swap
    Z2.prop[i] <- j
    Z2.prop[j.linkedto] <- i.linkedto
    reverse.move <- c(i, i.linkedto) # Reverse = double-swap back
  } else {
    # Always good practice to have an error "else"
    stop("Impossible condition in perform.Z2.step")
  }
  return(list(Z2=Z2.prop, rev=reverse.move))
}

# This function produces a matrix of probabilities for each step in an informed
# Z2 proposal.
# Returns: A matrix to hold, eventually, the proposal probability for a
# move i,j. Each row corresponds to an element of file 3. Each column
# corresponds to an element of the candidate set, in the order listed in `cand`.
# The value in M[i,j] is the probability of making the move i,j. Uses the
# Barker weights g(t) = t/(1+t)
# Parameters
# ivec - vector of indices of records in file 3 to consider for steps
# jvec - vector of indices of records in files 1/2 (candidates) to consider for steps
# (everything else) - needed to calculate likelihoods for informed step probs
calc.Z2.stepmatrix <- function(ivec, jvec,
                               n1, n2, n3, m, u, Z, Z2.curr, cmp.1to3, cmp.2to3,
                               aBM, bBM, trace=FALSE) {
  # Initialize matrix
  weights <- matrix(0, nrow=length(ivec), ncol=length(jvec))
  # Which likelihood function will we use?
  if (trace) {
    ell <- calc.log.lkl.tracing
  } else {
    ell <- calc.log.lkl
  }
  # Fill in the posterior at each resulting position
  for (fromidx in 1:length(ivec)) {
    for (candidx in 1:length(jvec)) {
      i <- ivec[fromidx]
      j <- jvec[candidx]
      # What Z2 would be proposed from modifying the pair (i,j)?
      Z2.prop <- perform.Z2.step(n1, n2, Z2.curr, i, j)$Z2
      # Calculate posterior for new Z2
      weights[fromidx, candidx] <- (ell(cmp.1to3, cmp.2to3, n1, n2, n3, m, u, Z, Z2.prop) + calc.log.Z2prior(n1, n2, n3, Z2.prop, Z, aBM, bBM))
    }
  }
  # Subtract the log posterior at the current state
  weights <- weights - (ell(cmp.1to3, cmp.2to3, n1, n2, n3, m, u, Z, Z2.curr) + calc.log.Z2prior(n1, n2, n3, Z2.curr, Z, aBM, bBM))
  # Bring out of log scale
  weights <- exp(weights)
  # Calculate g(t) for all elements of the matrix
  weights <- weights / (1 + weights)
  # Divide by the normalizing constant
  return(weights / sum(weights))
}
