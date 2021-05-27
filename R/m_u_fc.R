# This file contains functions related to the full conditional distribution
# of m and u.

# It depends on the package extraDistr for calculating the density and sampling
# from Dirichlet distributions

# Draw m and u from their full conditional distributions.
# Parameters:
#   cmpdata - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#   Z, Z2 - current values of the parameters Z and Z2
#   a, b - prior hyperparameters for the distributions of m and u, respectively.
#          See ?BRL::BRL for explanation
#   m.prev.pars - additions to hyperparameters for m from matches between
#          previous files. These are computed in a post-processing step of the
#          samples of Z, and can just be passed directly here.
#   u.prev.pars - same as m.prev.pars, but for u.
#   trace - whether to do link tracing.
# Return:
#   A list with two elements, m and u, which contain the fc sampled m and u.
r_m_u_fc <- function(cmpdata, Z, Z2,
                     a, b, m.prev.pars, u.prev.pars, trace=FALSE) {
  # 1. Create the additions to m and u hyperparameters derived from links with
  # the newest (third) file.
  counts <- disag.counts(cmpdata, Z, Z2, do.trace=trace)
  m.new.pars <- counts$match
  u.new.pars <- counts$nonmatch
  # 2. combine a, m.prev.pars and new m pars to sample m
  m <- rdirichlet.multi(alpha=m.new.pars + m.prev.pars + a, groups=cmpdata[[1]]$nDisagLevs)
  # 3. combine b, u.prev.pars and new u pars to sample u
  u <- rdirichlet.multi(alpha=u.new.pars + u.prev.pars + b, groups=cmpdata[[1]]$nDisagLevs)
  # 4. Put into list and return
  return(list(m=m, u=u))
}

# Draw m and u from their full conditional distributions in SMCMC context.
# Assumes cmpdata is a list of two lists: first, the cmpdata between files 1 and 2
# Next the cmpdata between files 1 and 3, and 2 and 3. This format could be expanded
# for multiple files.
# This file also assumes link tracing is occurring
r_m_u_fc_smcmc <- function(cmpdata, Z, Z2, a, b, directratio=FALSE) {
  if (directratio) {
    stop("Direct ratio ('fast') computation not implemented for r_m_u_fc_smcmc")
  } else {
    # Tally disagreement counts.
    tmp <- disag.counts(cmpdata[[1]], c(), Z, do.trace=TRUE)
    m.pars <- tmp$match
    u.pars <- tmp$nonmatch
    tmp <- disag.counts(cmpdata[[2]], Z, Z2, do.trace=TRUE)
    m.pars <- m.pars + tmp$match
    u.pars <- u.pars + tmp$nonmatch
    # DRaw m and u from dirichlet distributions
    m <- rdirichlet.multi(alpha = m.pars + a, groups = cmpdata[[1]][[1]]$nDisagLevs)
    u <- rdirichlet.multi(alpha = u.pars + b, groups = cmpdata[[1]][[1]]$nDisagLevs)
  }
  # Return both as a list
  list(m=m, u=u)
}

################################################################################
## DIRICHLET WRAPPERS ##########################################################
################################################################################

# Function to return the density of several dirichlet distributions. The
# motivation is dealing with the m and u parameters, which are several
# dirichlet r.v.'s stacked into a single vector.
# Parameters:
#   x - the value at which to calculate the density
#   alpha - the dirichlet parameters, stacked into one vector
#   groups - how to split the one vector into separate r.v.'s.
#            sum(groups) == length(x) == length(alpha). e.g. If groups = c(3,4,5),
#            then x[1:3] ~ dirichlet(alpha[1:3]), x[4:7] ~ dirichlet([4:7]), and
#            x[8:11] ~ dirichlet([8:11]).
#   log - whether to return the log density.
# Returns:
#   prod(ddirichlet(x[group_i], alpha[group_i]))
#   (or sum, if log == TRUE)
ddirichlet.multi <- function(x, alpha, groups, log=FALSE) {
  ddirichlet <- extraDistr::ddirichlet
  lp <- 0
  prvgrplen <- 0
  for (g in seq_along(groups)) {
    lp <- lp + ddirichlet(x[(prvgrplen+1):(prvgrplen+groups[g])],
                          alpha[(prvgrplen+1):(prvgrplen+groups[g])],
                          log=TRUE)
    prvgrplen <- prvgrplen + groups[g]
  }
  if (log == FALSE) {
    return(exp(lp))
  } else {
    return(lp)
  }
}

# Function to sample from several dirichlet distributions. The
# motivation is dealing with the m and u parameters, which are several
# dirichlet r.v.'s stacked into a single vector.
# Parameters:
#   alpha - the dirichlet parameters, stacked into one vector
#   groups - how to split the one vector into separate r.v.'s.
#            sum(groups) == length(alpha). e.g. If groups = c(3,4,5),
#            then x[1:3] ~ dirichlet(alpha[1:3]), x[4:7] ~ dirichlet([4:7]), and
#            x[8:11] ~ dirichlet([8:11]).
# Returns:
#   x, a stacked dirichlet rv according to the above scheme.
rdirichlet.multi <- function(alpha, groups) {
  rdirichlet <- extraDistr::rdirichlet
  prvgrplen <- 0
  x <- rep(0, sum(groups))
  for (g in seq_along(groups)) {
    x[(prvgrplen+1):(prvgrplen+groups[g])] <- rdirichlet(n=1, alpha[(prvgrplen+1):(prvgrplen+groups[g])])
    prvgrplen <- prvgrplen + groups[g]
  }
  return(x)
}
