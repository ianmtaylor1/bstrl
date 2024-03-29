# This file contains functions related to the full conditional distribution
# of m and u.

# It depends on the package extraDistr for calculating the density and sampling
# from Dirichlet distributions

# Draw m and u from their full conditional distributions.
# Parameters:
#   lastfilecmps - a list of comparison data objects. There are k-1 total objects
#             in the list. The first compares file 1 to file k, and so on, until
#             the last which compares file k-1 to file k. All objects should
#             therefore have equal n2 values.
#   sl - streaminglinks object defining current state of links between files.
#   a, b - prior hyperparameters for the distributions of m and u, respectively.
#          See ?BRL::BRL for explanation
#   m.prev.pars - additions to hyperparameters for m from matches between
#          previous files. These are computed in a post-processing step of the
#          samples of Z, and can just be passed directly here. Provides the
#          contribution of the comparisons between all k-1 previous files, i.e.
#          comparisons not included in lastfilecmps.
#   u.prev.pars - same as m.prev.pars, but for u.
# Return:
#   A list with two elements, m and u, which contain the fc sampled m and u.
r_m_u_fc_pprb <- function(lastfilecmps, sl, a, b, m.prev.pars, u.prev.pars) {
  # 1. Determine the contribution of the latest file's comparison data to the
  # full conditional distribution
  counts <- disag.counts.lastfile(lastfilecmps, sl)
  m.new.pars <- counts$match
  u.new.pars <- counts$nonmatch
  # 2. combine a, m.prev.pars and new m pars to sample m
  m <- rdirichlet.multi(alpha=m.new.pars + m.prev.pars + a, groups=lastfilecmps[[1]]$nDisagLevs)
  # 3. combine b, u.prev.pars and new u pars to sample u
  u <- rdirichlet.multi(alpha=u.new.pars + u.prev.pars + b, groups=lastfilecmps[[1]]$nDisagLevs)
  # 4. Put into list and return
  return(list(m=m, u=u))
}

# Draw m and u from their full conditional distributions.
# Parameters:
#   cmpdata - a list of lists of comparison data objects. (Standard triangular
#             comparison data format.)
#   sl - streaminglinks object defining current state of links between files.
#   a, b - prior hyperparameters for the distributions of m and u, respectively.
#          See ?BRL::BRL for explanation
# Return:
#   A list with two elements, m and u, which contain the fc sampled m and u.
r_m_u_fc_smcmc <- function(cmpdata, sl, a, b) {
  # 1. Determine the contribution of the latest file's comparison data to the
  # full conditional distribution
  counts <- disag.counts.allfiles(cmpdata, sl)
  m.new.pars <- counts$match
  u.new.pars <- counts$nonmatch
  # 2. combine a, m.prev.pars and new m pars to sample m
  m <- rdirichlet.multi(alpha=m.new.pars + a, groups=cmpdata[[1]][[1]]$nDisagLevs)
  # 3. combine b, u.prev.pars and new u pars to sample u
  u <- rdirichlet.multi(alpha=u.new.pars + b, groups=cmpdata[[1]][[1]]$nDisagLevs)
  # 4. Put into list and return
  return(list(m=m, u=u))
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
