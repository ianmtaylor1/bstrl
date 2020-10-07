# This file contains functions related to the full conditional distribution
# of m and u.

# It depends on the package extraDistr for calculating the density and sampling
# from Dirichlet distributions

# Draw m and u from their full conditional distributions.
# Parameters:
#   cmpdata.1to2, cmpdata.1to3 - comparisons with the "new" third file.
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
r_m_u_fc <- function(cmpdata.1to3, cmpdata.2to3, Z, Z2,
                     a, b, m.prev.pars, u.prev.pars, trace=FALSE) {
  # 1. Create the additions to m and u hyperparameters derived from links with
  # the newest (third) file.
  matchrows.2to3 <- matchrows(cmpdata.2to3, Z2, offset=cmpdata.1to3$n1)
  if (trace) {
    matchrows.1to3 <- matchrows(cmpdata.1to3,
                                trace(Z2, Z, steps=1, offset=cmpdata.1to3$n1),
                                offset=0)
    noncandrows.1to3 <- c()
  } else {
    matchrows.1to3 <- matchrows(cmpdata.1to3, Z2, offset=0)
    noncandrows.1to3 <- noncandrows(cmpdata.1to3, Z, offset=0)
  }
  m.new.pars <- colSums(cmpdata.1to3$comparisons[matchrows.1to3,,drop=FALSE]) +
    colSums(cmpdata.2to3$comparisons[matchrows.2to3,,drop=FALSE])
  u.new.pars <- attr(cmpdata.1to3$comparisons, "totals") +
    attr(cmpdata.2to3$comparisons, "totals") -
    m.new.pars -
    colSums(cmpdata.1to3$comparisons[noncandrows.1to3,,drop=FALSE])
  # 2. combine a, m.prev.pars and new m pars to sample m
  m <- rdirichlet.multi(alpha=m.new.pars + m.prev.pars + a, groups=cmpdata.1to3$nDisagLevs)
  # 3. combine b, u.prev.pars and new u pars to sample u
  u <- rdirichlet.multi(alpha=u.new.pars + u.prev.pars + b, groups=cmpdata.1to3$nDisagLevs)
  # 4. Put into list and return
  return(list(m=m, u=u))
}

################################################################################
## HELPER FUNCTIONS ############################################################
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
  for (g in 1:length(groups)) {
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
  for (g in 1:length(groups)) {
    x[(prvgrplen+1):(prvgrplen+groups[g])] <- rdirichlet(n=1, alpha[(prvgrplen+1):(prvgrplen+groups[g])])
    prvgrplen <- prvgrplen + groups[g]
  }
  return(x)
}
