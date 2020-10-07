# This file contains functions related to the full conditional distribution
# of m and u.

# It depends on the package extraDistr for calculating the density and sampling
# from Dirichlet distributions


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
