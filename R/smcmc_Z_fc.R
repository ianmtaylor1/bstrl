# This file contains componentwise full conditional samplers for Z and Z2
# 3/23/2021 - Currently used with the SMCMC transition and approximate
# jumping kernels. To do: make flexible enough to deal with more than 3 files


# Z, Z2, m, u - current parameter values
# cmpdata - a list-of-lists format of comparison data
# aBM, bBM - Z/Z2 prior parameters.
#
# Returns: the new value of Z
r_Z_fc_smcmc <- function(Z, Z2, m, u, cmpdata, aBM, bBM) {
  # For all indices in Z
    # For all possible values
      # Evaluate likelihood(s) and prior(s)
    # Normalize?
    # Select new value
  # Return
  Z
}

# Z, Z2, m, u - current parameter values
# cmpdata - a list-of-lists format of comparison data
# aBM, bBM - Z/Z2 prior parameters.
#
# Returns: the new value of Z2
r_Z2_fc_smcmc <- function(Z, Z2, m, u, cmpdata, aBM, bBM) {
  # For all indices in Z
    # For all possible values
      # Evaluate likelihood(s) and prior(s)
    # Normalize?
    # Select new value
  # Return
  Z2
}
