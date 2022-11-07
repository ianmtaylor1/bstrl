# This code generates the geco_small_result data object, a pre-run result of
# linking the geco_small dataset. Useful for examples and testing code to post-process
# results


# NOTE: since this requires the package to calculate the result, so this file should
# be run *after* the development version of bstrl has been installed locally,
# potentially with an updated geco_small dataset or updated RL functions. To avoid
# any conflicts with documentation or vignettes which use this object, install
# the development version with devtools::install(quick=TRUE)
library(bstrl)

data(geco_small)

# Names of the columns on which to perform linkage
fieldnames <- c("given.name", "surname", "age", "occup", "extra1", "extra2", "extra3", "extra4", "extra5", "extra6")

# How to compare each of the fields
types <- c("lv", "lv", # First name and last name use normalized edit distance
           "bi", "bi", "bi", "bi", "bi", "bi", "bi", "bi") # All others binary equal/unequal
breaks <- c(0, 0.25, 0.5) # Break continuous difference measures into 4 levels using these split points

# MCMC details
nIter <- 600
burn <- 100

# Steps mimicking the process in the vignette.
res.twofile <- bipartiteRL(geco_small[[1]], geco_small[[2]],
                           flds = fieldnames, types = types, breaks = breaks,
                           nIter = nIter, burn = burn,
                           seed = 0)
res.pprb3 <- PPRBupdate(res.twofile, geco_small[[3]], # Comparison details are stored with previous result
                        nIter = 600, burn = 100,
                        seed = 0,
                        refresh = 0.05)
res.pprb4 <- PPRBupdate(res.pprb3, geco_small[[4]], # Comparison details are stored with previous result
                        nIter = 600, burn = 100,
                        seed = 0,
                        refresh = 0.05)

# Final object will be 4-file PPRB result
geco_small_result <- res.pprb4

usethis::use_data(geco_small_result, overwrite = TRUE)
