% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMCMCupdate.R
\name{SMCMCupdate}
\alias{SMCMCupdate}
\title{Perform an SMCMC update of record linkage with a new file}
\usage{
SMCMCupdate(
  state,
  newfile,
  flds = NULL,
  nIter.jumping = 5,
  nIter.transition = 10,
  cores = 1,
  proposals.jumping = c("component", "LB"),
  proposals.transition = c("LB", "component"),
  blocksize = NULL,
  seed = 0
)
}
\arguments{
\item{state}{Existing record linkage state. Returned by either bipartiteRL,
PPRBupdate, SMCMCupdate, or multifileRL. This is used as the ensemble of
samples in SMCMC update}

\item{newfile}{A data.frame representing the new file: one row per record}

\item{flds}{Names of fields in the new file to use for comparison. Only used
if no global field names were specified in bipartiteRL initially.}

\item{nIter.jumping, nIter.transition}{Number of iterations to use in the
jumping kernel and transition kernel, respectively, for each ensemble
sample.}

\item{cores}{Number of cores to use for parallel execution. If cores == 1,
update is run sequentially. A cluster is created using
parallel::makeCluster().}

\item{proposals.jumping, proposals.transition}{Which kernel to use for Z
updates in the jumping and transition kernels, respectively.}

\item{blocksize}{Size of blocks to use for locally balanced proposals.
Default performs unblocked locally balanced proposals.}

\item{seed}{Random seed to set at the beginning of the MCMC run. This is
ignored if cores > 1.}
}
\value{
An object of class 'bstrlstate' containing posterior samples and
necessary metadata for passing to future streaming updates.
}
\description{
Perform an SMCMC update of record linkage with a new file
}
\examples{
data(geco_small)
data(geco_small_result)

# Add fifth file to previous four-file link result
filtered <- thinsamples(geco_small_result, 2) # Filter ensemble to 2 - very small for example
file5.result <- SMCMCupdate(filtered, geco_small[[5]],
                            nIter.jumping=1, nIter.transition=1, # Very small run for example
                            proposals.jumping="LB", proposals.transition="LB",
                            blocksize=5)

}
