% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PPRBupdate.R
\name{PPRBupdate}
\alias{PPRBupdate}
\title{Perform a PPRB update of record linkage with a new file}
\usage{
PPRBupdate(
  state,
  newfile,
  flds = NULL,
  nIter = NULL,
  burn = 0,
  blocksize = NULL,
  threestep = TRUE,
  seed = 0,
  refresh = 0.1
)
}
\arguments{
\item{state}{Existing record linkage state. Returned by either bipartiteRL,
PPRBupdate, SMCMCupdate, or multifileRL.}

\item{newfile}{A data.frame representing the new file: one row per record}

\item{flds}{Names of fields in the new file to use for comparison. Only used
if no global field names were specified in bipartiteRL initially.}

\item{nIter}{Number of iterations for which to run the PPRB sampler. By
default, this is the same as the number of samples present in 'state'.}

\item{burn}{Number of initial iterations to discard. The total number of
samples returned is nIter - burn.}

\item{blocksize}{Size of blocks to use for locally balanced proposals.
Default performs unblocked locally balanced proposals.}

\item{threestep}{Whether to perform three Gibbs sampling steps per iteration,
with past Z's updated with PPRB, m and u updated with full conditionals,
and the current Z updated with locally balanced proposals. If false, a two
step Gibbs sampler is used where past Z's, m and u are updated together
using PPRB and the current Z is updated with locally balanced proposals}

\item{seed}{Random seed to set at the beginning of the MCMC run}

\item{refresh}{How often to output an update including the iteration number
and percent complete. If refresh >= 1, taken as a number of iterations
between messages (rounded). If 0 < refresh < 1, taken as the proportion of
nIter. If refresh == 0, no messages are displayed.}
}
\value{
An object of class 'bstrlstate' containing posterior samples and
necessary metadata for passing to future streaming updates.
}
\description{
Perform a PPRB update of record linkage with a new file
}
\examples{
data(geco_small)
data(geco_small_result)

# Add fifth file to previous four-file link result
file5.result <- PPRBupdate(geco_small_result, geco_small[[5]],
                           nIter=2, burn=1) # Very small run for example

}
