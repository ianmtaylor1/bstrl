% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocessing.R
\name{precision}
\alias{precision}
\title{Calculate the precision of estimated links relative to true links}
\usage{
precision(sl.est, sl.true)
}
\arguments{
\item{sl.est}{streaminglinks object representing link estimates}

\item{sl.true}{streaminglinks object representing true links}
}
\value{
The precision of the estimated links.
}
\description{
Calculate the precision of estimated links relative to true links
}
\examples{
data(geco_small)
data(geco_small_result)

sl.true <- fromentities(geco_small[[1]]$entity, geco_small[[2]]$entity,
                        geco_small[[3]]$entity, geco_small[[4]]$entity)

posterior <- extractlinks(geco_small_result)
# Compare one posterior sample to previously computed known truth
class(sl.true)
precision(posterior[[42]], sl.true)

}
