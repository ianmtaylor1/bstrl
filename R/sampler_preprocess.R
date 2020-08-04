# This file contains functions for preprocessing data to create more useful
# structures for the sampler. E.g. precomputing commonly needed but static
# summaries of input data.

# This function takes an object comparing two files, as returned by
# BRL::compareRecords, and preprocesses it into a form that is most useful to
# the bstrl sampler.
preproc.cmpdata <- function(compared) {
  # Extract the binary comparison matrix
  m <- compared$comparisons
  # Attach the column sums as an attribute. Column sums are needed in likelihood
  # functions, and computing them over and over again adds up in terms of comp
  # time.
  attr(m, "totals") <- colSums(m)
  # Return the matrix with attached column sums
  return(m)
}
