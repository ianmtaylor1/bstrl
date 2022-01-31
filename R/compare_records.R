# Wrap the BRL compareRecords function and do preprocessing of the comparison
# data for efficiency later
compareRecords <- function(df1, df2, flds = NULL, flds1 = NULL, flds2 = NULL,
                           types = NULL, breaks = c(0, 0.25, 0.5)) {
  cmp <- BRL::compareRecords(df1, df2, flds = flds, flds1 = flds1, flds2 = flds2,
                             types = types, breaks = breaks)
  attr(cmp$comparisons, "totals") <- colSums(cmp$comparisons)
  cmp
}
