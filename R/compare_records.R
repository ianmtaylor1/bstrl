# Wrap the BRL compareRecords function and do preprocessing of the comparison
# data for efficiency later
compareRecords <- function(df1, df2, flds = NULL, flds1 = NULL, flds2 = NULL,
                           types = NULL, breaks = c(0, 0.25, 0.5)) {
  if (nrow(df1) >= nrow(df2)) {
    cmp <- BRL::compareRecords(df1, df2, flds = flds, flds1 = flds1, flds2 = flds2,
                               types = types, breaks = breaks)
  } else {
    # BRL assumes file 1 is at least as large as file 2. If that's not true,
    # Reverse the call then flip everything back around
    cmp <- BRL::compareRecords(df2, df1, flds = flds, flds1 = flds2, flds2 = flds1,
                               types = types, breaks = breaks)
    # Swap around the comparison matrix rows
    M <- matrix(seq_len(cmp$n1 * cmp$n2), nrow=cmp$n2, ncol=cmp$n1, byrow=T)
    reorder <- c(M)
    cmp$comparisons <- cmp$comparisons[reorder,]
    # Swap file 1 and 2 in compFields
    file1names <- cmp$compFields$file1
    cmp$compFields$file1 <- cmp$compFields$file2
    cmp$compFields$file2 <- file1names
    # Swap filesizes
    n1 <- cmp$n1
    cmp$n1 <- cmp$n2
    cmp$n2 <- n1
  }
  # Attach column sums
  attr(cmp$comparisons, "totals") <- colSums(cmp$comparisons)
  cmp
}
