# Functions in this perform bipartite linkage, by wrapping the BRL package
# and appending necessary extra information


#' Perform baseline bipartite record linkage before streaming updates
#'
#' This function establishes a baseline linkage between two files which can be
#' built upon with streaming updates adding more files. It outsources the linkage
#' work to the BRL package and appends information to the object which will allow
#' streaming record linkage to continue
#'
#' @param df1,df2 Files 1 and 2 as dataframes where each row is a record and
#'   each column is a field.
#' @param flds Names of the fields on which to compare the records in each file
#' @param flds1,flds2 Allows specifying field names differently for each file.
#' @param types Types of comparisons to use for each field
#' @param breaks Breaks to use for Levenshtein distance on string fields
#' @param nIter,burn MCMC run length parameters. The returned number of samples
#'   is nIter - burn.
#' @param a,b Prior parameters for m and u, respectively.
#' @param aBM,bBM Prior parameters for beta-linkage prior.
#' @param seed Random seed to set at beginning of MCMC run
#'
#' @return A list with class "bstrlstate" which can be passed to future streaming
#'   updates.
#'
#' @examples
#' data(geco_small)
#'
#' # Names of the columns on which to perform linkage
#' fieldnames <- c("given.name", "surname", "age", "occup",
#'                 "extra1", "extra2", "extra3", "extra4", "extra5", "extra6")
#'
#' # How to compare each of the fields
#' # First name and last name use normalized edit distance
#' # All others binary equal/unequal
#' types <- c("lv", "lv",
#'            "bi", "bi", "bi", "bi", "bi", "bi", "bi", "bi")
#' # Break continuous difference measures into 4 levels using these split points
#' breaks <- c(0, 0.25, 0.5)
#'
#' res.twofile <- bipartiteRL(geco_small[[1]], geco_small[[2]],
#'                            flds = fieldnames, types = types, breaks = breaks,
#'                            nIter = 10, burn = 5, # Very small number of samples
#'                            seed = 0)
#'
#' @export
bipartiteRL <- function(df1, df2,
                        flds=NULL, flds1=NULL, flds2=NULL, types=NULL, breaks=c(0,.25,.5),
                        nIter=1000, burn=round(nIter*.1), a=1, b=1, aBM=1, bBM=1,
                        seed=0) {

  # 'burn' is not a parameter of the other functions, so we control it here
  if( !is.numeric(burn) | (burn<0) | (burn>=nIter) )
    stop("burn should be an integer that satisfies 0 <= burn < nIter")

  # create comparison data
  cmpdata <- compareRecords(df1, df2, flds, flds1, flds2, types, breaks)

  # Perform linkage
  chain <- BRL::bipartiteGibbs(cmpdata, nIter, a, b, aBM, bBM, seed)

  # Store files and comparison data in format we will use later
  files <- list(df1, df2)
  comparisons <- list(list(cmpdata))

  # Store comparison and prior parameters
  priors <- list(a=a, b=b, aBM=aBM, bBM=bBM)
  cmpdetails <- list(fldsglobal=flds, types=types, breaks=breaks,
                     flds=list(flds1, flds2))

  # Store summary statistics of comparisons for use later
  m.fc.pars <- matrix(0, nrow=nrow(chain$m), ncol=nIter)
  u.fc.pars <- matrix(0, nrow=nrow(chain$m), ncol=nIter)
  for (i in 1:nIter) {
    tmp <- disag.counts.lastfile(list(cmpdata), streaminglinks(c(nrow(df1), nrow(df2)), chain$Z[,i]))
    m.fc.pars[,i] <- tmp$match
    u.fc.pars[,i] <- tmp$nonmatch
  }

  # Burn and return complete structure
  iterfilter <- setdiff(seq_len(nIter), seq_len(burn))
  structure(
    list(
      Z = chain$Z[,iterfilter,drop=F],
      m = chain$m[,iterfilter,drop=F],
      u = chain$u[,iterfilter,drop=F],
      files = files,
      comparisons = comparisons,
      cmpdetails = cmpdetails,
      priors = priors,
      m.fc.pars = m.fc.pars[,iterfilter,drop=F],
      u.fc.pars = u.fc.pars[,iterfilter,drop=F]
    ),
    class = "bstrlstate"
  )
}

