# This file contains functions to create, manipulate, and query an object
# representing the links between streaming files (i.e. all the Z^{(k)} vectors)

# Create an object of S3 class "streaminglinks", for files of the specified size
# filesizes should be a vector of length k, for k files, containing the sizes
# of each file in order
streaminglinks <- function(filesizes) {
  stopifnot(length(filesizes) >= 2)
  stopifnot(filesizes > 0)
  stopifnot(filesizes == floor(filesizes))
  # Create an initially unlinked configuration
  structure(
    list(
      Z = seq_len(sum(filesizes)),
      W = seq_len(sum(filesizes)),
      ns = filesizes
    ),
    class = "streaminglinks"
  )
}



#### Querying the link object


# Return the number of files linked by the given streaminglinks object
nfiles <- function(sl) {
  length(sl$ns)
}


# Return True/False whether the two record are coreferent
islinked <- function(sl, file1, record1, file2, record2) {
  stopifnot(file1 < file2)
  stopifnot(record1 <= sl$ns[file1], record1 >= 1)
  stopifnot(record2 <= sl$ns[file2], record2 >= 1)

  gidx1 <- local.to.global(sl$ns, file1, record1)
  gidx2 <- local.to.global(sl$ns, file2, record2)

  # Walk back from the second record until reaching the end of the chain,
  # reaching the first record, or passing the first record
  cur <- gidx2
  linked <- FALSE
  while ((sl$Z[cur] < cur) && (sl$Z[cur] >= gidx1)) {
    cur <- sl$Z[cur]
    if (cur == gidx1) linked <- TRUE
  }

  return(linked)
}

# Return a list of all linked pairs (directly or transitively)
alllinks <- function(sl, idx=c("global", "local")) {
  idx <- match.arg(idx)

  # Initial, direct links. Will form the starting point of returned links
  rhs <- which(sl$Z < seq_len(sum(sl$ns)))
  lhs <- sl$Z[rhs]
  allright <- rhs
  allleft <- lhs

  # We know that only records with a link in Z can ever be linked to anything
  # downstream, so we start with those and filter at each step
  while (length(rhs) > 0) {
    deadends <- (sl$Z[lhs] == lhs) # boolean vector indicating no further links
    rhs <- rhs[!deadends]
    lhs <- sl$Z[lhs[!deadends]]
    allright <- c(allright, rhs)
    allleft <- c(allleft, lhs)
  }
  # Now we have all pairs as global indexes

  # Optionally translate to local indices, and return
  if (idx == "local") {
    idx1 <- global.to.local(sl$ns, allleft)
    idx2 <- global.to.local(sl$ns, allright)
    return(list(file1=idx1$file, record1=idx1$idx,
                file2=idx2$file, record2=idx2$idx))
  } else {
    return(list(idx1 = allleft, idx2 = allright))
  }
}

#### Manipulating link objects


# In the given streaminglinks object, unlink the record at (file, record) from
# its next upstream (later file) link, returning the modified streaminglinks
# object. If the supplied record is not linked to anything upstream, do nothing.
unlink.up <- function(sl, file, record) {
  globalidx <- local.to.global(sl$ns, file, record)
  if (sl$W[globalidx] > globalidx) {
    sl$Z[sl$W[globalidx]] <- sl$W[globalidx]
    sl$W[globalidx] <- globalidx
  }
  return(sl)
}

# In the given streaminglinks object, unlink the record at (file, record) from
# its next downstream (earlier file) link, returning the modified streaminglinks
# object. If the supplied record is not linked to anything downstream, do nothing.
unlink.down <- function(sl, file, record) {
  globalidx <- local.to.global(sl$ns, file, record)
  if (sl$Z[globalidx] < globalidx) {
    sl$W[sl$Z[globalidx]] <- sl$Z[globalidx]
    sl$Z[globalidx] <- globalidx
  }
  return(sl)
}

# Add a link between the specified records
add.link <- function(sl, file1, record1, file2, record2,
                     conflict=c("null", "warn", "error", "overwrite")) {
  stopifnot(file1 < file2)
  stopifnot(record1 <= sl$ns[file1], record1 >= 1)
  stopifnot(record2 <= sl$ns[file2], record2 >= 1)
  conflict <- match.arg(conflict)
  # Get record global indices
  gidx1 <- local.to.global(sl$ns, file1, record1)
  gidx2 <- local.to.global(sl$ns, file2, record2)
  # Check if there are conflicts
  if ((sl$Z[gidx2] < gidx2) || (sl$W[gidx1] > gidx1)) {
    # If yes, we have to decide how to handle it
    if (conflict == "null") {
      return(NULL)
    } else if (conflict == "error") {
      stop("Cannot link file", file1, "record", record1, "to file", file2, "record", record2)
    } else {
      if(conflict == "warn") {
        warning("Conflicts in linking file", file1, "record", record1, "to file", file2, "record", record2,
                "- overwriting existing links")
      }
      # For either 'warn' or 'overwrite', unlink existing links and continue
      sl <- unlink.down(sl, file2, record2)
      sl <- unlink.up(sl, file1, record1)
    }
  }
  # If we get here, we can create the link and return the modified streaming link object
  sl$Z[gidx2] <- gidx1
  sl$W[gidx1] <- gidx2
  return(sl)
}


#### Index translation


local.to.global <- function(filesizes, file, record) {
  mapply(
    function(f, r, s) {sum(s[seq_len(f-1)]) + r},
    file, record,
    MoreArgs=list(filesizes)
  )
}

global.to.local <- function(filesizes, globalrec) {
  # Thresholds of the cumulative number of records through file k
  thresh <- cumsum(filesizes)

  # Find the file number of each record
  fileno <- findInterval(globalrec, thresh, left.open=T) + 1

  # Find the index within each file
  idxno <- globalrec - (thresh - filesizes)[fileno]

  return(list(file=fileno, idx=idxno))
}

