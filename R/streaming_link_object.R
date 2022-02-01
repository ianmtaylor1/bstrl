# This file contains functions to create, manipulate, and query an object
# representing the links between streaming files (i.e. all the Z^{(k)} vectors)

# Create an object of S3 class "streaminglinks", for files of the specified size
# filesizes should be a vector of length k, for k files, containing the sizes
# of each file in order
streaminglinks <- function(filesizes, Z=NULL) {
  stopifnot(length(filesizes) >= 2)
  stopifnot(filesizes > 0)
  stopifnot(filesizes == floor(filesizes))
  # Process supplied z, if any
  if (is.null(Z)) {
    Z <- seq_len(sum(filesizes))
    W <- seq_len(sum(filesizes))
  } else if (length(Z) == sum(filesizes)) {
    stopifnot(checkvalidZ(Z)) # Sanity check for valid Z
    W <- WfromZ(Z)
  } else if (length(Z) == sum(filesizes[-1])) {
    Z <- c(seq_len(filesizes[1]), Z)
    stopifnot(checkvalidZ(Z)) # Sanity check for valid Z
    W <- WfromZ(Z)
  } else {
    stop("Supplied Z must have length equal to total number of records, or total number of records except file 1.")
  }
  # Create an initially unlinked configuration
  structure(
    list(
      Z = Z,
      W = W,
      ns = filesizes
    ),
    class = "streaminglinks"
  )
}

# Take a concatenated Z vector and return the corresponding W (reversing the
# links to point up from earlier record to later record)
WfromZ <- function(Z) {
  W <- seq_along(Z)
  j <- which(Z <= seq_along(Z))
  i <- Z[j]
  W[i] <- j
  W
}

# Check if a concatenated Z vector is valid
checkvalidZ <- function(Z) {
  if (any(Z > seq_along(Z))) return(FALSE)
  links <- Z[Z < seq_along(Z)]
  if (length(unique(links)) != length(links)) return(FALSE)
  return(TRUE)
}

#### Querying the link object

# Return whether the streaming link object is valid: Z represents a valid
# link state and all its internal components are in agreement
# This should never be needed, but is here just in case.
isvalid <- function(sl) {
  (length(sl$Z) == sum(sl$ns)) && checkvalidZ(sl$Z) && all(sl$W == WfromZ(sl$Z))
}

# Return the number of files linked by the given streaminglinks object
nfiles <- function(sl) {
  length(sl$ns)
}

# Return a saveable Z vector representing this object, discarding redundant
# information (i.e. the self-links from file 1 to nowhere)
savestate <- function(sl) {
  sl$Z[seq(sl$ns[1] + 1, length(sl$Z))]
}

# Return True/False whether the two record are coreferent
islinked <- function(sl, file1, record1, file2, record2) {
  stopifnot(file1 < file2)
  stopifnot(record1 <= sl$ns[file1], record1 >= 1)
  stopifnot(record2 <= sl$ns[file2], record2 >= 1)

  gidx1 <- local.to.global(sl$ns, file1, record1)
  gidx2 <- local.to.global(sl$ns, file2, record2)

  islinked.gl(sl, gidx1, gidx2)
}

islinked.gl <- function(sl, gidx1, gidx2) {
  stopifnot(gidx1 < gidx2)
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

# Create blocks for blocked locally balanced proposals
# file = the file number whose links we intend to update. The blocks will
#   contain a subset of the records from this file, and a subset of the records
#   from all previous files
# blocksize = how
createblocks <- function(sl, file, blocksize=NULL,
                         method=c("approximate", "exact")) {
  method <- match.arg(method)
  stopifnot(file >= 2, file <= length(sl$ns))

  # File sizes
  nk <- sl$ns[file]
  nprev <- sum(sl$ns[seq_len(file-1)])
  # Candidates for block inclusion
  totalj <- nprev + seq_len(nk)
  totali <- seq_len(nprev)
  totali <- totali[(sl$W[totali] == totali) | (sl$W[totali] > nprev)]

  # First, correct blocksize
  if (is.null(blocksize)) {
    bsright <- nk
    bsleft <- length(totali)
  } else {
    bsright <- bsleft <- blocksize
  }
  if (bsright > nk) {
    bsright <- nk
  }
  if (bsleft > length(totali)) {
    bsleft <- length(totali)
  }

  iblock <- jblock <- c()
  if (method == "approximate") {
    # Now cycle until we have nonempty blocks
    while ((length(iblock) == 0) || (length(jblock) == 0)) {
      jblock <- totalj[sample(nk, bsright)]
      iblock <- totali[sample(length(totali), bsleft)]
      # Remove anything not linking within the block
      ikeep <- ((sl$W[iblock] == iblock) | (sl$W[iblock] %in% jblock))
      jkeep <- ((sl$Z[jblock] == jblock) | (sl$Z[jblock] %in% iblock))
      iblock <- iblock[ikeep]
      jblock <- jblock[jkeep]
    }
  } else if (method == "exact") {
    # TODO: algorithm to select blocks exactly matching blocksize
    stop("Exact block selection not implemented")
  }

  return(list(iblock=iblock, jblock=jblock))
}


#### Manipulating link objects


# In the given streaminglinks object, unlink the record at (file, record) from
# its next upstream (later file) link, returning the modified streaminglinks
# object. If the supplied record is not linked to anything upstream, do nothing.
unlink.up <- function(sl, file, record) {
  globalidx <- local.to.global(sl$ns, file, record)
  unlink.up.gl(sl, globalidx)
}

unlink.up.gl <- function(sl, globalidx) {
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
  unlink.down.gl(sl, globalidx)
}

unlink.down.gl <- function(sl, globalidx) {
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

  add.link.gl(sl, gidx1, gidx2, conflict)
}

add.link.gl <- function(sl, gidx1, gidx2,
                        conflict=c("null", "warn", "error", "overwrite")) {
  stopifnot(gidx1 < gidx2)
  conflict <- match.arg(conflict)
  # Check if there are conflicts
  if ((sl$Z[gidx2] < gidx2) || (sl$W[gidx1] > gidx1)) {
    # If yes, we have to decide how to handle it
    if (conflict == "null") {
      return(NULL)
    } else if (conflict == "error") {
      stop("Cannot link global record ", gidx1, " to global record ", gidx2)
    } else {
      if(conflict == "warn") {
        warning("Conflicts in linking global record ", gidx1, " to global record ", gidx2,
                "- overwriting existing links")
      }
      # For either 'warn' or 'overwrite', unlink existing links and continue
      sl <- unlink.down.gl(sl, gidx2)
      sl <- unlink.up.gl(sl, gidx1)
    }
  }
  # If we get here, we can create the link and return the modified streaming link object
  sl$Z[gidx2] <- gidx1
  sl$W[gidx1] <- gidx2
  return(sl)
}

# Swap out links among the first k-1 files of sl with the value supplied, Zpre.
# Zpre can be a vector of length n1 + ... + n(k-1) OR length n2 + ... + n(k-1).
# If the latter, the initial sequence is implied as seq_len(n1).
# This is used to perform the PPRB step which replaces previous Z's with
# resampled values from the prior stage.
swapprefix <- function(sl, Zpre,
                       conflict=c("null", "error")) {
  conflict <- match.arg(conflict)

  nfiles <- length(sl$ns)

  # Add redundant sequence to beginning if necessary
  if (length(Zpre) == sum(sl$ns[seq_len(nfiles - 2) + 1])) {
    Zpre <- c(seq_len(sl$ns[1]), Zpre)
  } else if (length(Zpre) != sum(sl$ns[seq_len(nfiles - 1)])) {
    stop("Zpre has incorrect length")
  }

  # Create new Z from prefix and last file
  Znew <- c(Zpre, sl$Z[seq(length(Zpre) + 1, sum(sl$ns))])

  # Construct into new sl object, catch any errors
  slnew <- tryCatch(
    streaminglinks(sl$ns, Znew),
    error = function(c) NULL
  )

  if (is.null(slnew)) {
    # If there is a conflict, we have to decide how to handle it. Too many
    # things potentially wrong to overwrite or warn, must fail hard. Either
    # return null or raise error.
    if (conflict == "null") {
      return(NULL)
    } else if (conflict == "error") {
      stop("Cannot link file", file1, "record", record1, "to file", file2, "record", record2)
    }
  }

  return(slnew)
}


# Perform either an add, delete, swap, or double-swap operation on index i and
# index j (globally indexed), for locally balanced proposals
performstep <- function(sl, i, j) {
  stopifnot(i >= 1, i < j, j <= length(sl$Z))
  # Save the current link from each for easy reference
  i.linkedto <- sl$W[i]
  j.linkedto <- sl$Z[j]
  # Check possibilities and perform steps
  slnew <- sl # New state
  reverse <- c(i,j) # Pair which can undo the move from the new state
  if (i.linkedto == j) { # Delete
    slnew <- unlink.down.gl(sl, j)
  } else if ((i.linkedto == i) && (j.linkedto == j)) { # Add
    slnew <- add.link.gl(sl, i, j, conflict="error") # Should be no conflict
  } else if ((i.linkedto == i) && (j.linkedto < j)) { # Single-swap 1
    slnew <- add.link.gl(sl, i, j, conflict="overwrite")
    reverse <- c(j.linkedto, j) # Reverse single-swap
  } else if ((i.linkedto > i) && (j.linkedto == j)) { # Single-swap 2
    slnew <- add.link.gl(sl, i, j, conflict="overwrite")
    reverse <- c(i, i.linkedto) # Reverse single-swap
  } else if ((i.linkedto > i) && (j.linkedto < j)) { # Double-swap
    slnew <- add.link.gl(add.link.gl(sl, i, j, conflict="overwrite"),
                         i.linkedto, j.linkedto, conflict="error")
    reverse <- c(i, i.linkedto)
  } else { # Catch-all
    stop("Unexpected error")
  }

  list(state=slnew, reverse=reverse)
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

