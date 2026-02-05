# ============================== Load balancing ================================

#' @title
#' Define chunks of vertices for analyses
#'
#' @description
#' Split a vector of vertex indices into chunks for memory-efficient
#' parallel processing. Each chunk is annotated with attributes that
#' facilitate progress reporting (chunk index and 25/50/75% milestones).
#'
#' @param iv Integer vector of vertex indices to use.
#' @param chunk_size Integer; target size of each chunk (default = 1000).
#'
#' @return A list of integer vectors, each containing vertex indices and 
#' progress marks attributes.
#' 
#' @author Serena Defina, 2024.
#' 
make_chunk_sequence <- function(iv, chunk_size = 1000) {

  chunk_seq <- split(iv, ceiling(seq_along(iv) / chunk_size))

  # Add chunk index as attribute (so it can be accessed for progress updates)
  # chunk_seq <- Map(function(chunk, i) {
  #   attr(chunk, "chunk_idx") <- i; chunk
  #   }, chunk_seq, seq_along(chunk_seq))
  progress_marks <- ceiling(chunk_size * c(0.25, 0.50, 0.75))

  for (i in seq_along(chunk_seq)) {

    attr(chunk_seq[[i]], "chunk_idx") <- i

    if (length(chunk_seq[[i]]) != chunk_size) { # e.g. last chunk
      this_chunk_size <- length(chunk_seq[[i]])
      progress_marks <- ceiling(this_chunk_size * c(0.25, 0.50, 0.75))
    }

    attr(chunk_seq[[i]], "chunk_25%") <- chunk_seq[[i]][progress_marks[1]]
    attr(chunk_seq[[i]], "chunk_50%") <- chunk_seq[[i]][progress_marks[2]]
    attr(chunk_seq[[i]], "chunk_75%") <- chunk_seq[[i]][progress_marks[3]]

  }

  return(chunk_seq)
}

#' @title
#' Run a chunked foreach loop sequentially or in parallel
#'
#' @description
#' Helper to execute a vertex-chunk loop either sequentially (e.g. for profiling)
#' or in parallel using \code{foreach} and \code{doParallel}. The expression
#' passed in \code{expr} is evaluated for each chunk in \code{chuck_seq}.
#'
#' @param n_cores Integer number of CPU cores to use. Values greater than
#'   1 enable parallel execution via a fork cluster; 1 forces sequential
#'   execution via \code{registerDoSEQ()}.
#' @param expr An expression to be evaluated in parallel or in sequence. 
#' @param progress_file Optional character string specifying a log file path. 
#'   If \code{NULL}, worker output is discarded.
#' @param seed Optional integer seed for reproducible \code{\%dopar\%} loops
#'   via \code{doRNG::registerDoRNG()}.
#' @param verbose Logical; if \code{TRUE}, prints basic messages about the
#'   selected execution mode.

#'
#' @details
#' TMP: On Unix-like systems, a fork cluster is used for parallel execution
#' (type \code{"FORK"}). On non-Unix systems this will fail anyay while 
#' FreeSurfer is required.
#'
#' @return
#' Invisibly returns \code{NULL}. Side effects are produced by \code{expr}.
#'
#' @author Serena Defina, 2026.
#' 
#' 
with_parallel <- function(n_cores,
                          expr,
                          progress_file = NULL,
                          seed = 3108,
                          verbose = TRUE) {
  # capture the expression and caller env
  expr <- substitute(expr)
  caller_env <- parent.frame()

  # save old backend info
  # was_registered <- foreach::getDoParRegistered()
  # old_name <- foreach::getDoParName()
  # old_workers <- foreach::getDoParWorkers()


  if (n_cores > 1) {
    # Check if a backend is already registered and clean it up
    if (foreach::getDoParRegistered()) {
      current_backend <- foreach::getDoParName()
      vw_message("   ! cleaning up existing parallel backend: ", current_backend,
                 verbose = verbose)
      # Unregister by switching to sequential
      foreach::registerDoSEQ()
      # Small delay to ensure cleanup completes
      Sys.sleep(0.1)
    }

    vw_message(" * preparing cluster of ", n_cores, " workers...", verbose = verbose)
    
    cl <- parallel::makeCluster(n_cores, type = "FORK", 
                                outfile = if (is.null(progress_file)) "" else progress_file)
    # NOTE: would not work on Windows anyway till freesurfer dependency is needed
    # type = ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK"))

    on.exit(parallel::stopCluster(cl), add = TRUE)

    doParallel::registerDoParallel(cl)

    # Perform %dopar% as %dorng% loops, for reproducible random numbers
    doRNG::registerDoRNG(seed)

    # Sync library paths: ensures all workers look for packages in the same place
    # as the main session
    # parallel::clusterCall(cluster, function(x) .libPaths(x), .libPaths())

  } else {
    vw_message(" * note: sequential processing...", verbose = verbose)
    # force sequential foreach
    foreach::registerDoSEQ()
  }

  # Evaluate expr exactly as if written in the caller
  eval(expr, envir = caller_env)

  return(invisible(NULL))
}


#' @title
#' Initialize per-chunk progress tracking
#'
#' @description
#' Prepare progress-tracking metadata for a given chunk of vertices. When
#' \code{verbose = TRUE}, prints a message indicating which chunk is being
#' processed and returns markers for 25\%, 50\%, and 75\% completion within
#' that chunk.
#'
#' @param chunk Integer vector of vertex indices for the current chunk.
#'   Expected to be an element of \code{chunk_seq} created by
#'   \code{\link{make_chunk_sequence}}.
#' @param chunk_seq Full list of chunks (as returned by
#'   \code{make_chunk_sequence}), used to compute the total number of chunks.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @return
#' If \code{verbose = TRUE}, a list with elements:
#' \itemize{
#'   \item \code{chunk_idx_str}: character label of the form
#'     \code{"<i>/<n_chunks>"}.
#'   \item \code{milestone_markers}: named list of integer vertex indices
#'     at 25\%, 50\%, and 75\% positions within the chunk.
#' }
#' If \code{verbose = FALSE}, returns \code{NULL} invisibly.
#'
#' @author Serena Defina, 2026.
#' 
init_progress_tracker <- function(chunk, chunk_seq, verbose) {

  if (verbose) {
    chunk_idx <- as.integer(attr(chunk, 'chunk_idx'))
    worker_id <- Sys.getpid() # useful for debugging

    chunk_idx_str <- sprintf("%d/%d", chunk_idx, length(chunk_seq))

    vw_message(sprintf("- Processing chunk %s (worker: %s)", 
                       chunk_idx_str, worker_id))
    utils::flush.console()
  
    update_freqs <- c("25%", "50%", "75%")

    milestone_markers <- sapply(update_freqs,
      function(f) as.integer(attr(chunk, paste0("chunk_", f))),
      simplify = FALSE, USE.NAMES = TRUE)
    
    return(list(chunk_idx_str, milestone_markers))
  }
  
  invisible(NULL)
}

#' @title
#' Update within-chunk progress for a (milestone) vertex
#'
#' @description
#' Check whether the current vertex index matches any pre-defined
#' within-chunk milestones (25\%, 50\%, 75\%), and print a progress
#' message if so.
#'
#' @param v Integer; current vertex index.
#' @param progress_tracker List returned by
#'   \code{\link{init_progress_tracker}}, containing \code{chunk_idx_str}
#'   and \code{milestone_markers}.
#' @param verbose Logical; if \code{TRUE}, prints progress messages.
#'
#' @return
#' Invisibly returns \code{NULL}. Used for its side effects (messages).
#'
#' @author Serena Defina, 2026.
#' 
update_progress_tracker <- function(v, progress_tracker, verbose) {

  if (verbose) {
    milestone <- names(which(progress_tracker$milestone_markers == v))
    if (length(milestone) == 1) {
      vw_message(sprintf("--- chunk %s is %s done.", progress_tracker$chunk_idx_str, milestone))
      utils::flush.console()
    }
  }
  
  invisible(NULL)
}


# ========================== External environment ==============================
# Temporarily disable parallel BLAS (and other) to prevent accidental implicit
# parallelism
# NOTE, not used. Leave that to the user insted 
with_tmp_sysenv <- function(tmp_sysenv = c(OMP_NUM_THREADS = 1,
                                           MKL_NUM_THREADS = 1,
                                           OPENBLAS_NUM_THREADS = 1,
                                           VECLIB_MAXIMUM_THREADS = 1,
                                           NUMEXPR_NUM_THREADS = 1),
                            code) {
  sysenv <- Sys.getenv(names(tmp_sysenv), unset = NA, names = TRUE)

  for (nm in names(tmp_sysenv)) Sys.setenv(nm = tmp_sysenv[[nm]])

  on.exit({
    for (nm in names(sysenv)) {
      if (is.na(sysenv[[nm]])) Sys.unsetenv(nm)
      else Sys.setenv(nm = sysenv[[nm]])
    }
  }, add = TRUE)

  force(code)
}