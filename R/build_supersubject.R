#' @title
#' Build "supersubject" by stacking all vertex data in one large file-backed
#' matrix with dimensions n_subjects x n_vertices.
#'
#' @param subj_dir : path to the FreeSurfer data, this expects a verywise structure.
#' @param folder_ids : the vector of observations to include. This holds the relative
#' path (from `subj_dir`) to the FreeSurfer data folder (e.g. "site1/sub-1_ses-01").
#' @param supsubj_dir : output path, where logs, backing files and the matrix itself
#' (if `save_rds == TRUE`) will be stored.
#' @param measure : vertex-wise measure, used to identify files.
#' @param hemi : hemisphere, used to identify files.
#' @param n_cores : number of cores to use for parallel processing.
#' @param fwhmc : (default = "fwhm10") full-width half maximum value, used to identify files.
#' @param fs_template : (default = "fsaverage") template on which to register vertex-wise data.
#' The following values are accepted:
#'  * fsaverage (default) = 163842 vertices (highest resolution),
#'  * fsaverage6 = 40962 vertices,
#'  * fsaverage5 = 10242 vertices,
#'  * fsaverage4 = 2562 vertices,
#'  * fsaverage3 = 642 vertices
#' Note that, at the moment, these are only used to downsample the brain map, for faster
#' model tuning. `verywise` expects the input data to be always registered on the "fsaverage"
#' template and the final analyses should also be run using `fs_template = "fsaverage"`
#' to avoid (small) imprecisions in vertex registration and smoothing.
#' @param backing : (default = `supsubj_dir`) location to save the matrix `backingfile`.
#' @param error_cutoff : (default = 20) how many missing directories or brain surface files
#' for the function to stop with an error. If < `error_cutoff` directories/files are not
#' found a warning is thrown and missing files are registered in the `issues.log` file.
# #' @param mask : (default = TRUE) only keep cortical vertices, according to FreeSurfer
# #' cortical map.
#' @param save_rds : (default = FALSE) save the supersubject file metadata for re-use
#' in other sessions.
#' @param verbose : (default = TRUE)
#'
#' @return A Filebacked Big Matrix with vertex data for all subjects (dimensions:
#' n_subjects x n_vertices)

#' @import bigstatsr
#'
#' @author Serena Defina, 2024.
#' @export
#'
build_supersubject <- function(subj_dir,
                               folder_ids,
                               supsubj_dir,
                               measure,
                               hemi,
                               n_cores,
                               fwhmc = "fwhm10",
                               fs_template = "fsaverage",
                               backing,
                               error_cutoff = 20,
                               save_rds = FALSE,
                               verbose = TRUE) {

  # Identify brain surface files to load ---------------------------------------
  vw_message(" * retrieving ", length(folder_ids), " brain surface files...",
             verbose = verbose)

  folder_list <- file.path(subj_dir, folder_ids)
  if(!dir.exists(supsubj_dir)) dir.create(supsubj_dir, showWarnings = FALSE)

  log_file <- file.path(supsubj_dir, paste0(hemi,".", measure,".issues.log"))

  # Check that the locations in phenotype exist ----------------
  folders_found <- folder_list[dir.exists(folder_list)]

  if (length(folders_found) < length(folder_list)) {

    folders_not_found <- setdiff(folder_list, folders_found)

    writeLines(c(paste("Attention:", length(folders_not_found),
                       "observations speficied in phenotype were not found:"),
                 gsub(" ", "^", folders_not_found), "\n"), log_file)

    # If many observations are missing, the folder id may be mispecified
    if (length(folders_not_found) > error_cutoff) {
      stop(length(folders_not_found),
      " observations specified in phenotype were not found in `subj_dir`.",
      "\nIs your `folder_id` correctly specified?")
    }

    # If not many observations are missing, notify user but keep at it
    warning(" ! ", length(folders_not_found),
    " observations specified in phenotype were not found in `subj_dir`.",
    "\n   See issues.log file for datails.")
  }

  # FreeSurfer output "w-g.pct" needs special treatment:
  if(measure == "w_g.pct") measure_file <- "w-g.pct" else measure_file <- measure

  mgh_file_name <- paste0(
    hemi, ".", measure_file,
    if (fwhmc != "") { paste0(".", fwhmc)}, ".", fs_template, ".mgh")

  mgh_files <- file.path(folders_found, "surf", mgh_file_name)

  # Check that the files exist ------------------------------
  files_found <- mgh_files[file.exists(mgh_files)]

  if (length(files_found) < length(mgh_files)) {

    files_not_found <- setdiff(mgh_files, files_found)

    writeLines(c(paste("Attention:", length(files_not_found),
                       "files were corrupt or missing:"),
                 files_not_found, "\n"), log_file)

    # If many files are missing, something may have gone wrong in the pre-processing
    if (length(files_not_found) > error_cutoff) {
      stop(length(files_not_found),
      " specified brain surface files were not found in `subj_dir`.",
      "\nAre you sure the FreeSurfer pre-processing ran successfully?")
    } else {
      # If not many files are missing, notify user but keep at it
      warning(" ! ", length(files_not_found),
              " brain surface files were corrupt or missing.",
              "\n   See issues.log file for datails.")
    }

  }

  vw_message(" * ", length(files_found),"/",length(folder_ids),
             " observations found.", verbose = verbose)

  # Build empty large matrix to store all vertex and subjects ------------------

  # Get dimensions
  n_files <- length(files_found) # Number of observations / subjects
  n_verts <- count_vertices(fs_template)

  if (fs_template != 'fsaverage') {
    vw_message(" ! NOTE: downsampling vertices induces (small) registration errors.",
               "\n   This is fine for model tuning but, in the final analysis, ",
               "we reccommend\n   using the high resolution `fsaverage` template.",
               verbose = TRUE)
  }

  # Define backing file for matrix
  if (missing(backing)) {
    backing <- file.path(supsubj_dir,
                         paste(hemi, measure, fs_template, "supersubject.bk",
                               sep="."))
  }
  if (file.exists(backing)) file.remove(backing) # TODO: warn the user

  # Initiate File-backed Big Matrix
  # Dimensions are set so I later access the data column by column and not row by row
  ss <- bigstatsr::FBM(
    nrow = n_files,
    ncol = n_verts,
    type = "float",
    backingfile = gsub(".bk$", "", backing),
    create_bk = !file.exists(backing)
  )

  vw_message(" * populating super-subject matrix...", verbose = verbose)

  # Temporarily disable parallel BLAS (and other) to prevent accidental implicit
  # parallelism
  with_tmp_sysenv(code = {
    failed_to_load <- bigstatsr::big_parallelize(
    X = ss,
    p.FUN = function(X, ind, files_found, n_verts, fs_template) {

      failure_log <- character(0)

      # Process each observation assigned to this worker
      for (i in ind){
        file_path <- files_found[i]

        tryCatch({
          # Write to FBM
          if (fs_template != 'fsaverage') {
            X[i, ] <- load.mgh(file_path)$x[1:n_verts]

          } else { # avoid sub-setting if not necessary
            X[i, ] <- load.mgh(file_path)$x
          }
        }, error = function(e) {
          # Track failed observations
          failure_log <<- c(failure_log, paste(file_path, e$message, '\n'))

          # Fill with NA to maintain shape
          X[i, ] <- NA_real_
        })

      }

      return(failure_log)
    },
    ncores = as.integer(n_cores),
    ind = seq_along(files_found),
    # Pass data to workers
    files_found = files_found,
    n_verts = n_verts,
    fs_template = fs_template,
    p.combine = "c"  # Combine logs from all cores
  )
  })

  # Write combined logs from main process (doing this outside parallel process
  # to avoid race conditions)
  if (length(failed_to_load) > 0) {
    writeLines(c(paste('Error:', length(failed_to_load),
                       'observations failed to load! Replacing with NA.'),
                 failed_to_load), log_file)
  }

  # Save row index names (i.e. observations) to ensure matching
  ss_rownames <- folder_ids[
    vapply(folder_ids,
           function(id) any(grepl(paste0("/", id, "/"), files_found, fixed = TRUE)),
           logical(1))
    ]

  utils::write.table(ss_rownames,
                     file = file.path(supsubj_dir,
                                      paste(hemi, measure, 'ss.rownames.csv',
                                            sep = '.')),
                     row.names = FALSE, col.names = FALSE, quote = FALSE,
                     sep = ",")

  # Save output
  if (save_rds) {
    vw_message(" * saving supersubject matrix to .rds file.", verbose = verbose)
    ss$save()
  }

  return(ss)
}

#' @title
#' Subset an existing supersubject matrix by matching folder IDs
#'
#' @description
#' This function subsets a \link[bigstatsr]{FBM} supersubject matrix that was
#' created using [build_supersubject()], retaining only the rows that
#' correspond to a set of folder (or row) IDs.
#' It reads row names from an associated \code{.csv} file, checks for missing IDs,
#' and writes logs if any are not found. The new subsetted matrix can be saved
#' for future use.
#'
#' @param supsubj_dir Character string indicating the path to the directory
#'   containing the supersubject files (i.e the supersubject matrix itself as a
#'   \code{.rds} file, and the associated \code{.bk} and \code{.rownames.csv}
#'   files).
#' @param supsubj_file Character string indicating the name of the supersubject
#'   \code{.rds} file. Must follow the naming pattern
#'   \code{"<hemi>.<measure>.<fs_template>.supersubject.rds"}.
#' @param folder_ids Character vector of folder IDs to retain in the new ss
#'   matrix. This should also be a column in the phenotype dataset.
#' @param error_cutoff Integer indicating the maximum number of missing IDs that
#'   is allowed before the function throws an error. If the number of missing IDs
#'   is \code{ <= error_cutoff }, a warning is issued instead.
#'   Default: 20.
#' @param new_supsubj_dir Character string indicating the path to the directory
#'   where the new supersubject files should be stored (either temporarily or
#'   permanently if \code{ save_rds == TRUE }. Created if it does not exist.
#' @param n_cores Integer indicating the number of cores to use for parallel
#'   processing.
#' @param save_rds Logical. If \code{TRUE}, the new ss is also saved to a
#'   \code{.rds} file inside \code{new_supsubj_dir}.
#' @param verbose Logical. Default: \code{TRUE}.
#'
#' @return
#' A \link[bigstatsr]{FBM} object containing the subsetted supersubject matrix.
#'
#' @details
#' The function performs the following steps:
#' 1. Reads row names from the supersubject's `.csv` file.
#' 2. Checks whether all `folder_ids` exist in the supersubject.
#' 3. Logs missing IDs to `issues.log` in `new_supsubj_dir`.
#' 4. If the number of missing IDs exceeds `error_cutoff`, stops with an error.
#' 5. Creates a new FBM with only the matching rows, writing it blockwise to
#'    avoid excessive RAM usage.
#' 6. Writes the filtered row names to `ss.rownames.csv` in `new_supsubj_dir`.
#'
# #' @examplesIf file.exists("path/to/original/ss/")
# #' # Subset a supersubject to a small set of IDs
# #' subset_supersubject(
# #'   supsubj_dir = "path/to/original/ss/",
# #'   supsubj_file = "<hemi>.<measure>.<fsaverage>.supersubject.rds",
# #'   folder_ids = pheno_data[, "folder_id"],
# #'   error_cutoff = 20,
# #'   new_supsubj_dir = "path/to/subsetted/ss/",
# #'   save_rds = TRUE
#'
#' @seealso \link{build_supersubject}
#'
#' @export
#'
subset_supersubject <- function(supsubj_dir,
                                supsubj_file,
                                folder_ids,
                                error_cutoff = 20,
                                new_supsubj_dir,
                                n_cores = 1,
                                save_rds = FALSE,
                                verbose = TRUE) {

  rownames_file <- file.path(supsubj_dir,
                             gsub('.fsaverage\\d*\\.supersubject.rds',
                                  '.ss.rownames.csv', supsubj_file))

  ss_rownames <- scan(file = rownames_file, what = character(), sep = "\n",
                      quiet = TRUE)

  ids_not_found <- setdiff(folder_ids, ss_rownames)
  n_ids_not_found <- length(ids_not_found)

  if (n_ids_not_found > 0) {

    log_file <- file.path(new_supsubj_dir,
                          gsub('.fsaverage\\d*\\.supersubject.rds',
                               '.issues.log', supsubj_file))

    writeLines(c(paste("Attention:", n_ids_not_found,
                       "observations speficied in phenotype were not found:"),
                 ids_not_found, "\n"), log_file)

    # If many observations are missing, the folder id may be mispecified
    if (n_ids_not_found > error_cutoff) {
      stop(n_ids_not_found,
           " observations specified in phenotype were not found in `subj_dir`.",
           "\nIs your `folder_id` correctly specified?",
           "\nSee issues.log file for datails.")
    }

    # If not many observations are missing, notify user but keep at it
    warning(" ! ", n_ids_not_found,
            " observations specified in phenotype were not found in `subj_dir`.",
            "\n   See issues.log file for datails.")
  }

  ss <- bigstatsr::big_attach(file.path(supsubj_dir, supsubj_file))

  rows_to_keep <- which(ss_rownames %in% folder_ids)

  if (length(rows_to_keep) < length(ss_rownames)) {
    # We gots to subset

    n_col <- ss$ncol # Number of vertices
    n_row <- length(rows_to_keep) # New number of rows

    # Define backing file for matrix
    if(!dir.exists(new_supsubj_dir)) {
      dir.create(new_supsubj_dir, showWarnings = FALSE)
    }

    backing <- file.path(new_supsubj_dir, gsub(".rds$", ".bk", supsubj_file))
    if (file.exists(backing)) file.remove(backing)

    new_ss <- bigstatsr::FBM(nrow = n_row,
                             ncol = n_col,
                             type = "float",
                             backingfile = gsub(".bk$", "", backing),
                             create_bk = !file.exists(backing))

    vw_message(" * subsetting super-subject matrix...",
               "\n   ", n_row, "/", ss$nrow, " rows matching data[,`folder_id`])",
               verbose = verbose)

    # Fill new ss matrix in chunks
    bigstatsr::big_apply(
        X = new_ss,
        a.FUN = function(X, ind, ss, rows_to_keep) {
          X[ind, ] <- ss[rows_to_keep[ind],]
          invisible(NULL)
        },
        ind = seq_along(rows_to_keep),
        rows_to_keep = rows_to_keep,
        ss = ss,
        block.size = 1000,
        ncores = n_cores
    )

    utils::write.table(ss_rownames[rows_to_keep],
                       file = file.path(new_supsubj_dir,
                                        basename(rownames_file)),
                       row.names = FALSE, col.names = FALSE, quote = FALSE,
                       sep = ",")

    # Save output
    if (save_rds) {
      vw_message(" * saving (subsetted) supersubject matrix to .rds file.",
                 verbose = verbose)
      new_ss$save()
    }

    return(new_ss)

  }
  return(ss)
}
