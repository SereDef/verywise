#' @title
#' Build "supersubject" by stacking all vertex data in one large file-backed
#' matrix with dimensions n_subjects x n_vertices.
#'
#' @param subj_dir : path to the FreeSurfer data, this expects a verywise structure.
#' @param folder_ids : the vector of observations to include. This holds the relative
#' path (from `subj_dir`) to the FreeSurfer data folder (e.g. "site1/sub-1_ses-01").
#' @param outp_dir : output path, where logs, backing files and the matrix itself
#' (if `save_rds == TRUE`) will be stored.
#' @param measure : vertex-wise measure, used to identify files.
#' @param hemi : hemisphere, used to identify files.
#' @param n_cores : number of cores for the parallel cluster
#' @param fwhmc : (default = "fwhm10") full-width half maximum value, used to identify files.
#' @param target : (default = "fsaverage"), used to identify files.
#' @param backing : (default = `outp_dir`) location to save the matrix \code{backingfile}.
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
                               outp_dir,
                               measure,
                               hemi,
                               n_cores,
                               fwhmc = "fwhm10",
                               target = "fsaverage",
                               backing,
                               error_cutoff = 20,
                               # mask = TRUE,
                               save_rds = FALSE, # dir_tmp,
                               verbose = TRUE) {

  # TODO: check measure names! measure2 <- measure
  # if(measure2 == "w_g.pct") measure2 <- "w-g.pct"

  # Identify brain surface files to load ---------------------------------------
  vw_message("Retrieving ", length(folder_ids), " brain surface files...", verbose=verbose)

  # files_list = list.dirs.till(subj_dir, n = 2)
  # files_list <- files_list[unlist(lapply(folder_ids, grep, files_list))]
  folder_list <- file.path(subj_dir, folder_ids)
  log_file <- file.path(outp_dir, paste0(hemi,".", measure,".issues.log"))

  # Check that the locations in phenotype exist ----------------
  folders_found <- folder_list[dir.exists(folder_list)]

  if (length(folders_found) < length(folder_list)) {

    folders_not_found <- setdiff(folder_list, folders_found)

    # If many observations are missing, the folder id may be mispecified
    if (length(folders_not_found) > error_cutoff) {
      stop(length(folders_not_found), " observations specified in phenotype were
           not found in `subj_dir`. Is your `folder_id` correctly specified?")
    }

    # If not many observations are missing, notify user but keep at it
    warning(length(folders_not_found), " observations specified in phenotype were
            not found in `subj_dir`. See .log file for datails.")

    writeLines(c(paste("Attention:",length(folders_not_found),
                       "observations speficied in phenotype were not found:"),
                       folders_not_found, "\n"), log_file)
  }

  mgh_file_name <- paste0(
    hemi, ".", measure,
    if (fwhmc != "") { paste0(".", fwhmc)},
    ".fsaverage.mgh" # TODO: handle other low-resolution registrations?
  )

  mgh_files <- file.path(folders_found, "surf", mgh_file_name)

  # Check that the files exist ------------------------------
  files_found <- mgh_files[file.exists(mgh_files)]

  if (length(files_found) < length(mgh_files)) {

    files_not_found <- setdiff(mgh_files, files_found)

    # If many files are missing, something may have gone wrong in the pre-processing
    if (length(files_not_found) > error_cutoff) {
      stop(length(files_not_found), " specified brain surface files were not found in `subj_dir`.
           Are you sure the FreeSurfer pre-processing ran successfully?")
    }

    # If not many files are missing, notify user but keep at it
    warning(length(files_not_found), " brain surface files were corrupt or missing.
            See .log file for datails.")

    writeLines(c(paste("Attention:", length(files_not_found),
                       "files were corrupt or missing:"),
                 files_not_found, "\n"), log_file)
  }

  vw_message(length(files_found),"/",length(folder_ids)," observations found.", verbose=verbose)

  # Build empty large matrix to store all vertex and subjects ------------------

  # Get dimensions
  n_files <- length(files_found) # Number of observations / subjects
  # n_verts <- load.mgh(files_found[1])$ndim1 # Number of vertices
  n_verts <- switch(target,
                    fsaverage = 163842,
                    fsaverage6 = 40962,
                    fsaverage5 = 10242,
                    fsaverage4 = 2562,
                    fsaverage3 = 642)

  if (target != 'fsaverage') {
    vw_message('Warining: downsampling vertices induces (small) registration errors.',
               'This is fine for model tuning but, in the final analysis,',
               'we reccommend using the high resolution `fsaverage` template.',
               verbose=TRUE)
  }

  # Mask non-cortical vertex data
  # vw_message("Applying cortical mask...", verbose=verbose)
  # if (mask) {
  #   is_cortex <- mask_cortex(hemi = hemi, target = target)
  #   if (length(is_cortex) != n_verts) stop("Length of cortical mask does not match number of vertices in the data.")
  #   n_verts <- sum(is_cortex)
  #   vw_message(n_verts, '/',length(is_cortex), 'vertices in cortex.', verbose=verbose)
  # } else {
  #   is_cortex <- rep(TRUE, n_verts)
  # }

  # Define backing file for matrix
  if (missing(backing)) {
    backing <- file.path(outp_dir, paste(hemi, measure, target,
                                         "supersubject.bk", sep="."))
  }
  if (file.exists(backing)) file.remove(backing) # TODO: warn the user

  vw_message("Building super-subject matrix...", verbose=verbose)

  # Initiate File-backed Big Matrix
  # Dimensions are set so I later access the data column by column and not row by row
  ss <- bigstatsr::FBM(
    nrow = n_files,
    ncol = n_verts,
    backingfile = gsub(".bk$", "", backing),
    create_bk = !file.exists(backing)
  )

  failed_to_load <- bigstatsr::big_parallelize(
    X = ss,
    p.FUN = function(X, ind, files_found, n_verts) {

      failure_log <- character(0)

      # Process each observation assigned to this worker
      file_path <- files_found[ind]

      tryCatch({
        # Write to FBM (+1 for vertex index row)
        if (target != 'fsaverage') {

          X[ind, ] <- load.mgh(file_path)$x[1:n_verts]

        } else {

          X[ind, ] <- load.mgh(file_path)$x }

      }, error = function(e) {
        # Track failed observations
        failure_log <<- c(failure_log, file_path)

        # Fill with NA to maintain shape
        X[ind, ] <- NA_real_
      })

      return(failure_log)
    },
    ncores = n_cores,
    ind = seq_along(files_found),
    # Pass data to workers
    files_found = files_found,
    # is_cortex = is_cortex,
    n_verts = n_verts,
    combine = "c"  # Combine logs from all cores
  )

  # Write combined logs from main process (doing this outside parallel process
  # to avoid race conditions)
  if (length(failed_to_load) > 0) {
    writeLines(c(paste('Error:', length(failed_to_load),
                       'observations failed to load! Replacing with NA.'),
                 failed_to_load), log_file)
  }

  # Save output
  if (save_rds) {
    vw_message("Saving supersubject matrix to .rds file.", verbose=verbose)
    ss$save()
  }

  return(ss)
}
