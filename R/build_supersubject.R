#' @title
#' Build "supersubject" by stacking all vertex data in one large file-backed
#' matrix with dimensions n_subjects x n_vertices.
#'
#' @param subj_dir : path to the FreeSurfer data, this expects a verywise structure.
#' @param folder_id : column in phenotype file that holds the IDs of the observations
#' to include. This also expects a verywise format (e.g. "site1/sub-1_ses-01")
#' @param files_list : (default = \code{\link{list.dirs.till}(subj_dir, n = 2)})
#' allows to specify manually which folders or files to include. You would typically
#' not use this parameter, but rather use \code{folder_id} to specify which observations
#' to include.
#' @param measure : (default = "thickness"), vertex-wise measure, used to identify files.
#' @param hemi : (default = "lh") hemisphere, used to identify files.
#' @param fwhmc : (default = "fwhm10") full-width half maximum value, used to identify files.
#' @param target : (default = "fsaverage"), used to identify files.
#' @param backing : (default = \code{subj_dir}) location to save the matrix \code{backingfile}.
#' @param n_cores : (default = 1) number of cores for parallel processing.
#' @param mask : (dafault = TRUE) only keep cortical vertices, according to FreeSurfer
#' cortical map.
#' @param save_rds : (default = FALSE) save the supersubject file metadata for re-use
#' in other sessions.
#' @param verbose : (default = TRUE)
#'
#' @return A Filebacked Big Matrix with vertex data for all subjects (dimensions:
#' n_subjects x n_vertices)

#' @import bigstatsr
#' @import doParallel
#' @import foreach
#'
#' @author Serena Defina, 2024.
#' @export
#'
build_supersubject <- function(subj_dir,
                               folder_id,
                               files_list = list.dirs.till(subj_dir, n = 2),
                               measure = "thickness",
                               hemi = "lh",
                               fwhmc = "fwhm10",
                               target = "fsaverage",
                               backing = file.path(subj_dir,
                                                   paste0(hemi, ".", measure,
                                                          "supersubject.bk")),
                               n_cores = 1,
                               mask = TRUE,
                               save_rds = FALSE, # dir_tmp,
                               verbose = TRUE) {
  # measure2 <- measure
  # if(measure2 == "w_g.pct") measure2 <- "w-g.pct"

  # Identify list of files to load ---------------------------------------------

  # Select only the ids in phenotype dataframe
  files_list <- files_list[unlist(lapply(folder_id, grep, files_list))]
  # Make sure the order of ids is the same as in the phenotype dataframe
  if (!identical(folder_id, gsub(paste0(subj_dir,'/'), '', files_list))) {
    warning("Participant order is not correct!")
  }

  mgh_file_name <- paste0(hemi, ".", measure,
                          if (fwhmc != "") {paste0(".", fwhmc)}, ".", target, ".mgh")

  mgh_files <- file.path(files_list, "surf", mgh_file_name)

  # Probably do not need this check...?
  # mgh_files_exist <- file.exists(mgh_files)
  # if (!all(mgh_files_exist)) {
  #   id_nonexist <- which(!mgh_files_exist)
  #   length_non <- length(id_nonexist)
  #   if (length_non < 20) stop("The following subjects do not have ",
  #                             mgh_file, ": ",
  #                             paste(files_list[id_nonexist], collapse = ", "))
  #   stop("20 or more subjects do not have the ", mgh_file, " file in their FreeSurfer output.")
  # }

  # Build empty large matrix to store all vertex and subjects ------------------

  # Get dimensions
  n_files <- length(mgh_files) # Number of observations / subjects
  n_verts <- load.mgh(mgh_files[1])$ndim1 # Number of vertices

  # Mask non-cortical vertex data
  if (mask) {
    cortex <- mask_cortex(hemi = hemi, target = target)
    if (length(cortex) != n_verts) stop("Length of cortical mask does not match number of vertices in the data.")
    n_verts <- sum(cortex)
  } else {
    cortex <- rep(TRUE, n_verts)
  }

  file.remove(backing) # TODO: TMP
  # Initiate Filebacked Big Matrix
  # Change this: so i can access the data column by column and not row by row
  ss <- bigstatsr::FBM(nrow = n_files, ncol = n_verts,
                       backingfile = gsub(".bk$", "", backing),
                       create_bk = !file.exists(backing))

  # Set up parallel processing
  # TODO: use bigparallelr instead?
  cl <- if(verbose) parallel::makeForkCluster(n_cores, outfile = "") else parallel::makeForkCluster(n_cores)
  doParallel::registerDoParallel(cl)
  # utils::capture.output(pb <- utils::txtProgressBar(0, n_files, style = 3), file = "/dev/null")
  i <- NULL
  foreach::foreach(i = seq_len(n_files)) %dopar% {
    # utils::setTxtProgressBar(pb, i)
    ss[i,] <- load.mgh(mgh_files[i])$x[cortex] # Populate row with participant info
    NULL # Don't want to return anything
  }
  parallel::stopCluster(cl)

  if (verbose) message("Supersubject object size: ", cat(utils::object.size(ss)))

  # empty_rows <- fbm_row_is_0(ss, n_cores = n_cores)
  # if (any(empty_rows)) message("Additionally removing ", sum(empty_rows), " vertices with constant 0 values.")
  # ss <- ss[cortex & !empty_rows,]

  # Save output
  if (save_rds) {
    message("Saving supersubject matrix to .rds file.")
    ss$save()
  }

  ss
}
