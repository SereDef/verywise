#' @title
#' Build "supersubject" by stacking all vertex data in one matrix
#'
#' @param subj_dir : path to the FreeSurfer data, expect verywise structure
#' @param folder_id : column in phenotype file that holds ids to include in
#' verywise folder format (e.g. "site1/sub-1_ses-01")
#' @param files_list : (default = \code{\link{list.dirs.till}(subj_dir, n = 2)})
#' but allows to specify manually.
#' @param measure : (default = "thickness"), measure, used to identify files.
#' @param hemi : (default = "lh") hemisphere, used to identify files.
#' @param fwhmc : (default = "fwhm10") full-width half maximum value, used to identify files.
#' @param target : (default = "fsaverage"), used to identify files.
#' @param backing : (default = \code{subj_dir}) location to save the matrix \code{backingfile}.
#' @param n_cores : (default = 1) number of cores for parallel processing.
#' @param verbose : (default = TRUE)
#'
#' @return A Filebacked Big Matrix with vertex data for all subjects (dimensions:
#' n_vertex * n_subject)

#' @importFrom bigstatsr FBM
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
                                                          "supersubject")),
                               n_cores = 1,
                               # mask, dir_tmp,
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
  n_files <- length(mgh_files)
  temp_mgh <- load.mgh(mgh_files[1])

  # Initiate Filebacked Big Matrix
  bigM <- bigstatsr::FBM(nrow = temp_mgh$ndim1, ncol = n_files, backingfile = backing)

  # Set up parallel processing
  # TODO: use bigparallelr instead?
  cl <- if(verbose) parallel::makeForkCluster(n_cores, outfile = "") else parallel::makeForkCluster(n_cores)
  doParallel::registerDoParallel(cl)
  utils::capture.output(pb <- utils::txtProgressBar(0, n_files, style = 3), file = "/dev/null")
  i <- NULL
  foreach::foreach(i = seq_len(n_files)) %dopar% {
    utils::setTxtProgressBar(pb, i)
    bigM[,i] <- load.mgh(mgh_files[i])$x
    NULL
  }
  parallel::stopCluster(cl)
  # TODO: mask cortex?
  # if (length(mask) != nrow(bigM)) stop("Length of mask does not equal number of vertices from data.")
  bigM
}
