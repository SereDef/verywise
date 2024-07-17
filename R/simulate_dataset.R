
#' @title
#' Simulate FreeSurfer data
#'
#' @description
#' Simulating FreeSurfer data for multiple cohorts/sites and multiple
#' timepoints/sessions. The funtions emulated the folder structure output obtained
#' by calling recon_all. This follows the folder names "sub-", participant ID,
#' "_ses-", session ID.
#'
#' @param path : Where should the data be created.
#' @param vw_resolution : (default = 163842) data dimension (how many vertices)
#' @param data_structure : A nested list, with top level determining the cohort/
#' dataset/site. Each site is itself a list with two items: \code{"sessions"}: a
#' vector of session names/numbers; and \code{"n_subjects"}: aninteger indicating
#' the number of subjects.
#' @param measure : (default = "thickness"), measure, used in file names.
#' @param hemi : (default = "lh") hemisphere, used in file names.
#' @param fwhmc : (default = "fwhm10") full-width half maximum value, used in
#' file names.
#' @param target : (default = "fsaverage"), used in file names.
#' @param vw_mean : (default = 6.5) mean of the simulated vertex-wise data.
#' @param vw_sd : (default = 0.5) standard deviation of the simulated vertex-wise data.
#' @param seed : (default = 31081996) seed used for randomization.
#'
#' @importFrom stats rnorm
#' @export
#'
simulate_freesurfer_data <- function(path,
                                     vw_resolution = 163842,
                                     data_structure = list(
                                       "cohort1" = list(
                                         "sessions" = c("01", "02"),
                                         "n_subjects" = 100
                                       ),
                                       "cohort2" = list(
                                         "sessions" = c("01", "02"),
                                         "n_subjects" = 150
                                       )
                                     ),
                                     measure = "thickness",
                                     hemi = "lh",
                                     fwhmc = "fwhm10",
                                     target = "fsaverage",
                                     vw_mean = 6.5,
                                     vw_sd = 0.5,
                                     seed = 31081996) {
  cat("Creating dataset...\n")
  set.seed(seed)

  # TODO: check path exists and create directory if not

  for (l in seq_along(data_structure)) {
    site <- names(data_structure[l])
    sess <- data_structure[[l]]$sessions
    n <- data_structure[[l]]$n_subjects

    cat("Cohort/site:", site, "\n")

    # Create site / cohort folder  --------------------------------------------
    site_dir <- file.path(path, site)
    dir.create(site_dir, showWarnings = FALSE)

    cat("Generating", n, "by", length(sess), "simulated cortical",
        measure, "datasets\n")

    for (i in seq_len(n)) {
      for (t in sess) {
        # Create subject folder structure -------------------------------------
        sub_dir <- file.path(site_dir, paste0("sub-", i, "_ses-", t))

        dir.create(sub_dir, showWarnings = FALSE)
        dir.create(file.path(sub_dir, "surf"), showWarnings = FALSE)
        dir.create(file.path(sub_dir, "stats"), showWarnings = FALSE)
        dir.create(file.path(sub_dir, "mri"), showWarnings = FALSE)

        # File name
        mgh_fname <- paste(hemi, measure, fwhmc, target, "mgh", sep = ".")

        # Randomly generated vertex-wise data
        vw_data <- stats::rnorm(vw_resolution, mean = vw_mean, sd = vw_sd)
        vw_data[vw_data <= 0] <- 0.1 # non-negative

        # Specify association cluster
        # aparc.annot <- readRDS("data/aparc.annot.rds") # load in the aparc annotation
        # isolate supe temp gyrus, precentral gyrus, middle temp gyrus
        # annot <- aparc.annot$vd_label %in% c(14474380, 14423100, 3302560)
        # data[annot] <- data[annot] - 0.03 * pheno$age[i]

        save.mgh(as.mgh(vw_data), file.path(sub_dir, "surf", mgh_fname))
        invisible(NULL)
      } # if (!(i%%100)) cat("Now at suject", i, "\n")
    }
  }
}
