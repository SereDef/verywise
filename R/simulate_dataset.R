#' @title
#' Simulate longitudinal brain surface dataset with phenotype
#'
#' @description
#' This function generates synthetic longitudinal datasets for
#' multiple sites/cohorts (as well as multiple timepoints/sessions per subject).
#' It generates:
#' \itemize{
#' \item the brain surface data in FreeSurfer format (.mgh files) organised in a
#'  verywise folder structure (see vignettes).
#' \item a mathing \code{pheno} dataframe with participant sex and age mock-data,
#'  saved as "phenotype.csv" file in the \code{path} directory.
#' }
#'
#' @param path Where should the dataset be created.
#' @param overwrite (default = TRUE) whether phenotype file should be overwitten.
#' @param verbose (deafult = TRUE)
#' @param ... Other arguments to be passed to \code{\link{simulate_freesurfer_data}}
#' @inheritParams simulate_freesurfer_data
#'
#' @seealso
#' \code{\link{simulate_freesurfer_data}}, \code{\link{simulate_long_pheno_data}}
#'
#' @author Serena Defina, 2024.
#'
#' @return NULL. Files are written to `path`.
#'
#' @export
#'
simulate_dataset <- function(path,
                             data_structure = list(
                               "cohort1" = list(
                                 "sessions" = c("01", "02", "03"),
                                 "n_subjects" = 10
                               ),
                               "cohort2" = list(
                                 "sessions" = c("01", "02"),
                                 "n_subjects" = 20
                               )
                             ),
                             simulate_association =  NULL,
                             overwrite = TRUE,
                             seed = 31081996,
                             verbose = TRUE,
                             ...) {

  # Create directory in path if it does not exist
  if (!dir.exists(path)) dir.create(path,
                                    showWarnings = FALSE, recursive = TRUE)

  pheno_file_path <- file.path(path, "phenotype.csv")

  if (!overwrite & file.exists(pheno_file_path)) {
    vw_message(" * re-using existing phenotype file created on: ",
               file.info(pheno_file_path)[,"mtime"], verbose = verbose)
    pheno <- utils::read.csv(pheno_file_path)

  } else {
    if (overwrite & file.exists(pheno_file_path)) {
      vw_message(" * overwriting existing phenotype file created on: ",
                 file.info(pheno_file_path)[,"mtime"], verbose = verbose)
    }

    # Simulate phenotype data
    pheno <- simulate_long_pheno_data(
      data_structure = data_structure,
      seed = seed,
      verbose = verbose
    )

    # TODO: make this also a multiple imputation object for testing
    utils::write.csv(pheno,
                     file = pheno_file_path,
                     row.names = FALSE
    )
  }

  if (!is.null(simulate_association)) {
    if (is.numeric(simulate_association)) {
      stopifnot(length(simulate_association) == nrow(pheno))
    } else if (is.character(simulate_association)) {
      parsed_ass <- strsplit(simulate_association, ' ', fixed = TRUE)[[1]]
      beta <- as.numeric(parsed_ass[1])
      var_name <- parsed_ass[3]
      simulate_association <- as.vector(beta * pheno[var_name])[[1]]
    } else {
      stop('`simulate_association` is not correctly specified.')
    }
  }

  # Simulate FreeSurfer data to match
  simulate_freesurfer_data(
    path = path,
    data_structure = data_structure,
    simulate_association = simulate_association,
    seed = seed,
    verbose = verbose,
    ...
  )
}


#' @title
#' Simulate phenotype data
#'
#' @description
#' Simulating phenotype data (in the long format) for multiple cohorts/sites and
#' multiple timepoints/sessions.
#'
#' @inheritParams simulate_dataset
#'
#' @return A dataframe (in long format) with id, sex and age data.
#'
#' @author Serena Defina, 2024.
#'
#' @export
#'
simulate_long_pheno_data <- function(data_structure = list(
                                       "cohort1" = list(
                                         "sessions" = c("01", "02"),
                                         "n_subjects" = 100
                                       ),
                                       "cohort2" = list(
                                         "sessions" = c("01", "02"),
                                         "n_subjects" = 150
                                       )
                                     ),
                                     seed = 31081,
                                     verbose = TRUE) {

  vw_message(" * creating phenotype file...", verbose = verbose)
  set.seed(seed)

  fake_data <- list()

  # Loop through cohorts/sites
  for (l in seq_along(data_structure)) {
    site <- names(data_structure[l])
    n_subjects <- data_structure[[l]]$n_subjects
    sessions <- data_structure[[l]]$sessions

    # Create baseline dataset (session 1)
    baseline <- data.frame(
      id = 1:n_subjects,
      time = sessions[1],
      sex = sample(c("Male", "Female"), n_subjects, replace = TRUE),
      age = round(stats::rnorm(n_subjects, mean = 40, sd = 10), 1)
    )

    if (length(sessions) == 1) {
      fake_data[[site]] <- baseline
    } else {
      t1 <- baseline
      for (s in seq(2, length(sessions))) {
        t2 <- data.frame( # NOTE: removed dependency to dplyr
          id = baseline$id,
          time = sessions[s],
          sex = baseline$sex,
          age = round(baseline$age + stats::rnorm(n_subjects, mean = s - 1, sd = 0.3))
        )
        t1 <- rbind(t1, t2)
      }
      fake_data[[site]] <- t1
    }
  }

  long_df <- dplyr::bind_rows(fake_data, .id = "site")

  long_df["folder_id"] <- file.path(
    long_df$site,
    paste0(
      "sub-", long_df$id,
      "_ses-", long_df$time
    )
  )
  return(long_df)
}


#' @title
#' Simulate FreeSurfer data
#'
#' @description
#' Simulating FreeSurfer data for multiple cohorts/sites and multiple
#' timepoints/sessions. This function emulates the folder structure output obtained
#' by calling FreeSurfer recon_all command. This follows the folder names "sub-", participant ID,
#' "_ses-", session ID.
#'
#' @param path Where should the data be created.
#' @param data_structure A nested list, with top level determining the cohort/
#' dataset/site. Each site is itself a list with two items: \code{"sessions"}: a
#' vector of session names/numbers; and \code{"n_subjects"}: an integer indicating
#' the number of subjects.
#' @param measure (default = "thickness"), measure, used in file names.
#' @param hemi (default = "lh") hemisphere, used in file names.
#' @param fwhmc (default = "fwhm10") full-width half maximum value, used in
#' file names.
#' @param fs_template Character string specifying the FreeSurfer template for
#'   vertex registration. Options:
#'   \itemize{
#'   \item \code{"fsaverage"} (default) = 163842 vertices (highest resolution),
#'   \item \code{"fsaverage6"} = 40962 vertices,
#'   \item \code{"fsaverage5"} = 10242 vertices,
#'   \item \code{"fsaverage4"} = 2562 vertices,
#'   \item \code{"fsaverage3"} = 642 vertices
#'   }
#' @param vw_mean (default = 6.5) mean of the simulated vertex-wise data.
#' @param vw_sd (default = 0.5) standard deviation of the simulated vertex-wise data.
#' @param simulate_association (default = NULL) simulate an association in the
#' format \code{"[beta] * [variable name]"}. For example \code{"0.05 * age"}.
#' This is by default isolated to three regions:
#' the superior temporal gyrus, precentral gyrus and middle temporal gyrus.
#' @param seed (default = 3108) seed used for randomization.
#' @param verbose (dafault = TRUE) verbosity.
#'
#' @author Serena Defina, 2024.
#'
#' @export
#'
simulate_freesurfer_data <- function(path,
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
                                     fs_template = "fsaverage",
                                     vw_mean = 6.5,
                                     vw_sd = 0.5,
                                     simulate_association = NULL,
                                     seed = 3108,
                                     verbose = TRUE) {

  if (hemi == "lh") hemi_name <- "left" else hemi_name <- "right"
  vw_message(" * creating FreeSurfer dataset (", hemi," hemisphere)...",
             verbose = verbose)
  set.seed(seed)

  check_path(path, create_if_not = TRUE)

  total_n_files = sum(sapply(
    data_structure,
    function(site) site$n_subjects * length(site$sessions)))

  vw_message(' * ', total_n_files, ' total files.', verbose = verbose)

  file_counter <- 1

  for (l in seq_along(data_structure)) {
    site <- names(data_structure[l])
    sess <- data_structure[[l]]$sessions
    n <- data_structure[[l]]$n_subjects

    # Create site / cohort folder  --------------------------------------------
    site_dir <- file.path(path, site)
    dir.create(site_dir, showWarnings = FALSE)

    vw_message(" * ", site, ": ", n, " (subjects) x ", length(sess),
      " (sessions) cortical ", measure, " files", verbose = verbose)

    for (i in seq_len(n)) {
      for (t in sess) {
        # Create subject folder structure -------------------------------------
        sub_dir <- file.path(site_dir, paste0("sub-", i, "_ses-", t))

        dir.create(file.path(sub_dir, "surf"), recursive = TRUE, showWarnings = FALSE)
        # Not creating "stats" and "mri" subfolders

        # File name
        mgh_fname <- paste(hemi, measure, fwhmc, fs_template, "mgh", sep = ".")

        n_verts <- switch(fs_template,
                          fsmicro = 10,
                          fsaverage = 163842,
                          fsaverage6 = 40962,
                          fsaverage5 = 10242,
                          fsaverage4 = 2562,
                          fsaverage3 = 642)

        # Start from a vector of type double, filled with 0s
        vw_data <- numeric(n_verts)

        # Randomly generated vertex-wise data

        # First isolate some regions from annotation file in R/sysdata.rda
        # Chose these because they are a small subset (~1.5%) of the surface
        # The rest of the values are left to 0 so they won't be analysed
        locate_verbosity <- if (file_counter == 1) verbose else FALSE
        roi_locs <- locate_roi(rois=c('temporalpole', 'frontalpole',
                                      'entorhinal'), n_verts = n_verts,
                               verbose = locate_verbosity)

        vw_data[roi_locs] <- stats::rnorm(sum(roi_locs),
                                          mean = vw_mean,
                                          sd = vw_sd)
        vw_data[vw_data < 0] <- 0.1 # make sure it is non-negative

        # Specify association cluster
        if (!is.null(simulate_association)) {
          stopifnot(is.numeric(simulate_association) &&
                      (length(simulate_association) == total_n_files))

          # Associations only in my fav region
          assoc_roi <- locate_roi(rois='entorhinal', n_verts = n_verts,
                                  verbose = locate_verbosity)

          vw_data[assoc_roi] <- vw_data[assoc_roi] -
            simulate_association[file_counter]
        }

        file_counter <- file_counter + 1

        save.mgh(as.mgh(vw_data), file.path(sub_dir, "surf", mgh_fname))
        invisible(NULL)
      }
    }
  }
}
