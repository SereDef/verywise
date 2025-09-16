#' @title
#' Simulate a longitudinal brain surface dataset with associated phenotype data
#'
#' @description
#' Generates a synthetic longitudinal dataset for multiple sites/cohorts,
#' each with multiple timepoints/sessions per subject. The function produces:
#'
#' \itemize{
#'   \item Brain surface data in FreeSurfer \code{.mgh} format, organised
#'   in a verywise folder structure (see vignettes for details).
#'   \item A matching \code{pheno} data frame with mock participant sex and
#'   age, saved as \file{"phenotype.csv"} in the \code{path} directory.
#' }
#'
#' This is useful for testing pipelines or demonstrations where realistic
#' FreeSurfer-style data and phenotypic information are required.
#'
#' @param path Character string. Directory where the dataset should be created.
#'   Will be created if it does not exist.
#' @param data_structure Named list defining cohorts/sites. Each element is a
#'   list with:
#'   \describe{
#'     \item{\code{"sessions"}}{Character vector of session labels.}
#'     \item{\code{"n_subjects"}}{Integer number of subjects.}
#'   }
#' @param fs_template Character (default = \code{"fsaverage"}). FreeSurfer
#'   template for vertex registration. This is used to determine the size of
#'   the synthetic brain surface data. Options:
#'   \itemize{
#'     \item \code{"fsaverage"} = 163842 vertices (highest resolution)
#'     \item \code{"fsaverage6"} = 40962 vertices
#'     \item \code{"fsaverage5"} = 10242 vertices
#'     \item \code{"fsaverage4"} = 2562 vertices
#'     \item \code{"fsaverage3"} = 642 vertices
#'   }
#' @param simulate_association Optional. If numeric, must be of length equal to
#'   the number of generated files; if character, must have the format
#'   \code{"<beta> * <variable_name>"}. Associations are injected into one small
#'   region (the entorhinal cortex).
#' @param overwrite Logical (default = \code{TRUE}). Whether to overwrite an
#'   existing phenotype file.
#' @param seed Integer (default = \code{3108}). Random seed.
#' @param verbose Logical (default = \code{TRUE}). If \code{TRUE}, print
#'   progress messages.
#' @param ... Additional arguments passed to \code{\link{simulate_freesurfer_data}}.
#'
#'
#' @seealso
#' \code{\link{simulate_freesurfer_data}},
#' \code{\link{simulate_long_pheno_data}}
#'
#' @author Serena Defina, 2024.
#'
#' @return
#' Invisibly returns \code{NULL}. Data and phenotype files are written to
#' \code{path}.
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
                             fs_template = "fsaverage",
                             simulate_association =  NULL,
                             overwrite = TRUE,
                             seed = 3108,
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
      var_z <- scale(pheno[, parsed_ass[3]])
      simulate_association <- as.vector(beta * var_z)
    } else {
      stop('`simulate_association` is not correctly specified.')
    }
  }

  # Simulate FreeSurfer data to match
  simulate_freesurfer_data(
    path = path,
    data_structure = data_structure,
    fs_template = fs_template,
    simulate_association = simulate_association,
    seed = seed,
    verbose = verbose,
    ...
  )
}


#' @title
#' Simulate (longitudinal) phenotype data
#'
#' @description
#' Generates synthetic phenotype data in long format for multiple cohorts/sites
#' and multiple timepoints/sessions per subject. Each record contains:
#' \itemize{
#'   \item Subject ID
#'   \item Session/timepoint
#'   \item Sex
#'   \item Age
#'   \item A \code{folder_id} field matching the FreeSurfer directory structure
#' }
#'
#' @inheritParams simulate_dataset
#'
#' @return
#' A \code{data.frame} in long format with columns:
#' \code{site}, \code{id}, \code{time}, \code{sex}, \code{age}, \code{folder_id}.
#'
#' @seealso
#' \code{\link{simulate_dataset}}, \code{\link{simulate_freesurfer_data}}
#'
#' @author
#' Serena Defina, 2024.
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
                                     seed = 3108,
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
      age = round(stats::rnorm(n_subjects, mean = 40, sd = 1.5), 1)
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
          age = round(baseline$age + stats::rnorm(n_subjects, mean = s+1, sd = 0.5))
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
#' Simulate longitudinal FreeSurfer vertex-wise data
#'
#' @description
#' Simulates FreeSurfer-formatted brain surface data for multiple cohorts/sites
#' and multiple timepoints/sessions. The output folder structure emulates that
#' produced by the FreeSurfer \code{recon-all} command, with directories named:
#' \code{sub-<ID>_ses-<SESSION>}.
#'
#' The vertex-wise data are stored as \code{.mgh} files, with optional simulation
#' of an association with a phenotype variable.
#'
#' @inheritParams simulate_dataset
#' @param measure Character (default = \code{"thickness"}). Surface measure to
#'   simulate. This is used for file names.
#' @param hemi Character (default = \code{"lh"}). Hemisphere code: \code{"lh"}
#'   (left) or \code{"rh"} (right). This is used for file names.
#' @param fwhmc Character (default = \code{"fwhm10"}). Full-width half maximum
#'   smoothing parameter. This is used for file names.
#' @param vw_mean Numeric (default = \code{6.5}). Mean of the simulated
#'   vertex-wise values.
#' @param vw_sd Numeric (default = \code{0.5}). Standard deviation of the
#'   simulated vertex-wise values.
#'
#' @author Serena Defina, 2024.

#' @return
#' Invisibly returns \code{NULL}. Vertex-wise data are written to \code{path}
#' in FreeSurfer-compatible format.
#'
#' @seealso
#' \code{\link{simulate_dataset}}, \code{\link{simulate_long_pheno_data}}
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
                                     fs_template = "fsaverage",
                                     measure = "thickness",
                                     hemi = "lh",
                                     fwhmc = "fwhm10",
                                     vw_mean = 6.5,
                                     vw_sd = 0.5,
                                     simulate_association = NULL,
                                     seed = 3108,
                                     verbose = TRUE) {

  hemi_name <- if (hemi == "lh") "left" else "right"
  vw_message(" * creating FreeSurfer dataset (", hemi_name," hemisphere)...",
             verbose = verbose)
  set.seed(seed)

  check_path(path, create_if_not = TRUE)

  total_n_files = sum(sapply(
    data_structure,
    function(site) site$n_subjects * length(site$sessions)))

  vw_message(' * ', total_n_files, ' total files.', verbose = verbose)

  file_counter <- 1

  n_verts <- count_vertices(fs_template)

  # First isolate some regions from annotation file in R/sysdata.rda
  # Chose these because they are a small subset (~1.5%) of the surface
  # The rest of the values are left to 0 so they won't be analysed
  roi_locs <- locate_roi(rois=c('temporalpole', 'frontalpole',
                                'entorhinal'), n_verts = n_verts,
                         verbose = verbose)
  # Associations only in my fav region
  assoc_roi <- locate_roi(rois='entorhinal', n_verts = n_verts,
                          verbose = verbose)

  # File name
  mgh_fname <- paste(hemi, measure, fwhmc, fs_template, "mgh", sep = ".")

  # Start from a vector of type double, filled with 0s
  vw_data <- numeric(n_verts)

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

        # Randomly generated vertex-wise data
        vw_data[roi_locs] <- stats::rnorm(sum(roi_locs),
                                          mean = vw_mean,
                                          sd = vw_sd)
        vw_data[vw_data < 0] <- 0.1 # make sure it is non-negative

        # Specify association cluster
        if (!is.null(simulate_association)) {
          stopifnot(is.numeric(simulate_association) &&
                      (length(simulate_association) == total_n_files))

          vw_data[assoc_roi] <- vw_mean + simulate_association[file_counter]
        }

        file_counter <- file_counter + 1

        save.mgh(as.mgh(vw_data), file.path(sub_dir, "surf", mgh_fname))
        invisible(NULL)
      }
    }
  }
}
