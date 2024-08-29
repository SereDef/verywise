#' @title
#' Simulate a dataset including phenotype and FreeSurfer data
#'
#' @description
#' Simulating a dataset with multiple cohorts/sites and multiple
#' timepoints/sessions. This generates:
#' \itemize{
#' \item a \code{pheno} data.frame with participant sex and age mock-data, saved
#' as "phenotype.csv" file in the \code{path} directory.
#' \item the FreeSurfer data saved as .mgh files organised in a verywise folder
#' structure.
#' }
#'
#' @param path Where should the dataset be created.
#' @param overwrite (default = TRUE) whether phenotype file should be overwitten.
#' @param ... Other arguments to be passed to \code{\link{simulate_freesurfer_data}}
#' @inheritParams simulate_freesurfer_data
#'
#' @seealso
#' \code{\link{simulate_freesurfer_data}}, \code{\link{simulate_long_pheno_data}}
#'
#' @author Serena Defina, 2024.
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
                             simulate_association =  0.05 * pheno$age,
                             overwrite = TRUE,
                             seed = 31081996,
                             ...) {

  # Create directory in path if it does not exist
  if (!dir.exists(path)) dir.create(path)

  pheno_file_path <- file.path(path, "phenotype.csv")

  if (!overwrite & file.exists(pheno_file_path)) {
    message("Re-using existing phenotype file created on: ",
            file.info(pheno_file_path)[,"mtime"])
    pheno <- utils::read.csv(pheno_file_path)

  } else {
    if (overwrite & file.exists(pheno_file_path)) {
      message("Overwriting existing phenotype file created on: ",
            file.info(pheno_file_path)[,"mtime"])
    }

    # Simulate phenotype data
    pheno <- simulate_long_pheno_data(
      data_structure = data_structure,
      seed = seed
    )

    # TODO: make this also a multiple imputation object for testing
    utils::write.csv(pheno,
                     file = pheno_file_path,
                     row.names = FALSE
    )
  }

  # Simulate FreeSurfer data to match
  simulate_freesurfer_data(
    path = path,
    data_structure = data_structure,
    simulate_association = simulate_association,
    seed = seed,
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
#' @inheritParams simulate_freesurfer_data
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
                                     seed = 31081996) {

  cat("Creating phenotype file...\n")
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
#' @param vw_resolution (default = 163842) data dimension (how many vertices)
#' @param measure (default = "thickness"), measure, used in file names.
#' @param hemi (default = "lh") hemisphere, used in file names.
#' @param fwhmc (default = "fwhm10") full-width half maximum value, used in
#' file names.
#' @param target (default = "fsaverage"), used in file names.
#' @param vw_mean (default = 6.5) mean of the simulated vertex-wise data.
#' @param vw_sd (default = 0.5) standard deviation of the simulated vertex-wise data.
#' @param simulate_association (default = NULL) simulate an association in the
#' format \code{beta * variable}. This is by default isolated to three regions:
#' the superior temporal gyrus, precentral gyrus and middle temporal gyrus.
#' @param seed (default = 31081996) seed used for randomization.
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
                                     vw_resolution = 163842,
                                     measure = "thickness",
                                     hemi = "lh",
                                     fwhmc = "fwhm10",
                                     target = "fsaverage",
                                     vw_mean = 6.5,
                                     vw_sd = 0.5,
                                     simulate_association = NULL,
                                     seed = 31081996) {
  cat("Creating FreeSurfer dataset...\n")
  set.seed(seed)

  # TODO: check path exists and create directory if not

  for (l in seq_along(data_structure)) {
    site <- names(data_structure[l])
    sess <- data_structure[[l]]$sessions
    n <- data_structure[[l]]$n_subjects

    # Create site / cohort folder  --------------------------------------------
    site_dir <- file.path(path, site)
    dir.create(site_dir, showWarnings = FALSE)

    cat(
      "Generated", n, "(subjects) x", length(sess), "(sessions) cortical",
      measure, "files in:", site, "\n"
    )

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
        if (!is.null(simulate_association)) {
          # Isolate some regions from annotation file in R/sysdata.rda
          annot <- aparc.annot$vd_label[1:vw_resolution] %in% c(
            14474380, # superior temporal gyrus
            14423100, # precentral gyrus
            3302560 # middle temporal gyrus
          )
          vw_data[annot] <- vw_data[annot] -
            simulate_association[i + (which(sess == t) - 1)]
        }

        save.mgh(as.mgh(vw_data), file.path(sub_dir, "surf", mgh_fname))
        invisible(NULL)
      }
    }
  }
}
