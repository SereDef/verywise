#' @title
#' Vertex-wise Random-Effects Meta-Analysis Across Studies
#'
#' @description
#' This function runs a vertex-wise random-effects meta-analysis across multiple
#' studies (sites or cohorts). Is expects individual study results as they are
#' outputted but `verywise::run_vw_lmm` or by `QDECR`.
#'
#' @param stack : Character or numeric identifier for the fixed term to meta-analyise.
#' @param hemi : Character string, hemisphere to analyze (`"lh"` or `"rh"`).
#' @param result_dirs : Character vector of directories containing the results files of each study.
#' @param study_names : Character vector of study names (must match the length of \code{result_dirs}).
# #' @param study_weights : Numeric vector of study weights (e.g. sample size)
#' @param fs_template : Character string specifying the FreeSurfer template surface.
#' The following values are accepted:
#'  * fsaverage (default) = 163842 vertices (highest resolution),
#'  * fsaverage6 = 40962 vertices,
#'  * fsaverage5 = 10242 vertices,
#'  * fsaverage4 = 2562 vertices,
#'  * fsaverage3 = 642 vertices
#' @param n_cores : Integer, number of CPU cores to use for parallel processing (default: 1).
#' @param verbose : Logical, verbose execution (default: TRUE)
#'
#' @return A named list with three \code{FBM} objects:
#'   \describe{
#'     \item{coef}{Meta-analytic effect size estimates at each vertex.}
#'     \item{se}{Meta-analytic standard errors at each vertex.}
#'     \item{p}{Meta-analytic p-values at each vertex.}
#'   }
#'   Additionally, these are exported as MGH files in the working directory.
#'
#' @details
#' For each vertex, the function loads the effect size and standard error from each study, computes the variance,
#' and runs a random-effects meta-analysis using \code{\link[metafor]{rma}}. Results are saved as
#' Filebacked Big Matrices (FBM) and as MGH files for downstream neuroimaging analysis.
#'
#' @seealso \code{\link[metafor]{rma}}, \code{\link[bigstatsr]{FBM}}
#'
#' @examples
#' \dontrun{
#' run_vw_meta(
#'   stack = 1,
#'   hemi = "lh",
#'   result_dirs = c("study1/results", "study2/results"),
#'   study_names = c("Study1", "Study2"),
#'   fs_template = "fsaverage",
#'   n_cores = 4
#' )
#' }
#'
#' @export
run_vw_meta <- function(stack, hemi, result_dirs,
                        study_names, # study_weights = NULL,
                        fs_template='fsaverage',
                        n_cores = 1, verbose = TRUE) {

  n_studies <- length(study_names)

  if(length(result_dirs) != n_studies) stop('Then number of result directories does not match the number of studies.')

  n_verts <- switch(fs_template,
                    fsaverage = 163842,
                    fsaverage6 = 40962,
                    fsaverage5 = 10242,
                    fsaverage4 = 2562,
                    fsaverage3 = 642,
                    stop("Unknown fs template"))

  # Effect size
  e_vw <- bigstatsr::FBM(n_studies, n_verts, init = 0,
                         backingfile='studies.coef')
  # Variance of effect size
  v_vw <- bigstatsr::FBM(n_studies, n_verts, init = 0,
                         backingfile = 'studies.se')

  for (s in 1:n_studies) {

    dir <- result_dirs[s]

    coef_file <- file.path(dir, paste0(hemi,'.stack',stack,'coef.mgh'))
    se_file <- file.path(dir, paste0(hemi,'stack',stack,'se.mgh'))

    ef_vw[s, ] <- load.mgh(coef_file)
    se_vw[s, ] <- load.mgh(se_file)

  }

  # result_path <- file.path(outp_dir, paste(hemi, measure, sep = "."))


  # Coefficients
  c_vw <- bigstatsr::FBM(1, n_verts, init = 0, backingfile = 'meta.coef')
  # Standard errors
  s_vw <- bigstatsr::FBM(1, n_verts, init = 0, backingfile = 'meta.se')
  # P values
  p_vw <- bigstatsr::FBM(1, n_verts, init = 1, backingfile = 'meta.p')

  if (n_cores > 1 ) vw_message("* preparing cluster of ", n_cores, " workers...",
                               verbose=verbose)
  if (n_cores > 124) {
    vw_message("** ", n_cores, " cores excedes R limit of 125 free connections for user operations, ",
               "reducing the number of parallel processes.", verbose=TRUE)
    n_cores <- min(n_cores, parallelly::freeConnections() - 1)
  }

  # Parallel processing setup
  cluster <- parallel::makeCluster(n_cores,
                                   type = ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK"))
  doParallel::registerDoParallel(cluster)
  on.exit({
    message(pretty_message("All done! :)"))
    parallel::stopCluster(cluster)
  }, add = TRUE)

  v = NULL
  foreach::foreach(v = seq_len(n_verts),
                   .packages = "metafor") %dopar% {

         yi = ef_vw[, v]
         sei = se_vw[, v]

         meta_res <- tryCatch(
           { # Random-effects meta-analysis
           metafor::rma(yi = yi, sei = sei, method = "REML")

           c_vw[, v] <- meta_res$beta[1]
           s_vw[, v] <- meta_res$se
           p_vw[, v] <- meta_res$pval

           },
         error = function(e) {
           # TODO: log error
             c_vw[, v] <- NA_real_
             s_vw[, v] <- NA_real_
             p_vw[, v] <- NA_real_
           })
  }

  out <- list(c_vw, s_vw, p_vw)
  names(out) <- c("coef", "se", "p")

  # Save to MGH files
  # Save model statistics into separate .mgh files

  # convert_to_mgh(out,
  #                result_path,
  #                stacks = seq_along(fixed_terms),
  #                stat_names = c(res_bk_names),
  #                verbose = verbose)
  #
  # for (stat in seq_len(length(out))) {
  #   fbm2mgh(fbm = out[[ stat ]],
  #           fname = paste(hemi, paste0("stack",stack),
  #                         "meta", names(out)[stat], "mgh", sep = "."))
  # }

  return(out)
}
