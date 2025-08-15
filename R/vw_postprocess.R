#' @title
#' Convert statistical result FBMs to FreeSurfer `.mgh` format
#'
#' @description
#' This function takes a list of \link[bigstatsr]{FBM} objects storing statistical
#' results and writes them to FreeSurfer-compatible `.mgh` files.
#' It supports coefficient, standard error, p-value, and residual maps, with
#' optional on-the-fly \eqn{-\log_{10}} transformation of p-values.
#'
#' @param vw_results A named list of \link[bigstatsr]{FBM} objects containing the
#'   statistical results.
#' @param result_path Character string indicating the base output path where the
#'   \code{.mgh} files will be written.
#' @param stacks Vector of stack identifiers (e.g. hemisphere stack IDs) to be
#'   included in output filenames (for non-residual stats).
#' @param stat_names Character vector of statistic names to process.
#'   Default: \code{c("coef","se","p","-log10p","resid")}.
#'   The special name \code{"-log10p"} triggers the on-the-fly p-value
#'   transformation.
#' @param verbose Logical. Default: \code{TRUE}
#'
#' @return Invisibly returns \code{NULL}.
#'   Side effects: \code{.mgh} files are written to disk.
#'
#' @details
#' For residuals, all rows are written into a single \code{.mgh} file. For other
#' stats, data is written to one file per row (i.e. term).
#'
#' When computing \eqn{-\log_{10}(p)} values, the transformation is applied
#' **in chunks of columns** to avoid loading the full FBM into memory.
#'
#' @examples
#' \dontrun{
#' library(bigstatsr)
#'
#' # Dummy vw_results list with small FBMs
#' vw_results <- list(
#'   coef = FBM(5, 10, init = rnorm(50)),
#'   se   = FBM(5, 10, init = runif(50, 0.1, 1)),
#'   p    = FBM(5, 10, init = runif(50, 0, 1)),
#'   resid= FBM(20, 10, init = rnorm(200))
#' )
#'
#' convert_to_mgh(
#'   vw_results = vw_results,
#'   result_path = "my_results/stat",
#'   stacks = 1:5
#' )
#' }
#'
#' @export
#'
convert_to_mgh <- function(vw_results,
                           result_path,
                           stacks,
                           stat_names = c("coef", "se", "p", "-log10p", "resid"),
                           verbose = TRUE){

  vw_message("Post-processing\n * converting coefficients, SEs, and p-values",
             " to .mgh format", verbose = verbose)

  lapply(stat_names, function(stat_name) {

    if (stat_name == "resid") {
      stat_mgh_paths <- paste(result_path, "residuals.mgh", sep = ".")
      mode <- "allrows.1file"
      vw_message(" * converting residuals to .mgh format...", verbose = verbose)
    } else {
      stat_mgh_paths <- paste(result_path, paste0("stack", stacks),
                              stat_name, "mgh", sep = ".")
      mode <- "1row.1file"
    }

    # Apply -log10 transformation
    if (stat_name == "-log10p") {
      vw_message(" * (-log10) tranforming pvalues...", verbose = verbose)
      vw_p <- vw_results[["p"]]
      fbm <- bigstatsr::FBM(nrow = vw_p$nrow,
                            ncol = vw_p$ncol,
                            type = vw_p$type_chr,
                            backingfile = gsub(".p.bk$", ".-log10p",
                                               vw_p$backingfile, fixed = TRUE))

      # fbm[] <- -1 * log10(vw_p[])

      # Memory safe transformation
      bigstatsr::big_apply(
        X = vw_p,
        a.FUN = function(X, ind) {
          fbm[, ind] <- -1 * log10(X[, ind])
          invisible(NULL)
        },
        ind = seq_len(vw_p$ncol),
        block.size = 1000
      )

    } else {
      fbm <- vw_results[[stat_name]]
    }

    fbm2mgh(fbm = fbm, fnames = stat_mgh_paths, mode = mode)

    invisible(NULL)
  })
}

#' @title
#' Move key result files from one directory to another (for sharing and visualization)
#'
#' @description
#' The function searches through a \code{verywise} results directory for all files
#' needed for visualization: i.e. clusters, coefficient maps, and "stack names".
#' It then copies them to a destination directory while preserving any
#' sub-directory structure in the source results folder.
#'
#' @param from_dir Character string indicating the path to the "source" results
#'   directory (i.e. \code{outp_dir} in the analysis call)
#' @param to_dir Character string indicating the path to the directory where
#'   matching files will be copied.
#'   Note: any required sub-directories will be created automatically.
#'
#' @details
#' Files are matched using these regular expression patterns:
#' \itemize{
#'   \item \code{"[a-z]{1}h.[a-z]+.stack[0-9]{1,2}.cache.th30.abs.sig.ocn.mgh"}:
#'         cluster files
#'   \item \code{"[a-z]{1}h.[a-z]+.stack[0-9]{1,2}.coef.mgh"}:
#'         coefficient (beta) maps
#'   \item \code{"stack_names.txt"}:
#'         stack names
#' }
#'
#' @return Invisibly returns a logical vector indicating whether each file was
#'   successfully copied.
#'   Side effects: \code{.mgh} files are copied to \code{to_dir}.
#'
#' @examples
#' \dontrun{
#' # Move matching files from "results" to "plots"
#' move_result_files("path/to/results", "path/to/plots")
#' }
#'
#' @export
#'
move_result_files <- function(from_dir, to_dir){

  # Fetch all files inside the folder
  all_files <- list.files(from_dir, recursive = TRUE)

  # Match only necessary files for plotting
  match_files <- c("[a-z]{1}h.[a-z]+.stack[0-9]{1,2}.cache.th30.abs.sig.ocn.mgh", # clusters
                   "[a-z]{1}h.[a-z]+.stack[0-9]{1,2}.coef.mgh", # betas
                   "stack_names.txt") # fixed terms names

  files_to_move <- grep(paste(match_files, collapse="|"),
                        all_files, value = TRUE)

  # Define destination paths, preserving sub-directory structure
  dest_files <- file.path(to_dir, files_to_move)

  # Create all required directories on the destination path
  dirs_needed <- unique(dirname(dest_files))
  for (dir_needed in dirs_needed){
    dir.create(dir_needed, recursive = TRUE, showWarnings = FALSE)
  }

  # Copy the files
  cp_output <- file.copy(file.path(from_dir, files_to_move), dest_files, overwrite = TRUE)

}

