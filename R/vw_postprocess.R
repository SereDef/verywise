#' @title
#' Convert statistical result FBMs to FreeSurfer `.mgh` format
#'
#' @description
#' This function takes a list of [bigstatsr::FBM] objects storing statistical
#' results and writes them to FreeSurfer-compatible `.mgh` files.
#' It supports coefficient, standard error, p-value, and residual maps, with
#' optional on-the-fly \eqn{-\log_{10}} transformation of p-values.
#'
#' @param vw_results A named list of [bigstatsr::FBM] objects containing the
#'   statistical results.
#' @param result_path Character string indicating the base output path where the
#'   `.mgh` files will be written.
#' @param fixed_terms Vector of fixed term names to be included in output filenames
#'   (as "stack1", "stack2"...).
#' @param random_terms Vector of random term names used to label ICC statistics.
#' @param stat_names Character vector of statistic names to process.
#'   Default: `c("coef", "se","p", "-log10p","resid")`.
#'   The special names: `"-log10p"` and `"fdr"` trigger the on-the-fly p-value 
#'   transformations.
#' @param save_resid Logical. Keep the backing around to then save residuals as fbm. 
#'   Default: `FALSE`.
#' @param verbose Logical. Default:`TRUE`
#'
#' @return Invisibly returns `NULL`. Side effects: `.mgh` files are written to disk.
#'
#' @details
#' For residuals, all rows are written into a single `.mgh` file. For other
#' stats, data is written to one file per row (i.e. term).
#' Note, by default, the p vector is cast down to float (single precision, 32‑bit), 
#' meaning about 7 significant decimal digits are stored (accurately). This is done 
#' to cut disk space needed for the results in half.
#'
#' When computing \eqn{-\log_{10}(p)} values, the transformation is applied
#' **in chunks of columns** to avoid loading the full FBM into memory.
#'
#'
#' @export
#'
convert_to_mgh <- function(vw_results,
                           result_path,
                           fixed_terms = NULL,
                           random_terms = NULL,
                           stat_names = c("coef", "se", "p", "-log10p", "resid"),
                           save_resid = FALSE,
                           verbose = TRUE){
  
  if (verbose) cli::cli_progress_step(
    "Convert coefficients, SEs, and p-values to .mgh format")

  lapply(stat_names, function(stat_name) {

    mode <- "1row.1file"
    rm_backing <- TRUE

    if (stat_name == "resid") {
      stat_mgh_paths <- paste(result_path, "residuals.mgh", sep = ".")
      mode <- "allrows.1file"
      rm_backing <- !save_resid
      
    } else if (stat_name == 'mfit') {
      stat_mgh_paths <- paste(result_path, 
        c('singular', 'aic', 'r2c', 'r2m', paste0('icc', seq_along(random_terms))), "mgh", sep = ".")
      
    } else if (stat_name == "cov" | is.null(fixed_terms)) {
      stat_mgh_paths <- paste(result_path, stat_name, "mgh", sep = ".")

    } else {
      stat_mgh_paths <- paste(result_path, 
        paste0("stack", seq_along(fixed_terms)), stat_name, "mgh", sep = ".")
    }

    # Apply -log10 transformation
    if (stat_name == "-log10p") {

      vw_p <- vw_results[["p"]]
      fbm <- bigstatsr::FBM(nrow = vw_p$nrow,
                            ncol = vw_p$ncol,
                            type = vw_p$type_chr,
                            backingfile = gsub(".p.bk$", ".-log10p",
                                               vw_p$backingfile, fixed = TRUE))

      # Memory safe transformation
      bigstatsr::big_apply(
        X = vw_p,
        a.FUN = function(X, ind) {
          fbm[, ind] <- as.single(-1 * log10(X[, ind]))  # explicit cast to float to avoid loss of precision warning
          invisible(NULL)
        },
        ind = seq_len(vw_p$ncol),
        block.size = 1000
      )

    } else if (stat_name == "fdr") {

     if (verbose) cli::cli_progress_step("FDR correction", spinner=TRUE)

      vw_p <- vw_results[["p"]]
      fbm <- bigstatsr::FBM(nrow = vw_p$nrow,
                            ncol = vw_p$ncol,
                            type = vw_p$type_chr,
                            backingfile = gsub(".p.bk$", ".fdr",
                                               vw_p$backingfile, fixed = TRUE))

      # TMP: not memory safe transformation
      fbm[] <- stats::p.adjust(vw_p[], method = 'BH')

      if (verbose) cli::cli_progress_done()

    } else {
      fbm <- vw_results[[stat_name]]
    }
    # Transform to .mgh
    fbm2mgh(fbm = fbm, fnames = stat_mgh_paths, mode = mode)

    # Remove .bk files
    if (rm_backing) file.remove(fbm$backingfile)

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

