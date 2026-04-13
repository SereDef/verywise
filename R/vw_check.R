check_formula <- function(formula) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("`formula` must be a valid R formula object.")
  }

  # TMP: Assume the brain measure is always the OUTCOME
  measure <- gsub("vw_", "", all.vars(formula)[1])

  check_measure(measure)

  return(measure)
}

check_measure <- function(measure, 
                          fs_metrics = c("thickness", "area", "area.pial", "curv", "jacobian_white", "pial", "pial_lgi", "sulc", "volume",
                                         "w_g.pct", "white.H", "white.K")) {
  
  if (!measure %in% fs_metrics) {
    stop(sprintf("The outcome in `formula` should be a brain surface metric. '%s' is not a valid metric. See `?run_vw_lmm`", measure))
  }

  return(invisible(NULL))
  
}

check_path <- function(dir_path, create_if_not = FALSE) { # file_exists = NULL

  param_name <- deparse(substitute(dir_path))

  if (param_name == "outp_dir" & is.null(dir_path)) {
    outp_dir <- file.path(getwd(), "verywise_results")
    vw_message(" ! WARNING: outpur directory unspecified, which is not recommended.",
                 " You can find the results at ", outp_dir)
    dir.create(outp_dir, showWarnings = FALSE)
    return(outp_dir)
  }

  if (!dir.exists(dir_path)) {
    if (create_if_not) {
      vw_message(sprintf(
        " * the `%s` specified ('%s') does not exist.\n   I'll try to create it.",
        param_name, dir_path))
      dir.create(dir_path, recursive=TRUE)
    } else {
      stop(sprintf("The `%s` specified ('%s') does not exist.", param_name, dir_path))
    }
  }

  # Check for files inside the folder [not used in the latest version]
  # if (!is.null(file_exists)) {
  #   file_exists <- as.vector(file_exists, mode = 'character')
  #   for (f in file_exists) {
  #     file_path <- file.path(dir_path, f)
  #     if (file.exists(file_path)) return(file_path)
  #   }
  # }

  return(dir_path)
}

check_data_list <- function(data_list, folder_id, formula) {

  if (!is.list(data_list) || length(data_list) == 0) {
    stop("`data_list` must be a non-empty list of data.frames.")
  }

  data1 <- data_list[[1]]

  check_data_frame(data1, folder_id, formula)

  # TODO: check that is in "long" format ...?

  # Ensure all data.frames are equivalent to data1
  if (!all(vapply(data_list,
                  function(df) { identical(df[,folder_id], data1[,folder_id]) &&
                      identical(names(df), names(data1)) }, logical(1)))) {
    stop("Each data.frame in `data_list` must have identical dimensions.")
  }

  invisible(NULL)
}

check_data_frame <- function(data, folder_id, formula) {

  stopifnot(inherits(data, "data.frame"))

  if (!(folder_id %in% names(data))) {
    stop(sprintf("Folder ID '%s' not found in data.", folder_id))
  }

  if (any(is.na(data[folder_id]))) {
    stop(sprintf("Folder ID '%s' contains missing values. This is not allowed, please remove any NAs.",
                 folder_id))
  }

  if (any(duplicated(data[folder_id]))) {
    stop(sprintf("Folder IDs must be unique.'%s' contains duplicates.", folder_id))
  }

  terms <- all.vars(formula)[-1] # first is the vw_measure that is not in the pheno
  missing_vars <- !(terms %in% names(data))

  if (any(missing_vars)) {
    stop(sprintf("Variables: '%s' are specified in the formula but not present in the data.",
                 terms[missing_vars]))
  }

  invisible(NULL)
}



check_stack_file <- function(fixed_terms, outp_dir) {

  stack_ids <- data.frame("stack_number" = seq_along(fixed_terms),
                          "stack_name" = as.character(fixed_terms),
                          stringsAsFactors = FALSE)

  stack_file <- file.path(outp_dir, "stack_names.txt")

  if (file.exists(stack_file)) {
    # Read existing file
    existing_stack_ids <- utils::read.table(stack_file,
                                            header = TRUE,
                                            sep = "\t",
                                            stringsAsFactors = FALSE)
    # Compare with new stacks
    if (!identical(existing_stack_ids, stack_ids)) {
      stop(sprintf("A `stack_names.txt` file already exists in `%s`", outp_dir),
           " but its content looks different than expected from your formula. ",
           "This is not allowed to avoid overwriting results by mistake. ",
           "Please use different output directories for different models. ",
           "If overwriting was intentional, please manually delete the existing",
           " `stack_names.txt` file.")
    }
  } else {
    # Write it
    utils::write.table(stack_ids, stack_file, sep = "\t", row.names = FALSE)
  }
}

check_ss_exists <- function(path, ss_file) {

  ss_path <- file.path(path, ss_file)

  if (!file.exists(ss_path)) {
    return(FALSE)
  }

  rn_file <- gsub('.fsaverage\\d*\\.supersubject.rds', '.ss.rownames.csv',
                  ss_file)
  rn_path <- file.path(path, rn_file)

  if (!file.exists(rn_path)) {
    param_name <- deparse(substitute(path))
    stop(sprintf("A `%s` file was found in `%s`, but its `ss.rownames.csv` is missing.",
                 ss_file, param_name))

  }
  return(TRUE)
}

check_weights <- function(weights, data1) {

  if (is.null(weights)) return(invisible(NULL))

  if (is.character(weights) && length(weights) == 1) {
    if (weights %in% names(data1)) return(invisible(NULL))
    stop(sprintf("Weights variable '%s' not found in data.", weights))
  }

  if (is.numeric(weights)) {
    if (length(weights) == nrow(data1)) return(invisible(NULL))
    stop("Weights vector has not the same dimensions as data.")
  }

  stop("The weights format you provided is not valid.")
}

check_cores <- function(n_cores){

  n_cores <- as.integer(n_cores)
  avail_cores <- parallelly::availableCores()
  avail_connections <- parallelly::freeConnections()

  if (n_cores < 1) stop("`n_cores` should be an integer that is 1 or higher.")

  if (n_cores > avail_cores) {
    vw_message(" * WARNING: You requested ", n_cores, " cores but only ",
               avail_cores," are available.\n   Resetting `n_cores` to ",
               avail_cores-1, ".", verbose = TRUE)
    n_cores <- as.integer(avail_cores-1)
  }

  if (n_cores > avail_connections) {
    vw_message(" * WARNING: ", n_cores, " cores excedes R limit of ",
               avail_connections, " free connections for user operations. ",
               "\n   Reducing the number of parallel processes to ",
               avail_connections-1, ".", verbose = TRUE)
    n_cores <- as.integer(min(n_cores, parallelly::freeConnections() - 1))
  }

  return(n_cores)
}

check_freesurfer_setup <- function(FS_HOME, verbose=TRUE) {

  # Set up the necessary FreeSurfer global variables
  if (is.null(FS_HOME) | FS_HOME == "") {

    if (Sys.getenv("FREESURFER_HOME") == "") {
      stop("FREESURFER_HOME needs to be specified or set up.")
    }
    # This is the default already but I recall that in case the user set this
    # after the package was compiled
    FS_HOME <- Sys.getenv("FREESURFER_HOME")

  } else {

    # vw_message("Setting up FreeSurfer environment...", verbose=verbose)
    exit_code <- system(paste("source",
                              file.path(FS_HOME, "SetUpFreeSurfer.sh")))
    if (exit_code != 0) { # --> sourcing failed
      stop(FS_HOME, " is not a FREESURFER_HOME directory.")
    }
    Sys.setenv(FREESURFER_HOME = FS_HOME)
  }

  vw_message("Using FreeSurfer version {.field { basename(FS_HOME) }}", 
             type = 'note', verbose = verbose)

  if (Sys.getenv("SUBJECTS_DIR") == "") {
    Sys.setenv(SUBJECTS_DIR = file.path(FS_HOME, "subjects"))
  }
  # Add $FREESURFER_HOME/bin to $PATH (if not there already) so that FreeSurfer
  # commands can also be called from RStudio
  if (!grepl(file.path(FS_HOME, "bin"), Sys.getenv("PATH"))) {
    Sys.setenv(PATH = paste(Sys.getenv("PATH"),
                            file.path(FS_HOME, "bin"), sep=":"))
  }
  return(invisible(NULL))
}

check_numeric_param <- function(param, lower = -Inf, upper = Inf, 
  integer = FALSE, set = NULL) {
  
  param_name <- deparse(substitute(param))

  if (!is.numeric(param) || length(param) != 1 || is.na(param)) {
    cli::cli_abort("{.strong {param_name}} must be a single numeric value.")
  }
  if (integer && (param %% 1 != 0)) {
    cli::cli_abort("{.strong {param_name}} must be an integer.")
  }
  if (param < lower || param > upper) {
    cli::cli_abort("{.strong {param_name}} must be between:
       {.val {lower}} and {.val {upper}}.")
  }
  if (!is.null(set) && !param %in% set) {
    cli::cli_abort("{.strong {param_name}} must be one of: {.or {.val {set}}}.")
  }
}

check_row_match <- function(ss_file, pheno, folder_ids) {

  rownames_file <- gsub('.fsaverage\\d*\\.supersubject.bk',
                        '.ss.rownames.csv', ss_file)

  ss_rownames <- scan(file = rownames_file,
                      what = character(), sep = "\n", quiet = TRUE)

  if (!identical(ss_rownames, folder_ids)) {

    vw_message(' * matching phenotype with brain surface data')

    # There should be no NA in match(ss_rownames, folder_ids) because
    # there should be no ss_rownames that are not in folder_id
    # (see build_supersubject and subset_supersubject functions)
    matching_key <- match(ss_rownames, folder_ids)
    if (any(is.na(matching_key))) stop("Some rows in ss are not in `folder_id`.")
    
    # Match the row names in ss
    if (inherits(pheno, "data.frame")) {
       pheno <- pheno[matching_key, ]
    } else {
      # folder_id should be identical in all dfs (see check_data_list)
      pheno <- lapply(pheno, function(df) { return(df[matching_key, ])})
    }

    # Note the user was already notified and an error should have already occurred
    # if the drop was above the cutoff
    obs_drop <- length(folder_ids) - length(matching_key)
    if (obs_drop > 0) {
      vw_message('   ', obs_drop, ' observations were dropped from phenotype.')
    }
  }

  return(pheno)
}

# validate_hemi_vw_lmm_inputs <- function(formula, data_list, subj_dir, outp_dir,
#                                         folder_id, hemi, fwhm, mcz_thr, cwp_thr, seed, n_cores) {
#   check_formula(formula)
#   check_data_list(data_list, folder_id)
#   check_directory(subj_dir, "subj_dir")
#   if (!is.null(outp_dir)) check_directory(dirname(outp_dir), "outp_dir parent")
#   check_hemi(hemi)
#   check_numeric_param(fwhm, "fwhm", lower = 0)
#   check_numeric_param(mcz_thr, "mcz_thr", lower = 0)
#   check_numeric_param(cwp_thr, "cwp_thr", lower = 0, upper = 1)
#   check_numeric_param(seed, "seed", integer = TRUE)
#   check_numeric_param(n_cores, "n_cores", lower = 1, integer = TRUE)
# }

check_lambda_grid <- function(lambda_grid) {

  # Auto grid: wide but not absurd; assumes lambda mostly in [1e-6, 1e2]
  if (is.null(lambda_grid)) {
    lambda_grid <- c(0, 10^seq(-6, 2, length.out = 40))
  } else {
    lambda_grid <- as.numeric(lambda_grid)
    if (any(!is.finite(lambda_grid)) || any(lambda_grid < 0)) stop("lambda_grid must be finite and >= 0.")
  }

  lambda_grid
}
