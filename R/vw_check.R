check_formula <- function(formula) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("`formula` must be a valid R formula object.")
  }

  # TMP: Assume the brain measure is always the OUTCOME
  measure <- gsub("vw_", "", all.vars(formula)[1])

  fs_metrics <- c("thickness", "area", "area.pial", "curv",
                  "jacobian_white", "pial", "pial_lgi", "sulc", "volume",
                  "w_g.pct", "white.H", "white.K")

  if (!measure %in% fs_metrics) {
    stop(sprintf("The outcome in `formula` should be a brain surface metric.",
                 "'%s' is not a valid metric. See `?run_vw_lmm`", measure))
  }

  return(measure)
}

check_data_list <- function(data_list, folder_id, formula) {

  if (!is.list(data_list) || length(data_list) == 0) {
    stop("`data_list` must be a non-empty list of data.frames.")
  }

  data1 <- data_list[[1]]

  if (!all(vapply(data_list,
                  function(df) { identical(dim(df), dim(data1)) &&
                      identical(names(df), names(data1)) }, logical(1)))) {
    stop("Each data.frame in `data_list` must have the same dimentions.")
  }

  if (!(folder_id %in% names(data1))) {
    stop(sprintf("Folder ID '%s' not found in data.", folder_id))
  }

  if (any(is.na(data1[folder_id]))) {
    stop(sprintf("Folder ID '%s' contains missing values. ",
                 "This is not allowed, please remove any NAs.", folder_id))
  }

  if (any(duplicated(data1[folder_id]))) {
    stop(sprintf("Folder IDs must be unique.'%s' contains duplicates.", folder_id))
  }

  terms <- all.vars(formula)[-1] # first is the vw_measure that is not in the pheno
  missing_vars <- !(terms %in% names(data1))

  if (any(missing_vars)) {
    stop(sprintf("Variables: '%s' are specified in the formula but not present in the data.",
                 terms[missing_vars]))
  }
  # TODO: check that is in long format ?
}

check_path <- function(dir_path, create_if_not=FALSE, file_exists=NULL) {
  param_name <- deparse(substitute(dir_path))

  if (!dir.exists(dir_path)) {
    if (create_if_not) {
      vw_message(sprintf("The `%s` specified ('%s') does not exist. I'll try to create it.",
                         param_name, dir_path))
      dir.create(dir_path, recursive=TRUE)
    } else {
      stop(sprintf("The `%s` specified ('%s') does not exist.", param_name, dir_path))
    }
  }

  if (!is.null(file_exists)) {
    file_path <- file.path(dir_path, file_exists)
    if (file.exists(file_path)) {
    warning(sprintf("A `%s` file already exists inside '%s'. Be careful, you may be overwriting results.",
                    file_exists, dir_path))
    }
  }
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


check_weights <- function(weights, data_list) {

  if (is.null(weights)) return(invisible(NULL))

  data1 <- data_list[[1]]

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

  if (n_cores < 1) stop("`n_cores` shoudl be an integer that is 1 or higher.")

  if (n_cores > avail_cores) {
    vw_message("WARNING: You requested ", n_cores, " cores but only ", avail_cores,
               " are available. Resetting `n_cores` to ", avail_cores-1, ".",
               verbose=TRUE)
    n_cores <- avail_cores-1
  }

  if (n_cores > avail_connections) {
    vw_message("WARNING: ", n_cores, " cores excedes R limit of ",
               avail_connections, " free connections for user operations. ",
               "Reducing the number of parallel processes to ", avail_connections-1,
               ".", verbose=TRUE)
    n_cores <- min(n_cores, parallelly::freeConnections() - 1)
  }

  return(n_cores)
}



check_freesurfer_setup <- function(FS_HOME, verbose=TRUE) {
  # Set up the necessary FreeSurfer global variables
  if (Sys.getenv("FREESURFER_HOME") == "") {

    if (is.null(FS_HOME) | FS_HOME == "") stop("FREESURFER_HOME needs to be specified or set up.")

    vw_message("Setting up FreeSurfer evironment...", verbose=verbose)
    Sys.setenv(FREESURFER_HOME = FS_HOME)
    system(paste("source", file.path(FS_HOME,"SetUpFreeSurfer.sh")))
  }

  if (Sys.getenv("SUBJECTS_DIR") == "") Sys.setenv(SUBJECTS_DIR = file.path(FS_HOME,'subjects'))

  # Add $FREESURFER_HOME/bin to $PATH (if not there already) so that FreeSurfer
  # commands can also be called from RStudio
  if (!grepl(file.path(FS_HOME, "bin"), Sys.getenv("PATH"))) {
    Sys.setenv(PATH=paste(Sys.getenv("PATH"), file.path(FS_HOME,"bin"), sep=":"))
  }

}

check_numeric_param <- function(param, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(param) || length(param) != 1 || is.na(param)) {
    stop(sprintf("`%s` must be a single numeric value.", name))
  }
  if (integer && (param %% 1 != 0)) {
    stop(sprintf("`%s` must be an integer.", name))
  }
  if (param < lower || param > upper) {
    stop(sprintf("`%s` must be between %s and %s.", name, lower, upper))
  }
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
