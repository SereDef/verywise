check_formula <- function(formula) {
  if (missing(formula) || !inherits(formula, "formula")) {
    stop("`formula` must be a valid R formula object.")
  }
}

check_data_list <- function(data_list, folder_id) {

  if (!is.list(data_list) || length(data_list) == 0) {
    stop("`data_list` must be a non-empty list of data.frames.")
  }

  if (!(folder_id %in% names(data_list[[1]]))) {
    stop(sprintf("Folder ID '%s' not found in data.", folder_id))
  }

  # TODO:
  # Check there are no duplicates in folder_id
  # check that is in long format
  # check the the formula are columns in dataset
}

check_path <- function(dir_path) {
  param_name <- deparse(substitute(dir_path))
  if (!dir.exists(dir_path)) {
    stop(sprintf("The `%s` specified ('%s') does not exist.", param_name, dir_path))
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

# check_formula(formula)
# check_data_list(data_list, folder_id)
# check_directory(subj_dir, "subj_dir")
# if (!is.null(outp_dir)) check_directory(dirname(outp_dir), "outp_dir parent")
# check_hemi(hemi)
# check_numeric_param(fwhm, "fwhm", lower = 0)
# check_numeric_param(mcz_thr, "mcz_thr", lower = 0)
# check_numeric_param(cwp_thr, "cwp_thr", lower = 0, upper = 1)
# check_numeric_param(seed, "seed", integer = TRUE)
# check_numeric_param(n_cores, "n_cores", lower = 1, integer = TRUE)
#
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
