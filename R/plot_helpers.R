#' Validate and load hemisphere data
#'
#' @description
#' Checks if the input is a valid file path or numeric vector. If a `.mgh` file 
#' path is provided, it automatically loads it. Validates that the number of 
#' vertices matches the specified FreeSurfer template.
#'
#' @param hemi A numeric vector, a character string representing a file path 
#'   (e.g., `.mgh`, `.gii`), or `NULL`.
#' @param fs_template A character string specifying the FreeSurfer template 
#'   (e.g., `"fsaverage5"`, `"fsaverage"`).
#'
#' @return A numeric vector representing the surface data, or `NULL` if the input is `NULL`.
#' 
check_hemi <- function(hemi, fs_template) {
  
  if (is.null(hemi)) return(hemi)  # empty or file path: skip check
  
  hemi_name <- deparse(substitute(hemi))
  
  if (is.character(hemi)) {
    if (!file.exists(hemi)) vw_error("{hemi_name} file not found: {.file {hemi}}")
    
    if (grepl('\\.mgh$', hemi, ignore.case = TRUE)) {
      hemi <- load.mgh(hemi)
    } else {
      # nilearn.surface.load_surf_data()
      nilearn_surf_ext <- c("\\.mgz$", "\\.gii$", "\\.nii$",
                            "\\.nii\\.gz$", "\\.npy$", "\\.txt$", "\\.csv$")
    
      if (!any(grepl(paste(nilearn_surf_ext, collapse = "|"), hemi, ignore.case = TRUE)))
          vw_message(c("!" = "{hemi_name}: unrecognised file extension in {.file {basename(hemi)}}.",
                      " " = "nilearn will attempt to load it anyway but things may get weird.",
                      ">" = "try reading it in yourself and providing a vector instead, or using 
                      one of the supported extensions (e.g. .mgh, .mgz, .csv... see `nilearn.surface.load_surf_data()`)"
      ))

      return(hemi)
    }
  }

  hemi <- as.numeric(hemi)

  n_vert <- count_vertices(fs_template)

  if (length(hemi) < n_vert) {
    vw_error("{hemi_name} vector length ({.warn {length(hemi)}}) does not match 
      {fs_template} template ({n_vert}), set `fs_template` to the correct resolution.")
  }

  if (length(hemi) > n_vert) {
      vw_message("!" = "{hemi_name} vector length ({.warn {length(hemi)}}) does not match 
      {fs_template} template ({n_vert}), I will try to subset it.")
  }

  hemi
}

#' Validate a logical surface mask
#'
#' @description
#' Ensures the provided mask is a logical vector and matches the vertex 
#' resolution of the specified FreeSurfer template.
#'
#' @param mask A logical vector used for thresholding, or `NULL`.
#' @param fs_template A character string specifying the FreeSurfer template.
#'
#' @return A logical vector, or `NULL` if the input is `NULL`.
#' 
check_mask <- function(mask, fs_template) {
  
  if (is.null(mask)) return(mask)  # empty or file path: skip check
  
  mask_name <- deparse(substitute(mask))

  if (!is.logical(mask)) {
     vw_error("{.arg {mask_name}} must be a logical vector or NULL.")
  }

  n_vert <- count_vertices(fs_template)

  if (length(mask) != n_vert) {
      vw_message("!" = "{mask_name} vector length ({.warn {length(mask)}}) does not match 
      {fs_template} template ({n_vert}), I will try to subset it.")
  }

  mask
}

#' Resolve FreeSurfer surface mesh location
#'
#' @description
#' Determines whether to use a local FreeSurfer installation for surface meshes 
#' or fall back to `nilearn`'s automatic download and caching mechanism.
#'
#' @param fs_template A character string specifying the template (e.g., `"fsaverage"`).
#' @param fs_home A character string specifying the local `FREESURFER_HOME` path, 
#'   or `NULL` to read from environment variables.
#' @param verbose Logical. Whether to print informational messages. Default `TRUE`.
#'
#' @return The FreeSurfer home path when the template is found locally,
#'   or `NULL` to signal Python to use the `nilearn` download+cache path (via fetch_surf_fsaverage)
#' 
resolve_mesh <- function(fs_template, fs_home, verbose = TRUE) {
  
  if (is.null(fs_home)) {
    fs_home <- Sys.getenv("FREESURFER_HOME")
  }
  
  if (nzchar(fs_home)) {
    surf_dir <- file.path(fs_home, "subjects", fs_template, "surf")
    if (dir.exists(surf_dir)) {
      vw_message("Using local FreeSurfer mesh from {fs_home}", verbose = verbose, type = 'note')
      return(fs_home)
    }
    vw_message("{fs_template} surface mesh not found in {.path $FREESURFER_HOME/subjects/}
         {cli::symbol$arrow_right} falling back to nilearn.")
  }

  # Warn only when the template is not yet cached
  cache <- file.path(path.expand("~"), "nilearn_data", fs_template)
  if (!dir.exists(cache) || length(list.files(cache, recursive = TRUE)) == 0L) {
     vw_message(c("i" = "Will download {fs_template} surface mesh. This may take up to a minute on the first run.",
                  " " = "Meshes will be then cached in {.path ~/nilearn_data/{fs_template}}",
                  " " = "Alternatively, provide a path to FreeSurfer (via `fs_home` argument or environment variables)"))
  }

  NULL 
}

# ===========================================================================================
#' Load and threshold coefficient maps
#'
#' @description
#' Loads a vertex-wise coefficient map from a specified results directory and 
#' generates a significance mask based on a threshold method (e.g., Cluster-Wise 
#' Significance, FDR, uncorrected p-values, or an absolute numeric threshold).
#'
#' @param hemi Hemisphere identifier (`"lh"` or `"rh"`).
#' @param measure The morphological measure (e.g., `"thickness"`, `"area"`).
#' @param stack The stack or contrast identifier.
#' @param res_dir Path to the directory containing the `.mgh` results files.
#' @param threshold The threshold criteria. Can be `"cws"` for cluster-wise 
#'   significance, a string expression like `"fdr < 0.05"` or `"p < 0.01"`, 
#'   a numeric absolute value, or `NULL` for no thresholding.
#'
#' @return A named list containing:
#'   * `coef`: The numeric vector of loaded coefficients.
#'   * `mask`: A logical vector indicating significant vertices, or `NULL`.
#' @export
#' 
load_and_mask_coef <- function(hemi, measure, stack, res_dir, threshold) {

  coef_file <- file.path(res_dir, paste(hemi, measure, stack, "coef.mgh", sep = "."))

  if (!file.exists(coef_file)) {
    vw_message(c("!" = "Coefficient file not found for {hemi}, skipping: {.file {coef_file}}"))
    return(list(coef = NULL, mask = NULL))
  }
  
  coef <- load.mgh(coef_file)

  if (is.null(threshold)) return(list(coef = coef, mask = NULL))
  
  mask <- if (identical(threshold, "cws")) {

    ocn <- .load_threshold_file(type='OCN', res_dir=res_dir, hemi=hemi, measure=measure, stack=stack)
    
    # keep only vertices belonging to a significant cluster (ocn is a
    # positive integer label; non-significant vertices are 0) # coef[ocn==0] <- NA

    if (!is.null(ocn)) ocn != 0 else NULL
  } else if (is.character(threshold)) {

    if (startsWith(threshold, "fdr")) {
      fdr <- .load_threshold_file(type='FDR', res_dir=res_dir, hemi=hemi, measure=measure, stack=stack)
      if (!is.null(fdr)) eval(str2lang(threshold)) else NULL
    } else if (startsWith(threshold, "p")) {
      p <- .load_threshold_file(type='P', res_dir=res_dir, hemi=hemi, measure=measure, stack=stack)
      if (!is.null(p)) eval(str2lang(threshold)) else NULL
    } else {
      vw_message("Invalid threshold specified: {threshold}. We currently support fdr or p-based thresholds. No mask will be applied.")
      NULL
    }
  } else if (is.numeric(threshold)) {
    abs(coef) >= threshold
  } else { 
    vw_message('Invalid threshold specified: {threshold}. No mask will be applied.')
    NULL 
  }

  list(coef = coef, mask = mask)
}

#' Load threshold significance maps
#'
#' Internal helper to locate and load specific threshold `.mgh` files 
#' (OCN, FDR, or P-value maps).
#'
#' @keywords internal
#' @noRd
.load_threshold_file <- function(type, res_dir, hemi, measure, stack) {

  file_type <- switch(type, 
    OCN = "\\.cache\\..*\\.sig\\.ocn\\.mgh$",
    FDR = "\\.fdr.mgh$",
    P = "\\.p.mgh$"
  )

  file <- list.files(res_dir,
      pattern = paste0("^", hemi, "\\.", measure, "\\.", stack, file_type),
      full.names = TRUE)
      
  if (length(file) == 0) {
    vw_message(c("!" = "No {type} file found for {hemi}.", "i" = "Plotting unmasked coefficients."))
    return(NULL)
  } 
      
  if (length(file) > 1) {
    vw_message("!" = "Multiple {type} files found for {hemi}, using: {.file {basename(file[1])}}")
    file <- file[1]
  }

  return(load.mgh(file))
}

# =============================================================================================
#' Compute vertex-wise differences between two maps
#'
#' @description
#' Calculates the absolute difference (`a - b`) between two cortical surface maps 
#' for both the left and right hemispheres. Prints statistical summaries of the 
#' difference distribution.
#'
#' @param lh_a,lh_b Left hemisphere data maps (numeric vectors or file paths) to contrast.
#' @param rh_a,rh_b Right hemisphere data maps (numeric vectors or file paths) to contrast.
#' @param label_a Character string representing the name of the first map (`a`).
#' @param label_b Character string representing the name of the second map (`b`).
#'
#' @return A named list containing `lh` and `rh` numeric vectors representing 
#'   the calculated differences (`a - b`).
#' @export
#' 
vw_diff <- function(lh_a = NULL, lh_b = NULL,
                    rh_a = NULL, rh_b = NULL,
                    label_a = "a", label_b = "b") {
  
  lh_diff <- .hemi_diff(.resolve_map(lh_a), .resolve_map(lh_b), "lh")
  rh_diff <- .hemi_diff(.resolve_map(rh_a), .resolve_map(rh_b), "rh")

  all_diff <- c(lh_diff, rh_diff)
  qs <- stats::quantile(all_diff, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

  min_diff <- qs[[1]]
  max_diff <- qs[[5]]

  stats <- round(
    c(Min = min_diff, `Q1` = qs[[2]], Median = qs[[3]], 
      Mean = mean(all_diff, na.rm = TRUE), `Q3` = qs[[4]], Max = max_diff), 3)

  stat_names <- paste(sprintf("%7s", names(stats)), collapse = "")
  stat_value <- paste(sprintf("%7s", format(stats, trim = TRUE)),collapse = "")

  n_finite <- sum(is.finite(all_diff))
  n_negative <- sum(all_diff < 0, na.rm = TRUE) # b > a
  p_negative <- round(100 * n_negative / n_finite, 1)

  vw_message(c(
    "i" = "Difference map: {label_a} - {label_b} ({n_finite} vertices)",
    " " = "{stat_names}",
    " " = "{stat_value}",
    " " = "\u00a0",
    "*" = "{.strong {label_a} {'<'} {label_b}} for {n_negative} vertices ({p_negative}%)"
  ))

  return(list('lh'=lh_diff, 'rh'=rh_diff))
}

#' Subtract hemisphere maps safely
#'
#' Internal helper to subtract map `b` from `a` with sanity and length checks.
#'
#' @keywords internal
#' @noRd
.hemi_diff <- function(a, b, hemi) {
    if (is.null(a) && is.null(b)) return(NULL)

    if (is.null(a) || is.null(b)) {
      which_null <- if (is.null(a)) "a" else "b"
      vw_error(c("{hemi}_{null_side} is NULL", 
      "i" = "You mush supply both {hemi}_a and {hemi}_b, for a diff to be calculated."))
    }
    
    if (length(a) != length(b)) {
      vw_error(c("{hemi}_a and {hemi}_b have different lengths ({length(a)} vs {length(b)}).",
        "i" = "Both maps must be on the same surface template."
      ))
    }

    a - b
}

#' Resolve map input to numeric vector
#'
#' Internal helper to load an MGH file if a character path is provided, 
#' or return the object directly if it's already a numeric vector.
#'
#' @keywords internal
#' @noRd
.resolve_map <- function(x) {
  if (is.character(x)) return(load.mgh(x))
  if (is.null(x) || is.numeric(x)) return(x)
  vw_error("{.arg {deparse(substitute(x))}} must be a numeric vector or a file path, not {.cls {class(x)}}.")
}
