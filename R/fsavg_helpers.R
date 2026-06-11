#' Count Vertices for a FreeSurfer Template Surface
#'
#' Returns the number of vertices per hemisphere for a given FreeSurfer
#' fsaverage surface template. Templates are icosphere subdivisions of the
#' fsaverage standard surface, with vertex counts decreasing at lower
#' resolutions.
#'
#' @param fs_template A character string specifying the FreeSurfer template.
#'   Must be one of:
#'   \describe{
#'     \item{\code{"fsaverage"}}{Full-resolution standard template (163,842 vertices per hemisphere).}
#'     \item{\code{"fsaverage6"}}{Medium-resolution template (40,962 vertices per hemisphere).}
#'     \item{\code{"fsaverage5"}}{Low-resolution template (10,242 vertices per hemisphere).}
#'     \item{\code{"fsaverage4"}}{Very low-resolution template (2,562 vertices per hemisphere).}
#'     \item{\code{"fsaverage3"}}{Minimal-resolution template (642 vertices per hemisphere).}
#'     \item{\code{"fsmicro"}}{Micro test template (100 vertices per hemisphere).}
#'   }
#'
#' @return An integer giving the number of vertices per hemisphere for the
#'   specified template. Returns \code{NULL} invisibly if \code{fs_template}
#'   does not match any known template.
#'
#' @examples
#' count_vertices("fsaverage")   # 163842
#' count_vertices("fsaverage5")  # 10242
#' count_vertices("fsmicro")     # 100
#'
#' @export
#' 
count_vertices <- function(fs_template) {

  n_verts <- switch(fs_template,
                    fsmicro = 100,
                    fsaverage = 163842,
                    fsaverage6 = 40962,
                    fsaverage5 = 10242,
                    fsaverage4 = 2562,
                    fsaverage3 = 642)
  return(n_verts)
}


#' @title
#' Locate vertices belonging to specific cortical ROIs
#'
#' @description
#' This function identifies which vertices in a cortical surface mesh correspond
#' to one or more specified regions of interest (ROIs), based on the 36 regions of
#' the Desikan-Killiany atlas distributed with FreeSurfer (\code{aparc.annot}).
#' If no ROIs are specified, it returns a full lookup table of all ROIs with vertex
#' counts, proportions, and lobe assignments.
#'
#' @param rois Optional character vector of ROI names to locate (e.g.,
#'   \code{c("insula", "precentral")}). If \code{NULL} (default), the function
#'   returns a full ROI lookup table.
#' @param n_verts Integer. Total number of surface vertices to consider
#'   (default: \code{163842} "fsaverage").
#' @param hemi String: "lh" or "rh". Used to point to the hemisphere-specific
#'   aparc annotation file.
#' @param verbose Logical (default: \code{TRUE}).
#'
#' @return
#' \itemize{
#'   \item If \code{rois} is \code{NULL}, returns a \code{data.frame} with:
#'     \itemize{
#'       \item \code{roi_id}: ROI numeric ID
#'       \item \code{roi_label}: ROI name
#'       \item \code{vw_count}: Number of vertices in the ROI
#'       \item \code{vw_prop}: Proportion of total vertices in the ROI
#'       \item \code{roi_lobe}: Anatomical lobe classification
#'       \item \code{lobe_count}: Total number of vertices in the lobe
#'       \item \code{lobe_prop}: Proportion of vertices in the lobe
#'     }
#'   \item If \code{rois} is provided, returns a logical vector of length
#'     \code{n_verts}, where \code{TRUE} indicates the vertex belongs to one of
#'     the selected ROIs.
#' }
#'
#' @examples
#' # Get the full ROI lookup table
#' roi_table <- locate_roi()
#'
#' # Get a logical vector of vertices in the insula
#' insula_mask <- locate_roi(rois = "insula")
#'
#' @export
#'
locate_roi <- function(rois = NULL, n_verts = 163842, hemi = c('lh','rh'), verbose = TRUE) {
  # using Desikan-Killiany 36 regions

  hemi <- match.arg(hemi)
  annot <- aparc.annot[[hemi]]

  all_rois <- annot$label_codes[1:n_verts]
  roi_counts <- as.data.frame(table(all_rois))
  names(roi_counts) <- c('roi_id', 'vw_count')
  roi_counts['vw_prop'] <- round(roi_counts[,'vw_count']/n_verts, 3)

  notations <- annot$colortable_df[, c('struct_name', 'code', 'hex_color_string_rgb')]
  names(notations) <- c('roi_label', 'roi_id', 'roi_color')

  roi_lookup <- merge(notations, roi_counts, by = 'roi_id',
                      all.x = FALSE, all.y = TRUE)
  roi_lookup <- roi_lookup[order(-roi_lookup$vw_prop), ]

  # Group by lobe (vertex count noted next to each region)
  lobe_map <- list(
    frontal = c('superiorfrontal',      # 12.179
                'precentral',           # 10.740
                'rostralmiddlefrontal', #  7.243
                'lateralorbitofrontal', #  4.188
                'caudalmiddlefrontal',  #  3.736
                'paracentral',          #  3.294
                'parsopercularis',      #  3.119
                'medialorbitofrontal',  #  2.653
                'parstriangularis',     #  2.046
                'parsorbitalis',        #    956
                'frontalpole'),         #    272
    parietal = c('superiorparietal', # 10.456
                'postcentral',       #  9.519
                'supramarginal',     #  8.600
                'inferiorparietal',  #  7.871
                'precuneus'),        #  7.308
    temporal = c('superiortemporal',  # 7.271
                'fusiform',           # 4.714
                'middletemporal',     # 4.452
                'inferiortemporal',   # 4.415
                'bankssts',           # 2.137
                'parahippocampal',    # 1.838
                'entorhinal',         # 1.102
                'transversetemporal', # 1.064
                'temporalpole'),      #   839
    occipital = c('lateraloccipital', # 6.379
                'lingual',            # 4.205
                'pericalcarine',      # 1.912
                'cuneus'),            # 1.630
    cingulate = c('posteriorcingulate',      # 3.266  # (Parietal)
                'isthmuscingulate',          # 2.531  # (Parietal)
                'caudalanteriorcingulate',   # 1.439  # (Frontal)
                'rostralanteriorcingulate'), # 1.350. # (Frontal)
    insula = c("insula") # 5.229
  )

  # Named lookup vector: roi_label -> lobe
  roi_to_lobe <- stats::setNames(
    rep(names(lobe_map), lengths(lobe_map)), unlist(lobe_map, use.names = FALSE))

  roi_lookup$roi_lobe <- roi_to_lobe[roi_lookup$roi_label]
  roi_lookup$roi_lobe[is.na(roi_lookup$roi_lobe)] <- "none"

  lobe_stats <- stats::aggregate(
    cbind(lobe_count = vw_count, lobe_prop = vw_prop) ~ roi_lobe, data = roi_lookup, FUN = sum)

  # If no input ROIs are provided just return the entire lookup
  if (is.null(rois)) return(roi_lookup)

  # If only one ROI is provided make sure this is a vector
  rois <- as.vector(rois, mode = "character")

  selected_rois <- roi_lookup[roi_lookup[,'roi_label'] %in% rois, ]
  vw_message(
    c("i" = "Selected {.val {sum(selected_rois$vw_count)}} vertices ({sum(selected_rois$vw_prop)*100}%) 
      from {length(rois)} region{?s}: {.val2 {selected_rois$roi_label}}"), 
    verbose = verbose)

  roi_ids <- selected_rois[, 'roi_id']

  roi_locs <- all_rois %in% as.integer(roi_ids)

  return(roi_locs)
}

#' Format a Human-Readable Outcome Label
#'
#' Constructs a descriptive label for a FreeSurfer surface-based morphometric
#' outcome by combining hemisphere and measure names. Intended for use in
#' plot titles, model output tables, and report generation.
#'
#' @param hemi A character string specifying the hemisphere. Must be one of
#'   \code{"lh"} (left hemisphere) or \code{"rh"} (right hemisphere).
#'
#' @param measure A character string specifying the FreeSurfer surface measure.
#'   Must be one of:
#'   \describe{
#'     \item{\code{"thickness"}}{Cortical thickness.}
#'     \item{\code{"area"}}{Cortical surface area (white surface).}
#'     \item{\code{"area.pial"}}{Cortical surface area (pial surface).}
#'     \item{\code{"curv"}}{Mean curvature.}
#'     \item{\code{"jacobian_white"}}{Jacobian determinant of the white surface
#'       deformation field, reflecting distortion applied during spherical
#'       registration (\code{-surfreg}).}
#'     \item{\code{"pial"}}{Pial surface coordinates.}
#'     \item{\code{"pial_lgi"}}{Local gyrification index (pial surface).}
#'     \item{\code{"sulc"}}{Sulcal depth.}
#'     \item{\code{"volume"}}{Cortical gray matter volume.}
#'     \item{\code{"w_g.pct"}}{White/gray matter intensity ratio.}
#'     \item{\code{"white.H"}}{Mean curvature of the white surface.}
#'     \item{\code{"white.K"}}{Gaussian curvature of the white surface.}
#'   }
#'
#' @param sep A character string used to separate the hemisphere label from
#'   the measure name. Defaults to \code{" - "}.
#'
#' @return A single character string combining the hemisphere and measure
#'   labels, e.g. \code{"Left hemisphere - Cortical thickness"}.
#'   Returns a string with \code{NA} as the measure component if
#'   \code{measure} does not match any known value (due to \code{switch()}
#'   returning \code{NULL}, coerced by \code{paste0}).
#'
#' @details
#' Hemisphere labels map \code{"lh"} → \code{"Left"} and
#' \code{"rh"} → \code{"Right"}. Measure labels follow FreeSurfer's
#' naming conventions as produced by \code{mris_preproc} and
#' \code{mrisurf_thickness}. The \code{sep} argument allows flexible
#' formatting for different output contexts (e.g., \code{sep = ": "}
#' for axis labels, \code{sep = "\n"} for multi-line plot titles).
#'
#' @seealso \code{\link{count_vertices}} for template resolution lookup.
#'
#' @examples
#' outcome_name("lh", "thickness")
#' # "Left hemisphere - Cortical thickness"
#'
#' outcome_name("rh", "sulc", sep = ": ")
#' # "Right hemisphere: Sulcal depth"
#'
#' @export
#' 
outcome_name <- function(hemi, measure, sep =' - ') {

  hemi_name <- if (hemi == "lh") "Left" else "Right"

  meas_name <- switch(measure,
    thickness = 'Cortical thickness',
    area = 'Cortical surface area (white surface)',
    area.pial = 'Cortical surface area (pial surface)',
    curv = 'Mean curvature',
    jacobian_white = 'Jacobian determinant (white surface)', # how much the white surface was distorted in order to register to the spherical atlas during the -surfreg step
    pial = 'Pial surface coordinates',
    pial_lgi = 'Local gyrification index (pial surface)', # not in QDECR?
    sulc = 'Sulcal depth',
    volume = 'Cortical gray matter volume',
    w_g.pct = 'White/gray matter intensity ratio',
    white.H = 'Mean curvature (white surface)',
    white.K = 'Gaussian curvature (white surface)'
  )

  paste0(hemi_name, ' hemisphere', sep, meas_name)

}

probe_data_resolution <- function(subj_surf_dir,
                                  hemi,
                                  measure_file,
                                  fwhmc,
                                  fs_template) {

  mgh_file_pattern <- paste0("^", hemi, "\\.", measure_file,
                             ".*\\.fsaverage.*\\.mgh$")
  hits <- list.files(subj_surf_dir, pattern = mgh_file_pattern)

  if (length(hits) == 0L) {
    stop("No fsaverage*.mgh surface files found in: ", subj_surf_dir,
         "\nPlease check your measure, hemi, and fwhm settings.")
  }

  # Extract template (fsaverageX) and fwhm tags
  found_templates <- sub(".*\\.(fsaverage[^.]*)\\.mgh$", "\\1", hits)
  found_fwhm <- ifelse(
    grepl("fwhm[0-9]+", hits),
    sub(".*(fwhm[0-9]+).*", "\\1", hits),
    NA_character_
  )

  # Template selection ---------------------------------------------------------
  if (fs_template %in% found_templates) {
    tmpl <- fs_template
  } else {
    # Define resolution ranking (lower index = higher resolution)
    tmpl_order <- paste0('fsaverage', c('', 6:3))
    rank_tmpl <- function(x) match(x, tmpl_order)

    req_rank <- rank_tmpl(fs_template)

    found_ranks <- rank_tmpl(found_templates)
    non_na <- !is.na(found_ranks)

    # If requested is known and all found are lower-res (higher rank number)
    if (!is.na(req_rank) && any(non_na)) {
      min_found_rank <- min(found_ranks[non_na])

      if (req_rank < min_found_rank) {
        cli::cli_abort(c(
          "Could not find requested {.field {fs_template}} template in (some of) the files in:",
          " " = "{.path {dirname(dirname(subj_surf_dir))}}",
          " " = "Only lower-resolution templates {.pkg {unique(found_templates)}} found.",
          ">" = "Please either:",
          " " = "- Resample with {.code mri_surf2surf}",
          " " = "- Rerun FreeSurfer with {.field {fs_template}}",
          " " = "- Adjust {.arg fs_template} to one of the available templates"))
      }
    }

    # Otherwise, fall back to first available
    tmpl <- found_templates[1L]

    cli::cli_warn(c(
      "Data was not processed using requested {.field {fs_template}} template.",
      "i" = "Downsampling from {.field {tmpl}} instead. Note:",
      "!" = "This is fine for model tuning, but it introduces (small) registration
      errors. We therefore recommend using {.field fsaverage} in your final analysis.",
      "!" = "If you want cluster-wise pvalues for this template you will need to
      run {.code mri_glmfit --sim} to obtain the correct CDS files."))
    }

  # Restrict hits to chosen template
  idx_tmpl   <- found_templates == tmpl
  tmpl_hits  <- hits[idx_tmpl]
  tmpl_fwhm  <- found_fwhm[idx_tmpl]

  # FWHM selection ------------------------------------------------------------
  fh <- fwhmc
  available_fwhm <- unique(stats::na.omit(tmpl_fwhm))

  if (length(available_fwhm) == 0L) {
    if (!is.null(fwhmc) && nzchar(fwhmc)) {
      warning("No fwhm* tag found for template '", tmpl,
              "', ignoring requested fwhmc='", fwhmc, "'.")
    }
    fh <- ""
  } else {
    if (!is.null(fwhmc) && nzchar(fwhmc) && fwhmc %in% tmpl_fwhm) {
      fh <- fwhmc
    } else {
      fh <- available_fwhm[1L]
      if (!is.null(fwhmc) && nzchar(fwhmc) && fwhmc != fh) {
        warning("Requested fwhm '", fwhmc,
                "' not found. Using '", fh, "' instead.")
      }
    }
  }
  list(fs_template = tmpl, fwhmc = paste0(fh,'.'))
}
