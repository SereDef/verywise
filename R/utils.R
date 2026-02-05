


# ============================= Folder navigators =============================
#' @title
#' List sub-directories till depth \code{n}
#'
#' @description
#' Small recursive utility to list directories given complex folder structure.
#' This is used in \code{\link{build_supersubject}} for listing all sub-directories
#' until the correct level.
#'
#' @param path : the path where to begin search
#' @param n : the depth of the recursive process (i.e., how many levels of sub-directories to list)
#'
#' @author Serena Defina, 2024.
#' @return A list of directory/sub-directory names.
#'
list.dirs.till <- function(path, n) {
  res <- list.dirs(path, recursive = FALSE)

  if (n > 1) {
    add <- list.dirs.till(res, n - 1)
    return(add)
  } else {
    return(res)
  }
}

# ======================== Simulation helpers ==================================

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

  # Groups
  roi_lookup <- roi_lookup %>%
    dplyr::mutate(roi_lobe = dplyr::case_when(
      roi_label %in% c('superiorfrontal',      # 12.179
                       'precentral',           # 10.740
                       'rostralmiddlefrontal', #  7.243
                       'lateralorbitofrontal', #  4.188
                       'caudalmiddlefrontal',  #  3.736
                       'paracentral',          #  3.294
                       'parsopercularis',      #  3.119
                       'medialorbitofrontal',  #  2.653
                       'parstriangularis',     #  2.046
                       'parsorbitalis',        #    956
                       'frontalpole'           #    272
                       ) ~ "frontal",
      roi_label %in% c('superiorparietal', # 10.456
                       'postcentral',      #  9.519
                       'supramarginal',    #  8.600
                       'inferiorparietal', #  7.871
                       'precuneus'         #  7.308
                       ) ~ "parietal",
      roi_label %in% c('superiortemporal',   # 7.271
                       'fusiform',           # 4.714
                       'middletemporal',     # 4.452
                       'inferiortemporal',   # 4.415
                       'bankssts',           # 2.137
                       'parahippocampal',    # 1.838
                       'entorhinal',         # 1.102
                       'transversetemporal', # 1.064
                       'temporalpole'        #   839
                       ) ~ "temporal",
      roi_label %in% c('lateraloccipital', # 6.379
                       'lingual',          # 4.205
                       'pericalcarine',    # 1.912
                       'cuneus'            # 1.630
                       ) ~ "occipital",
      roi_label %in% c('posteriorcingulate',       # 3.266  # (Parietal)
                       'isthmuscingulate',         # 2.531  # (Parietal)
                       'caudalanteriorcingulate',  # 1.439  # (Frontal)
                       'rostralanteriorcingulate'  # 1.350. # (Frontal)
                       ) ~ "cingulate",
      roi_label %in% c('insula' # 5.229
                       ) ~ "insula",

      TRUE ~ "none"
    )) %>%
    dplyr::group_by(roi_lobe) %>%
    dplyr::mutate(lobe_count = sum(.data$vw_count),
                  lobe_prop = sum(.data$vw_prop)) %>%
    as.data.frame()

  # If no input ROIs are provided just return the entire lookup
  if (is.null(rois)) {
    return(roi_lookup)
  }

  # If only one ROI is provided make sure this is a vector
  rois <- as.vector(rois, mode = "character")

  selected_rois <- roi_lookup[roi_lookup[,'roi_label'] %in% rois, ]
  vw_message(' * selected ', sum(selected_rois$vw_count), ' vertices (',
             sum(selected_rois$vw_prop)*100, '%)', verbose = verbose)

  roi_ids <- selected_rois[, 'roi_id']

  roi_locs <- all_rois %in% as.integer(roi_ids)

  return(roi_locs)
}

# ==============================================================================

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


# ==============================================================================
# Read and save annotation files for internal use 
# fs_home = "/Applications/freesurfer/7.4.1"
# hemis <- c("lh", "rh")
# aparc.annot <- lapply(hemis, function(hemi){
#   annot_file <- file.path(fs_home, "subjects/fsaverage/label", 
#                           paste0(hemi, ".aparc.annot"))
#   freesurferformats::read.fs.annot(annot_file)})
# names(aparc.annot) <- hemis
# usethis::use_data(aparc.annot, internal = TRUE, overwrite = TRUE)
# ==============================================================================