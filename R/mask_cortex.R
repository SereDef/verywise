#' @title
#' Clean out vertices that are not on the cortex
#'
#' @description
#' Clean the super-subject matrix of the vertices that have 0 value or are not
#' located on the cortex, according to FreeSurfer \code{*cortical.mask.mgh} file
#' distributed in \code{inst/extdata}. This function is used in
#' \code{\link{build_supersubject}}.
#'
#' @param hemi : (default = "lh") hemisphere ("lh" or "rh")
#' @param fs_template : (default = "fsaverage")
#'
#' @return A (large) logical vector for cortical (or technical, i.e. vertex = 0) mask.
#'
mask_cortex <- function(hemi = "lh", fs_template = "fsaverage") {
  # Load mask file
  mask_file <- file.path(
    system.file("extdata", package = "verywise"), "cortex_mask",
    # paste0(hemi, ".", fs_template, ".cortex.mask.mgh")
    paste0(hemi, ".", fs_template, ".cortex.mask.rds")
  )
  mask <- readRDS(mask_file)
  # Turn into logical vector
  # mask <- as.logical(load.mgh(mask_file)$x) # TRUE where cortex
}


#'  @title
#' Save logical cortex mask from FreeSurfer cortex.label files
#'
#' @param freesurfer_home : path to freesurfer home
#' @param fsaverage_template : one of "fsaverage" (163842 vertices),
#' "fsaverage6" (40962), "fsaverage5" (10242), "fsaverage4" (2562),
#' "fsaverage3" (642)
#' @param hemi : "lh" of "rh"
#' @param outp_path : path where to save the logical mask
#'
create_cortex_mask <- function(freesurfer_home='',
                               fsaverage_template='fsaverage',
                               hemi='lh', outp_path=getwd()) {

  # Expecting FreeSurfer subject folder structure
  label_file <- file.path(freesurfer_home, 'subjects', fsaverage_template, 'label',
                          paste0(hemi, '.cortex.label'))

  tot_verts <- switch(fsaverage_template,
                      fsaverage = 163842,
                      fsaverage6 = 40962,
                      fsaverage5 = 10242,
                      fsaverage4 = 2562,
                      fsaverage3 = 642,
                      stop(sprintf("Unknown fsaverage template '%s'",
                                   fsaverage_template)))

  # The first line is a comment, and the 2nd one contains a single number: the number of vertex lines following.
  n_cortex_verts = utils::read.table(label_file, skip=1L, nrows=1L,
                                     col.names = c('num_verts'),
                                     colClasses = c("integer"))$num_verts[1]

  cortex_verts = utils::read.table(label_file, skip=2L,
                                   col.names = c('vertex_index', 'coord1', 'coord2', 'coord3', 'value'),
                                   colClasses = c("integer", "numeric", "numeric", "numeric", "numeric"))

  cortex_verts_idx = cortex_verts$vertex_index

  # Check input
  if(length(cortex_verts_idx) != n_cortex_verts) {
    stop(sprintf("Expected %d vertex rows in label file '%s' from header, but received %d.\n",
                 n_cortex_verts, label_file, length(cortex_verts_idx)))
  }
  if (n_cortex_verts > tot_verts) {
    stop(sprintf("%s should not contain more than %d vertices (%d found).\n",
                 fsaverage_template, tot_verts, n_cortex_verts))

  }

  # Create the binary mask
  cortex_mask <- logical(tot_verts)
  # NOTE: FreeSurfer 0-based index --> R 1-based indexing
  cortex_mask[cortex_verts_idx + 1L] <- TRUE

  print(stats::addmargins(table(cortex_mask)))

  saveRDS(cortex_mask, file.path(outp_path,
                                 paste0(hemi,'.', fsaverage_template,'.cortex.mask.rds')))

  NULL
}

