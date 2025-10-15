plot_vw_results <- function(res_dir, term, hemi = c("lh", "rh"), measure = 'area', outline_rois, fs_home) {

  # Hemisphere 
  hemi <- match.arg(hemi)
  hemi_name <- if (hemi == "lh") "left" else "right"

  # Retrieve fsaverage surface meshes
  surf_dir <- file.path(fs_home, "subjects/fsaverage/surf")
  surf_mesh <- list(
    pial = file.path(surf_dir, paste0(hemi, ".pial")),
    inflated = file.path(surf_dir, paste0(hemi, ".inflated"))
  )

  # Find term stack
  stack_file <- file.path(res_dir, "stack_names.txt")
  if (file.exists(stack_file)) {
    # Read existing file
    stack_ids <- utils::read.table(stack_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # Extract stack
    stack <- paste0('stack', stack_ids[ stack_ids$stack_name == term, 'stack_number'])

  } else {
    stop(sprintf("Cannot find the `stack_names.txt` in `%s`.", res_dir),
           " Results folder is incorrect or corrupted.")
  }

  # Set up python env 
  reticulate::py_require(c("nibabel", "nilearn", "matplotlib", "plotly"))
  # py_config()
  # ensure the packages go into the Python used by reticulate.
  # py_install(c("nibabel", "nilearn", "matplotlib"), pip = TRUE)


  # OR create / point to a conda environment
  # conda_create("wv_env")
  # conda_install("wv_env", c("nibabel", "nilearn"))
  # use_condaenv("wv_env", required = TRUE)
  
  # import Python modules
  nb <- reticulate::import("nibabel")
  nl <- reticulate::import("nilearn.plotting")
  np <- reticulate::import("numpy")

  # Load coefficients .MGH file
  coef_map_file <- file.path(res_dir, paste(hemi, measure, stack, "coef.mgh", sep="."))
  coef_map <- nb$load(coef_map_file)$get_fdata()

  # Load cluster info
  ocn_map_file <- file.path(res_dir, paste(hemi, measure, stack, "cache.th30.abs.sig.ocn.mgh", sep = "."))
  ocn_map <- nb$load(ocn_map_file)$get_fdata()

  clusters <- table(ocn_map)
  print(clusters)

  # roi_lookup <- verywise::locate_roi()
  roi_mask <- lapply(seq_along(outline_rois), function(roi_n) {
    roi_name <- outline_rois[roi_n]
    mask_n <- ifelse(verywise::locate_roi(roi_name, verbose = FALSE), roi_n, 0)
    return(mask_n)
  })
  roi_mask <- Reduce(`+`, roi_mask)

  surf_fig <- nl$plot_surf_stat_map(
    surf_mesh = surf_mesh$inflated,
    stat_map = coef_map,
    cmap = "coolwarm",
    colorbar = TRUE,
    hemi = hemi_name,
    engine = "plotly",
    vmin = min(coef_map), 
    vmax = max(coef_map),
    threshold = 0.0001,
    cbar_tick_format = ".3f",
    view = 'medial',
    title = paste('Effect of', term, 'on', measure, '-', hemi_name, 'hemisphere')
  )

  surf_fig$add_contours(
    roi_map = np$array(roi_mask),
    levels = as.list(seq_along(outline_rois)),
    labels = as.list(outline_rois),
    lines = list(list(width = 5))
  )

  # Launch browser 
  surf_fig$show() 

  return(surf_fig)
}

# library(freesurferformats)
# coef <- read.fs.mgh(file.path(fs7_data_path, "results", "lh.area.stack3.coef.mgh"))
# surf <- read.fs.surface(file.path(fs_home, "subjects/fsaverage/surf/lh.inflated"))


# Simple HTML plot (no ROI contours)
# surf_plot <- nl$view_surf(
#   surf_mesh = surf_mesh$inflated,
#   surf_map = coef_map,
#   cmap = "coolwarm",   # diverging colormap
#   colorbar = TRUE,
#   # threshold = 0.01     # optional: only show |values| > 0.01
# )

for (term in c("wisdom", "age", "sexMale", "(Intercept)")) {

  p = plot_vw_results(res_dir = "tests/testthat/fixtures/fs7/results",
                  term = term,
                  hemi = "lh",
                  measure = "area",
                  outline_rois = c('temporalpole', 'frontalpole','entorhinal'), 
                  fs_home = "/Applications/freesurfer/7.4.1")
}
