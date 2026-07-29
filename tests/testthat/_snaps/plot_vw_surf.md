# plot_vw_surf validates missing inputs

    Code
      plot_vw_surf()
    Condition
      Error in `plot_vw_surf()`:
      ! At least one of `lh` or `rh` must be supplied.

# plot_vw_surf validates output extension

    Code
      plot_vw_surf(lh = lh_coef_file, to_file = "plot.pdf")
    Condition
      Error in `plot_vw_surf()`:
      ! Output file name must end in '.png'.

# plot_vw_surf validates views and ROIs

    Code
      plot_vw_surf(lh = lh_coef_file, views = c("lateral", "invalid_view"))
    Condition
      Error in `plot_vw_surf()`:
      ! Invalid view: invalid_view
      i Please choose from lateral, dorsal, anterior, medial, ventral, or posterior

---

    Code
      plot_vw_surf(lh = lh_coef_file, roi_outline = "fake_roi", fs_home = subj_dir)
    Condition
      Error in `plot_vw_surf()`:
      ! Invalid ROI: fake_roi
      i Please choose from superiorfrontal, precentral, superiorparietal, postcentral, supramarginal, inferiorparietal, precuneus, superiortemporal, rostralmiddlefrontal, lateraloccipital, insula, fusiform, middletemporal, inferiortemporal, lateralorbitofrontal, lingual, caudalmiddlefrontal, posteriorcingulate, ..., frontalpole, or unknown

