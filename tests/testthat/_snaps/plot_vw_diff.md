# vw_diff errors when maps have mismatched lengths

    Code
      plot_vw_diff(lh_a = lh_a_dummy, lh_b = c(1, 2, 3))
    Condition
      Error in `.hemi_diff()`:
      ! lh_a and lh_b have different lengths (10242 vs 3).
      i Both maps must be on the same surface template.

# vw_diff errors when one side is NULL but the other is provided

    Code
      plot_vw_diff(lh_a = lh_a_dummy, lh_b = NULL)
    Condition
      Error:
      ! Could not evaluate cli `{}` expression: `null_side`.
      Caused by error:
      ! object 'null_side' not found

# vw_diff errors with invalid object types

    Code
      plot_vw_diff(lh_a = data.frame(a = 1), lh_b = lh_b_dummy)
    Condition
      Error in `.resolve_map()`:
      ! `lh_a` must be a numeric vector or a file path, not <data.frame>.

# computes accurate difference maps from numeric vectors

    Code
      res <- plot_vw_diff(lh_a = lh_a_dummy, lh_b = lh_b_dummy, rh_a = rh_a_dummy,
        rh_b = rh_b_dummy, label_a = "Condition A", label_b = "Condition B",
        fs_template = "fsaverage5")
    Message
      i Difference map: Condition A - Condition B (20484 vertices)
             Min      Q1  Median    Mean      Q3     Max
             0.5     0.5    0.75    0.75       1       1
      * Condition A < Condition B in 0 vertices (0%)
      * Condition A > Condition B in 20484 vertices (100%)
      i Using local FreeSurfer mesh from /Applications/freesurfer/7.4.1
      ✔ Interactive brain map opened

