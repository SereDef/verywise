# errors if res_dir does not exist

    Code
      plot_vw_map("/nonexistent/path", term = "age")
    Condition
      Error in `plot_vw_map()`:
      ! Results directory does not exist: '/nonexistent/path'

# errors if term is absent from stack_names.txt

    Code
      plot_vw_map(res_dir, term = "__no_such_term__")
    Condition
      Error in `plot_vw_map()`:
      ! Term "__no_such_term__" not found in 'stack_names.txt'.
      i Available terms: "(Intercept)", "sexMale", "age", or "wisdom"

# errors when coef file is missing for requested hemisphere

    Code
      plot_vw_map(dir, term = "age", hemi = "lh", threshold = NULL)
    Message
      ! Coefficient file not found for lh, skipping:
        '<TEMP_DIR>/lh.area.stack3.coef.mgh'
    Condition
      Error in `plot_vw_map()`:
      ! No coefficient MGH files found for term "age" / measure "area".
      i Expected e.g. 'lh.area.stack3.coef.mgh' in '<TEMP_DIR>'

# warns and continues when stack_names.txt is missing (meta-analysis mode)

    Code
      res <- plot_vw_map(dir, term = "age", hemi = "lh", threshold = NULL)
    Message
      i Cannot find 'stack_names.txt' in
        '<TEMP_DIR>'.
        Assuming this is a meta-analysis.
      i Using local FreeSurfer mesh from /Applications/freesurfer/7.4.1
      ✔ Interactive brain map opened

