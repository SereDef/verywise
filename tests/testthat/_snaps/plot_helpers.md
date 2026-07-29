# check_hemi errors on nonexistent files

    Code
      check_hemi("nonexistent_file.mgh", "fsaverage")
    Condition
      Error in `check_hemi()`:
      ! "nonexistent_file.mgh" file not found: 'nonexistent_file.mgh'

# check_hemi warns on mismatched vector lengths

    Code
      res <- check_hemi(long_vec, "fsaverage3")
    Message
      ! long_vec vector length (10000) does not match fsaverage3 template (642), I
        will try to subset it.

# check_mask validates logical inputs

    Code
      check_mask(c(1, 2, 3), "fsaverage")
    Condition
      Error in `check_mask()`:
      ! `c(1, 2, 3)` must be a logical vector or NULL.

# load_and_mask_coef warns and skips missing files

    Code
      res <- load_and_mask_coef("lh", "area", "stack999", res_dir, threshold = NULL)
    Message
      ! Coefficient file not found for lh, skipping:
        'fixtures/fs7/results/lh.area.stack999.coef.mgh'

# load_and_mask_coef applies CWS thresholds appropriately

    Code
      res <- load_and_mask_coef("lh", "area", stack_id, res_dir, threshold = "cws")

# resolve_mesh finds local FreeSurfer dir if it exists

    Code
      res <- resolve_mesh("fsaverage", mock_fs_home)
    Message
      i Using local FreeSurfer mesh from <FS_HOME>

# resolve_mesh returns NULL (Python download flag) if mesh is missing

    Code
      res <- resolve_mesh("fsaverage", mock_fs_home)
    Message
      fsaverage surface mesh not found in '$FREESURFER_HOME/subjects/' > falling back
      to nilearn.

