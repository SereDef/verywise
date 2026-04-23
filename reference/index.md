# Package index

## All functions

- [`as.mgh()`](https://seredef.github.io/verywise/reference/as.mgh.md) :
  Create an object structured like an MGH object

- [`barnard.rubin()`](https://seredef.github.io/verywise/reference/barnard.rubin.md)
  : Pooled degrees of freedom calculation

- [`build_supersubject()`](https://seredef.github.io/verywise/reference/build_supersubject.md)
  : Build "supersubject" by stacking all vertex data in one large
  file-backed matrix with dimensions n_subjects x n_vertices.

- [`check_pheno_obj()`](https://seredef.github.io/verywise/reference/check_pheno_obj.md)
  : Check that object exists in the global environment

- [`chunk_Ymats()`](https://seredef.github.io/verywise/reference/chunk_Ymats.md)
  : Compute the sum and sum of squares of Y

- [`compress_local()`](https://seredef.github.io/verywise/reference/compress_local.md)
  : Compress a site's results into a tar.gz archive

- [`compute_clusters()`](https://seredef.github.io/verywise/reference/compute_clusters.md)
  : Compute significant clusters of vertices

- [`convert_to_mgh()`](https://seredef.github.io/verywise/reference/convert_to_mgh.md)
  :

  Convert statistical result FBMs to FreeSurfer `.mgh` format

- [`create_cortex_mask()`](https://seredef.github.io/verywise/reference/create_cortex_mask.md)
  : @title Save logical cortex mask from FreeSurfer cortex.label files

- [`estimate_fwhm()`](https://seredef.github.io/verywise/reference/estimate_fwhm.md)
  : Estimate full-width half maximum (FWHM)

- [`fbm2mgh()`](https://seredef.github.io/verywise/reference/fbm2mgh.md)
  : Save file backed matrix (FBM) to MGH file

- [`fbm_col_has_0()`](https://seredef.github.io/verywise/reference/fbm_col_has_0.md)
  : Check if all row elements are 0 in FBM

- [`imp2list()`](https://seredef.github.io/verywise/reference/imp2list.md)
  : Convert imputation object to a list of dataframes

- [`init_progress_tracker()`](https://seredef.github.io/verywise/reference/init_progress_tracker.md)
  : Initialize per-chunk progress tracking

- [`list.dirs.till()`](https://seredef.github.io/verywise/reference/list.dirs.till.md)
  :

  List sub-directories till depth `n`

- [`load.mgh()`](https://seredef.github.io/verywise/reference/load.mgh.md)
  : Load an MGH file into memory

- [`load_pheno_file()`](https://seredef.github.io/verywise/reference/load_pheno_file.md)
  : Load "phenotype" file into R based on its extension

- [`locate_roi()`](https://seredef.github.io/verywise/reference/locate_roi.md)
  : Locate vertices belonging to specific cortical ROIs

- [`make_chunk_sequence()`](https://seredef.github.io/verywise/reference/make_chunk_sequence.md)
  : Define chunks of vertices for analyses

- [`mask_cortex()`](https://seredef.github.io/verywise/reference/mask_cortex.md)
  : Clean out vertices that are not on the cortex

- [`move_result_files()`](https://seredef.github.io/verywise/reference/move_result_files.md)
  : Move key result files from one directory to another (for sharing and
  visualization)

- [`plot_vw_map()`](https://seredef.github.io/verywise/reference/plot_vw_map.md)
  : Plot vertex-wise coefficient maps on a 3D cortical surface

- [`plot_vw_surf()`](https://seredef.github.io/verywise/reference/plot_vw_surf.md)
  : Plot vertex-wise maps on brain surfaces

- [`plot_vw_surf_r()`](https://seredef.github.io/verywise/reference/plot_vw_surf_r.md)
  : Plot vertex-wise maps on brain surfaces

- [`prettify_message()`](https://seredef.github.io/verywise/reference/prettify_message.md)
  : Print pretty messages to console

- [`run_voxw_lmm()`](https://seredef.github.io/verywise/reference/run_voxw_lmm.md)
  :

  Run voxel-wise linear mixed model using
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)

- [`run_vw_fed_aggr()`](https://seredef.github.io/verywise/reference/run_vw_fed_aggr.md)
  : Run federated aggregation for vertex-wise distributed LMM

- [`run_vw_fed_local()`](https://seredef.github.io/verywise/reference/run_vw_fed_local.md)
  : Run vertex-wise linear model on local site data

- [`run_vw_lmm()`](https://seredef.github.io/verywise/reference/run_vw_lmm.md)
  :

  Run vertex-wise linear mixed model using
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)

- [`run_vw_meta()`](https://seredef.github.io/verywise/reference/run_vw_meta.md)
  : Vertex-wise Random-Effects Meta-Analysis Across Studies

- [`save.mgh()`](https://seredef.github.io/verywise/reference/save.mgh.md)
  : Save an MGH file from memory

- [`significant_cluster_stats()`](https://seredef.github.io/verywise/reference/significant_cluster_stats.md)
  : Calculate Significant Cluster Statistics

- [`simulate_distrib_dataset()`](https://seredef.github.io/verywise/reference/simulate_distrib_dataset.md)
  : Simulate a distributed (multi-site) cross-sectional brain surface
  dataset

- [`simulate_freesurfer_data()`](https://seredef.github.io/verywise/reference/simulate_freesurfer_data.md)
  : Simulate longitudinal FreeSurfer vertex-wise data

- [`simulate_long_pheno_data()`](https://seredef.github.io/verywise/reference/simulate_long_pheno_data.md)
  : Simulate (longitudinal) phenotype data

- [`simulate_longit_dataset()`](https://seredef.github.io/verywise/reference/simulate_longit_dataset.md)
  : Simulate a longitudinal brain surface dataset with associated
  phenotype data

- [`single_lmm()`](https://seredef.github.io/verywise/reference/single_lmm.md)
  : Run a single linear mixed model and extract statistics

- [`subset_supersubject()`](https://seredef.github.io/verywise/reference/subset_supersubject.md)
  : Subset an existing supersubject matrix by matching folder IDs

- [`unpack_formula()`](https://seredef.github.io/verywise/reference/unpack_formula.md)
  : Unpack R formula

- [`update_progress_tracker()`](https://seredef.github.io/verywise/reference/update_progress_tracker.md)
  : Update within-chunk progress for a (milestone) vertex

- [`vw_error()`](https://seredef.github.io/verywise/reference/vw_error.md)
  : Themed error

- [`vw_message()`](https://seredef.github.io/verywise/reference/vw_message.md)
  : Print a styled message to the console

- [`vw_pool()`](https://seredef.github.io/verywise/reference/vw_pool.md)
  :

  Pool [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) model
  output across imputed datasets

- [`vw_pretty_message()`](https://seredef.github.io/verywise/reference/vw_pretty_message.md)
  : Print decorated message to console if verbose = TRUE

- [`with_parallel()`](https://seredef.github.io/verywise/reference/with_parallel.md)
  : Run a chunked foreach loop sequentially or in parallel
