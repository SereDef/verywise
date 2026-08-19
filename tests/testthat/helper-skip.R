skip_if_no_freesurfer <- function() {
  fs_home <- Sys.getenv("FREESURFER_HOME")
  
  # Skip if environment variable is empty or the directory doesn't exist
  if (!nzchar(fs_home) || !dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not found in FREESURFER_HOME")
  }
}