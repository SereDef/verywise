skip_if_no_freesurfer <- function() {
  fs_home_bk <- "/Applications/freesurfer/7.4.1"  # mac only fallback

  fs_home <- if (dir.exists(fs_home_bk)) fs_home_bk else {
    Sys.getenv("FREESURFER_HOME")
  }

  # Skip if we ended up with an empty path or a non-existent directory
  if (!nzchar(fs_home) || !dir.exists(fs_home)) {
    testthat::skip("FreeSurfer not found in FREESURFER_HOME")
  }

  invisible(fs_home)
}