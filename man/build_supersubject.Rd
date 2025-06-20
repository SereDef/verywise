% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_supersubject.R
\name{build_supersubject}
\alias{build_supersubject}
\title{Build "supersubject" by stacking all vertex data in one large file-backed
matrix with dimensions n_subjects x n_vertices.}
\usage{
build_supersubject(
  subj_dir,
  folder_ids,
  outp_dir,
  measure,
  hemi,
  n_cores,
  fwhmc = "fwhm10",
  fs_template = "fsaverage",
  backing,
  error_cutoff = 20,
  save_rds = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{subj_dir}{: path to the FreeSurfer data, this expects a verywise structure.}

\item{folder_ids}{: the vector of observations to include. This holds the relative
path (from `subj_dir`) to the FreeSurfer data folder (e.g. "site1/sub-1_ses-01").}

\item{outp_dir}{: output path, where logs, backing files and the matrix itself
(if `save_rds == TRUE`) will be stored.}

\item{measure}{: vertex-wise measure, used to identify files.}

\item{hemi}{: hemisphere, used to identify files.}

\item{n_cores}{: number of cores to use for parallel processing.}

\item{fwhmc}{: (default = "fwhm10") full-width half maximum value, used to identify files.}

\item{fs_template}{: (default = "fsaverage") template on which to register vertex-wise data.
The following values are accepted:
 * fsaverage (default) = 163842 vertices (highest resolution),
 * fsaverage6 = 40962 vertices,
 * fsaverage5 = 10242 vertices,
 * fsaverage4 = 2562 vertices,
 * fsaverage3 = 642 vertices
Note that, at the moment, these are only used to downsample the brain map, for faster
model tuning. `verywise` expects the input data to be always registered on the "fsaverage"
template and the final analyses should also be run using `fs_template = "fsaverage"`
to avoid (small) imprecisions in vertex registration and smoothing.}

\item{backing}{: (default = `outp_dir`) location to save the matrix \code{backingfile}.}

\item{error_cutoff}{: (default = 20) how many missing directories or brain surface files
for the function to stop with an error. If < `error_cutoff` directories/files are not
found a warning is thrown and missing files are registered in the `issues.log` file.}

\item{save_rds}{: (default = FALSE) save the supersubject file metadata for re-use
in other sessions.}

\item{verbose}{: (default = TRUE)}
}
\value{
A Filebacked Big Matrix with vertex data for all subjects (dimensions:
n_subjects x n_vertices)
}
\description{
Build "supersubject" by stacking all vertex data in one large file-backed
matrix with dimensions n_subjects x n_vertices.
}
\author{
Serena Defina, 2024.
}
