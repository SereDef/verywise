% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mask_cortex.R
\name{create_cortex_mask}
\alias{create_cortex_mask}
\title{@title
Save logical cortex mask from FreeSurfer cortex.label files}
\usage{
create_cortex_mask(
  freesurfer_home = "",
  fsaverage_template = "fsaverage",
  hemi = "lh",
  outp_path = getwd()
)
}
\arguments{
\item{freesurfer_home}{: path to freesurfer home}

\item{fsaverage_template}{: one of "fsaverage" (163842 vertices),
"fsaverage6" (40962), "fsaverage5" (10242), "fsaverage4" (2562),
"fsaverage3" (642)}

\item{hemi}{: "lh" of "rh"}

\item{outp_path}{: path where to save the logical mask}
}
\description{
@title
Save logical cortex mask from FreeSurfer cortex.label files
}
