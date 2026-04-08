# Start-up message
.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("Welcome, ", pkgname," user!\nThis is version: ", version,
  "\nFor questions, issues, and bug reports, please see https://github.com/SereDef/verywise")

  # Use FORK for bigstatsr on Unix to avoid PSOCK connection stacking
  if (.Platform$OS.type == "unix") {
    op <- getOption("bigstatsr.cluster.type")
    if (is.null(op) || op == "PSOCK") {
      options(bigstatsr.cluster.type = "FORK")
    }
  }
}

utils::globalVariables(c("roi_lobe", "vw_count", "vw_prop", "chunk", "i"))
