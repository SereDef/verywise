# Start-up message
.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("Welcome, ", pkgname," user!\nThis is version: ", version,
  "\nFor questions, issues, and bug reports, please see https://github.com/SereDef/verywise")
}
