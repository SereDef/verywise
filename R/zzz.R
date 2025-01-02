# Start-up message
.onAttach <- function(libname, pkgname) {
  version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
  packageStartupMessage("Welcome, ", pkgname," user!\nThis is version: ", version,
  "\nFor questions, issues, and bug reports, please see https://github.com/SereDef/verywise")
}


# Find FreeSurfer location on load
GLOBALS <- new.env()

.onLoad <- function(libname, pkgname) {

  FS_HOME <- Sys.getenv('FREESURFER_HOME')

  if (FS_HOME == "") {
    packageStartupMessage("ATTENTION! I could not find a FreeSurfer installation.",
    "This is currently necessary for computing significant clusters.")
  } else {
    assign('FS_HOME', FS_HOME, GLOBALS)
    packageStartupMessage(paste("FS_HOME:", FS_HOME))
  }

}
