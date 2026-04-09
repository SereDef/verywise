# Start-up message
.onAttach <- function(libname, pkgname) {

  version <- utils::packageVersion(pkgname)

  # Be quiet in non-interactive / batch contexts (CRAN, R CMD check, SLURM)
  if (!interactive()) {
    packageStartupMessage(sprintf("Loaded %s %s", pkgname, version))
    return(invisible())
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    packageStartupMessage(
      sprintf("Welcome, %s user! This is version: %s", pkgname, version)
    )
    return(invisible())
  }

  welcome <- function() {
    cli::cli_inform(cli::cli({
      cli::cli_div(id = "welcome", 
                   theme = list(".pkg" = list(color = "cyan", `font-weight` = "bold"),
                                ".v" = list(color = "cyan"),
                                ".url" = list(color = "silver", `font-style` = "italic")))
      cli::cli_text("Welcome {.pkg verywise} user! This is version {.v {version}} {cli::symbol$smiley}")
      cli::cli_alert_info("For docs, questions and bug reports: {.url https://github.com/SereDef/verywise}")
      cli::cli_end()
    }),
    class = "packageStartupMessage")
  }

  welcome()

  # cli::cli_rule(
  #   left  = cli::style_bold(cli::col_cyan("verywise")),
  #   right = cli::col_silver(paste0("v", version))
  # )

  # Use FORK for bigstatsr on Unix to avoid PSOCK connection stacking
  if (.Platform$OS.type == "unix") {
    op <- getOption("bigstatsr.cluster.type")
    if (is.null(op) || op == "PSOCK") {
      options(bigstatsr.cluster.type = "FORK")
    }
  }

}

utils::globalVariables(c("roi_lobe", "vw_count", "vw_prop", "chunk", "i"))
