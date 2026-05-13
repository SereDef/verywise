#' Print a plotting environment sit-rep
#'
#' Runs a full diagnostic of the Python/Kaleido/Xvfb environment used for
#' static brain map export. Useful for debugging \code{Error 525} and related
#' Kaleido/Chromium issues on HPC clusters.
#'
#' Mirrors the same Python setup as \code{\link{plot_vw_surf}}: calls
#' \code{py_require()} to ensure the uv environment is ready, then
#' \code{.vw_surf_init_py()} to load \code{patch_kaleido.py}, and finally
#' runs the Python sitrep. If Python cannot be initialised at all, a
#' partial R-side sitrep is printed instead.
#'
#' @return Invisible \code{NULL}. Output is printed to the console.
#'
#' @examples
#' \dontrun{
#' plot_sitrep()
#' }
#'
#' @export
plot_sitrep <- function() {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    vw_error(c(
      "Surface plotting requires {.pkg reticulate}.",
      "i" = "Install it with {.code install.packages('reticulate')}."
    ))
  }

  sitrep_file <- system.file("python", "plot_sitrep.py", package = "verywise")
  if (!nzchar(sitrep_file)) {
    stop(
      "verywise: could not find inst/python/plot_sitrep.py - ",
      "try reinstalling the package.",
      call. = FALSE
    )
  }

  # --- Step 1: initialise the Python environment -----------------------
  # Mirrors plot_vw_surf(): py_require() first (triggers uv env creation),
  # then .vw_surf_init_py() (loads patch_kaleido + renderer).
  # If either fails, fall back to the R-only partial sitrep.
  py_ready <- tryCatch({
    reticulate::py_require(c("nilearn", "numpy", "matplotlib",
                             "plotly", "kaleido", "choreographer", "logistro"))
    .vw_surf_init_py()
    TRUE
  }, error = function(e) FALSE)

  if (!py_ready) {
    message("\n[verywise] Could not initialise Python environment.")
    message("  Partial R-side checks follow - fix any issues shown,")
    message("  then re-run plot_sitrep() for the full report.\n")
    .r_sitrep()
    return(invisible(NULL))
  }

  # --- Step 2: run the full Python sitrep ------------------------------
  reticulate::py_run_string(paste0(
    "import importlib.util as _u\n",
    "_spec = _u.spec_from_file_location('plot_sitrep', r'", sitrep_file, "')\n",
    "_mod  = _u.module_from_spec(_spec)\n",
    "_spec.loader.exec_module(_mod)\n",
    "_mod.plot_sitrep()\n"
  ))

  invisible(NULL)
}


# --- R-only partial sitrep --------------------------
# Called only when Python cannot be initialised at all.

.r_sitrep <- function() {

  ok   <- function(msg) paste0("\u001b[32;1m\u2714  ", msg, "\u001b[0m")
  warn <- function(msg) paste0("\u001b[33;1m\u26a0  ", msg, "\u001b[0m")
  err  <- function(msg) paste0("\u001b[31;1m\u2718  ", msg, "\u001b[0m")
  info <- function(msg) paste0("\u001b[36m\u2139  ", msg, "\u001b[0m")
  row  <- function(label, status) message(sprintf("  %-30s %s", label, status))

  message("  \u001b[1mPlatform\u001b[0m")
  row("R version",  info(R.version$version.string))
  row("OS",         info(paste(Sys.info()[["sysname"]], Sys.info()[["release"]])))
  row("Node",       info(Sys.info()[["nodename"]]))
  row("reticulate", info(as.character(utils::packageVersion("reticulate"))))

  message("\n  \u001b[1mX Display (Xvfb)\u001b[0m")
  display <- Sys.getenv("DISPLAY")
  if (nzchar(display)) row("DISPLAY", ok(display)) else
    row("DISPLAY", err("not set - run: export DISPLAY=:99"))

  if (file.exists("/tmp/.X99-lock")) {
    row("Xvfb :99 lock file", ok("/tmp/.X99-lock exists"))
  } else {
    row("Xvfb :99 lock file",
        err("missing - run: Xvfb :99 -screen 0 1280x1024x24 &"))
  }

  dbus <- Sys.getenv("DBUS_SESSION_BUS_ADDRESS")
  if (nzchar(dbus)) row("DBUS_SESSION_BUS_ADDRESS", ok(substr(dbus, 1, 50))) else
    row("DBUS_SESSION_BUS_ADDRESS", warn("not set (Chromium may print D-Bus noise)"))

  message("\n  \u001b[1mChrome binary\u001b[0m")
  choreo <- path.expand(
    "~/.local/share/choreographer/deps/chrome-linux64/chrome")
  if (file.exists(choreo)) row("choreographer bundled", ok(choreo)) else
    row("choreographer bundled", warn("not found"))
  for (bin in c("chromium-browser", "chromium", "google-chrome")) {
    p <- Sys.which(bin)
    if (nzchar(p)) row(bin, ok(p)) else row(bin, warn("not found"))
  }

  message("\n  \u001b[1mKaleido env vars\u001b[0m")
  chrome_args <- Sys.getenv("KALEIDO_CHROME_ARGS")
  if (nzchar(chrome_args)) {
    row("KALEIDO_CHROME_ARGS", ok(chrome_args))
    display <- Sys.getenv("DISPLAY")
    if (nzchar(display) &&
        (grepl("--disable-gpu", chrome_args) ||
         grepl("--use-gl=swiftshader", chrome_args))) {
      row("Flag conflict",
          err("--disable-gpu or --use-gl=swiftshader + DISPLAY -> Error 525"))
    }
  } else {
    row("KALEIDO_CHROME_ARGS", warn("not set"))
  }

  message("")
}