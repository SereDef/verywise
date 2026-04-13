#' @title
#' Print pretty messages to console
#'
#' @description
#' Adds fillers to messages
#'
#' @param content : The message to print
#' @param fill : (default "=")
#' @param tot_nchar : (default: 60)
#' @param center : (default: TRUE) center the content.
#'
#' @return A prettyfyed console message.
#'
#' @author Serena Defina, 2025.
#'
prettify_message <- function(content, fill = "=", tot_nchar = 80, center = TRUE) {

  # Add space around the content
  content <- paste0(" ", content, " ")

  # Content length
  content_nchar <- nchar(content)

  # Repeat filling pattern
  rep_fill <- function(fill, n) {
    return (paste(rep(fill, n), collapse=""))
  }

  # Paste fillers and message
  if(content_nchar > tot_nchar) {
    tot_nchar <- content_nchar
    pretty_message <- paste0("\n", fill, content, fill, "\n")
  } else {
    fill_space <- tot_nchar - content_nchar
    if (center) {
      # If not even, add a space...?
      filler <- rep_fill(fill, ceiling(fill_space/2))
      pretty_message <- paste0("\n", filler, content, filler, "\n")
    } else {
      filler <- rep_fill(fill, fill_space)
      pretty_message <- paste0("\n", fill, content, filler, "\n")
    }
  }
  # TODO: make sure this always has the same nchar...?

  return(pretty_message)
}


#' @title
#' Print decorated message to console if verbose = TRUE
#'
#' @param msg Content to print
#' @param ... Other arguments to pass to prettify_message
#' @param verbose (default = TRUE)
#'
vw_pretty_message <- function(msg, ..., verbose = TRUE) {
  if (verbose) message(prettify_message(msg, ...))
  invisible()
}


vw_init_message <- function(model_type, verbose = TRUE) {
  cli::cli_rule(left=model_type, 
    right=paste0('verywise v', utils::packageVersion('verywise')))
}

# #' @title
# #' Print message to console if verbose = TRUE
# #'
# #' @param ... : content to print
# #' @param verbose (default = TRUE)
# #'
# vw_message <- function(..., verbose = TRUE) {
#   if (verbose) message(...)
#   invisible()
# }


VW_THEME <- list(
  ".pkg"  = list(color = "cyan", "font-weight" = "bold"),
  ".val"  = list(color = "blue", "font-weight" = "bold"),
  ".ok"   = list(color = "green"),
  ".warn" = list(color = "orange"),
  ".err"  = list(color = "red"),
  ".time" = list(color = "silver", "font-style" = "italic"),
  ".note" = list(color = "blue")
)

#' Print a styled message to the console
#'
#' Routes each call to the appropriate `cli` alert type based on the
#' conventional leading-character prefixes used in the existing codebase:
#'
#' | Prefix in `msg` | Rendered as |
#' |-----------------|-------------|
#' | `" * "` / `"* "` | bullet (`cli_bullets`) |
#' | `" ! "` / `"   !"` | warning alert |
#' | `"NOTE:"` | info alert |
#' | anything else | `cli_text` |
#'
#' Use the explicit typed variants (`vw_step`, `vw_bullet`, …) for new call
#' sites rather than relying on auto-detection.
#'
#' @param ...    Concatenated via `paste0()` to form the message string.
#' @param verbose Logical (default `TRUE`).
#' @param type  Optional override: `"info"`, `"step"`, `"bullet"`,
#'   `"warning"`, `"success"`, `"note"`, `"sub"`.
#' @export
#' 
vw_message <- function(..., verbose = TRUE, type = NULL) {
  if (!verbose) return(invisible(NULL))

  # caller environment for inline `{}` evaluation
  env <- parent.frame()

  # apply theme around a single message
  div_id <- cli::cli_div(theme = VW_THEME)
  on.exit(cli::cli_end(div_id), add = TRUE)

  msg <- paste0(...)
  msg <- gsub("^\n+", "", msg)
  msg <- gsub("\n+$", "", msg)

  if (!is.null(type)) {
    .vw_emit(type, msg, .envir = env)
    return(invisible(NULL))
  }

  if (grepl("^\\s{0,3}[*]\\s", msg)) {
    clean <- trimws(sub("^\\s*[*]\\s*", "", msg))
    .vw_emit("bullet", clean, .envir = env)
  } else if (grepl("^\\s{0,4}[!]\\s", msg)) {
    clean <- trimws(sub("^\\s*[!]\\s*", "", msg))
    .vw_emit("warning", clean, .envir = env)
  } else if (grepl("^NOTE:", msg)) {
    .vw_emit("note", msg, .envir = env)
  } else {
    .vw_emit("info", msg, .envir = env)
  }

  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Internal dispatcher
# ---------------------------------------------------------------------------

#' @keywords internal
.vw_emit <- function(type, msg, .envir) {
  switch(
    type,
    step    = {cli::cli_text("") # one blank line before the rule
               cli::cli_rule(msg)},  # rules don't use inline markup
    success = cli::cli_alert_success(msg, .envir = .envir),
    warning = cli::cli_alert_warning(msg, .envir = .envir),
    danger  = cli::cli_alert_danger(msg, .envir = .envir),
    note    = cli::cli_alert_info(msg, .envir = .envir),
    table   = cli_table(msg),
    bullet  = cli::cli_bullets(c("*" = msg), .envir = .envir),
    sub     = cli::cli_bullets(c(" " = msg), .envir = .envir),
    info    = cli::cli_text(msg, .envir = .envir)
  )
  invisible(NULL)
}

cli_table <- function(tab, title = NULL) {
  lev <- names(tab)
  cnt <- as.integer(tab)

  # column widths
  w1 <- max(nchar(c("Value", lev)))
  w2 <- max(nchar(c("N", as.character(cnt))))

  if (!is.null(title)) {
    cli::cli_h2(title)
    cli::cli_rule()
  }

  # header (bold)
  cli::cli_text(
    sprintf("{.strong %-*s}  {.strong %*s}", w1, "Value", w2, "N")
  )

  # rows
  for (i in seq_along(lev)) {
    cli::cli_text(
      sprintf("%-*s  %*d", w1, lev[i], w2, cnt[i])
    )
  }

  invisible(tab)
}

# VW_THEME <- list(
#   ".pkg"  = list(color = "cyan",    "font-weight" = "bold"),
#   ".step" = list(color = "silver"),
#   ".val"  = list(color = "yellow",  "font-weight" = "bold"),
#   ".ok"   = list(color = "green"),
#   ".warn" = list(color = "orange"),
#   ".err"  = list(color = "red"),
#   ".time" = list(color = "silver",  "font-style" = "italic"),
#   "h1"    = list(color = "cyan",    "font-weight" = "bold"),
#   "h2"    = list(color = "silver",  "font-weight" = "bold"),
#   ".note" = list(color = "blue")
# )

# vw_message <- function(..., verbose = TRUE, type = "info") {
#   if (!verbose) return(invisible(NULL))
#   switch(type,
#     info    = cli::cli_alert_info(...),
#     success = cli::cli_alert_success(...),
#     warning = cli::cli_alert_warning(...),
#     step    = cli::cli_h2(...)
#   )
# }

# vw_print_model <- function(formula, data, measure, hemi, verbose) {
#   if (!verbose) return(invisible(NULL))

#   # pretty-format formula (see format_formula() from prior fix)
#   f_str <- format_formula(formula)

#   cli::cli_div(theme = VW_THEME)
#   cli::cli_rule(left = "Model")
#   cli::cli_verbatim(f_str)
#   cli::cli_rule(left = "Data")
#   cli::cli_bullets(c(
#     "*" = "Measure   : {.val {measure}}",
#     "*" = "Hemisphere: {.val {hemi}}",
#     "*" = "Subjects  : {.val {length(unique(data$IDC))}}",
#     "*" = "Rows      : {.val {nrow(data)}}",
#     "*" = "Sessions  : {.val {paste(levels(data$ses), collapse = ' | ')}}"
#   ))
#   cli::cli_end()
# }
# ── Model ──────────────────────────────────────────────────────
# vw_area ~
#         dep_std * age_MRI +
#         assigned_sex +
#         edu_cat +
#         (1 | IDC)
# ── Data ───────────────────────────────────────────────────────
# • Measure   : thickness
# • Hemisphere: lh
# • Subjects  : 4 209
# • Rows      : 9 966
# • Sessions  : F09 | F13 | F17

# vw_print_sites <- function(site_list, verbose) {
#   if (!verbose) return(invisible(NULL))
#   cli::cli_div(theme = VW_THEME)
#   cli::cli_rule(left = "Sites")
#   cli::cli_bullets(setNames(
#     paste0("{.val ", names(site_list), "} — ",
#            sapply(site_list, nrow), " rows"),
#     rep("*", length(site_list))
#   ))
#   cli::cli_end()
# }
# ── Sites ───────────────────────────────────────────────────────
# • GenR   — 9 966 rows
# • ABCD   — 21 956 rows
# • ABCD2  — 5 201 rows

# vw_warn_summary <- function(warn_counts, verbose) {
#   if (!verbose || sum(warn_counts) == 0) return(invisible(NULL))

#   cli::cli_div(theme = VW_THEME)
#   cli::cli_rule(left = "Fit Warnings", col = "orange")

#   if (warn_counts[["singular"]] > 0)
#     cli::cli_alert_warning(
#       "{.warn {warn_counts[['singular']]}} singular fit(s) — consider simplifying random effects"
#     )
#   if (warn_counts[["convergence"]] > 0)
#     cli::cli_alert_warning(
#       "{.warn {warn_counts[['convergence']]}} convergence warning(s)"
#     )
#   if (warn_counts[["na_coef"]] > 0)
#     cli::cli_alert_warning(
#       "{.warn {warn_counts[['na_coef']]}} vertex(es) returned NA coefficients — skipped"
#     )

#   cli::cli_end()
# }
# # ── Fit Warnings ─────────────────────────────────── ⚠
# # ! 312 singular fit(s) — consider simplifying random effects
# # ! 4 vertex(es) returned NA coefficients — skipped

# vw_done <- function(output_path, elapsed, verbose) {
#   if (!verbose) return(invisible(NULL))

#   cli::cli_div(theme = VW_THEME)
#   cli::cli_rule()
#   cli::cli_alert_success(
#     "{.ok Done} in {.time {round(elapsed, 1)}}s — results saved to {.path {output_path}}"
#   )
#   cli::cli_end()
# }
# # ────────────────────────────────────────────────────────────────
# # ✓ Done in 87.3s — results saved to ./results/lh_thickness_dep/