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

# vw_setup_cli_output <- function() {
#   if (!interactive()) {
#     old <- options(
#       cli.ansi = FALSE,
#       cli.unicode = FALSE,
#       cli.progress_clear = FALSE
#     )
#     return(old)
#   }
#   invisible(NULL)
# }

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
  ".val"  = list(color = "blue", "font-weight" = "bold"),
  ".val2" = list(color = "blue"),
  ".val3" = list(color = "cyan"),
  
  ".warn" = list(color = "orange"),
  ".err"  = list(color = "red"),
  ".time" = list(color = "silver", "font-style" = "italic")
  # '.file' = list() use default
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

  # Capture input — preserve named vectors, collapse plain strings
  msg <- c(...)

  # ── msg is a named vector: forward directly to cli_bullets / cli_warn ─────
  if (length(msg) > 1 || !is.null(names(msg))) {
    if (!is.null(type)) {
      switch(type,
        warning = cli::cli_warn(msg, .envir = env),
        danger  = cli::cli_abort(msg, .envir = env),
        note    = cli::cli_inform(msg, .envir = env),
        cli::cli_inform(msg, .envir = env)
      )
    } else {
      cli::cli_bullets(msg, .envir = env)
    }
    return(invisible(NULL))
  }

  # ── msg is a single string: ───────────────────────────────────────────────
  msg <- gsub("^\n+", "", msg)
  msg <- gsub("\n+$", "", msg)

  if (!is.null(type)) {
    .cli_format(type, msg, .envir = env)
    return(invisible(NULL))
  }

  # Quick patters for common message types
  if (grepl("^\\s{0,3}[*]\\s", msg)) {
    clean <- trimws(sub("^\\s*[*]\\s*", "", msg))
    .cli_format("bullet", clean, .envir = env)
  } else if (grepl("^\\s{0,4}[!]\\s", msg)) {
    clean <- trimws(sub("^\\s*[!]\\s*", "", msg))
    .cli_format("warning", clean, .envir = env)
  } else if (grepl("^NOTE:", msg)) {
    .cli_format("note", msg, .envir = env)
  } else {
    .cli_format("info", msg, .envir = env)
  }

  invisible(NULL)
}

#' @title Themed error
#'
#' @description
#' Wrapper around \code{cli::cli_abort()} that applies \code{VW_THEME} and
#' reports the call from the correct frame. Accepts the same \code{msg}
#' conventions as \code{vw_message()}: a plain string, a named character
#' vector (for bullet-style bodies), or inline cli markup.
#'
#' @param msg  A character string or named character vector passed to
#'   \code{cli::cli_abort()}. Use names \code{"i"}, \code{"x"}, \code{"v"},
#'   \code{">"} for bullet formatting, e.g.
#'   \code{c("Something failed.", "i" = "Check {.arg x}.")}.
#' @param ...  Additional arguments forwarded to \code{cli::cli_abort()}
#'   (e.g. \code{class}, \code{call}, \code{.internal}).
#' @param call Caller enviroment (`rlang::caller_env()`)
#'
#' @author Serena Defina, 2026.
#'
vw_error <- function(msg, ..., call = rlang::caller_env()) {
  div_id <- cli::cli_div(theme = VW_THEME)
  on.exit(cli::cli_end(div_id), add = TRUE)

  cli::cli_abort(msg, ..., call = call, .envir = rlang::caller_env())
}

# ---------------------------------------------------------------------------
# Internal dispatcher
# ---------------------------------------------------------------------------

#' @keywords internal
.cli_format <- function(type, msg, .envir) {
  switch(
    type,
    step    = {cli::cli_text("") # one blank line before the rule
               cli::cli_rule(msg)},  # rules don't use inline markup
    success = cli::cli_alert_success(msg, .envir = .envir),
    warning = cli::cli_alert_warning(msg, .envir = .envir),
    danger  = cli::cli_alert_danger(msg, .envir = .envir),
    note    = cli::cli_alert_info(msg, .envir = .envir),
    table   = .cli_table(msg),
    bullet  = cli::cli_bullets(c("*" = msg), .envir = .envir),
    sub     = cli::cli_bullets(c(" " = msg), .envir = .envir),
    info    = cli::cli_text(msg, .envir = .envir)
  )
  invisible(NULL)
}

#' @keywords internal
.cli_table <- function(tab, title = NULL) {
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

