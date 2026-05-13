  .print_median_range <- function(mat, idx, name, digits = 2, pad = 15, note = '',
     verbose = verbose) {

    summ <- summary(mat[idx, ], digits = digits)

    form <- paste0("%.",digits,"f")

    vmed <- sprintf(form, summ['Median'])
    vmin <- sprintf(form, summ['Min.'])
    vmax <- sprintf(form, summ['Max.'])

    n_space <- pad - cli::ansi_nchar(name, type = "width")
    filler  <- strrep("\u00a0", max(n_space, 1))

    vw_message("* {.strong {name}}:{filler}{vmed} [{vmin}, {vmax}] {.time {note}}")

    invisible(NULL)

  }

.print_cluster_stats <- function(ocn, coef, idx, name, digits = 2, pad = 15, note = '',
     verbose = verbose) {

  clust <- ocn[idx, ]
  n_clusters <- max(clust, na.rm = TRUE)

  n_space <- pad - cli::ansi_nchar(name, type = "width")
  filler  <- strrep("\u00a0", max(n_space, 1))

  msg <- "* {.strong {name}}:{filler}{n_clusters} clusters |"

  if (n_clusters > 0) {
    sign_coef <- coef[idx, which(clust > 0)]

    summ <- summary(sign_coef, digits = digits)

    form <- paste0("%.",digits,"f")

    vmed <- sprintf(form, summ['Median'])
    vmin <- sprintf(form, summ['Min.'])
    vmax <- sprintf(form, summ['Max.'])

    msg <- paste(msg, " {vmed} [{vmin}, {vmax}] {.time {note}}")
    
  }
  
  vw_message(msg)

  invisible(NULL)
}

vw_summarize_model_fit <- function(fitstats, verbose = TRUE){

  if (!verbose) return(invisible(NULL))

  # Row 1: singular fits -------------------------------------------------------------
  singular_fits <- table(fitstats[1, ], useNA = 'no')

  total_ran <- sum(singular_fits)

  singular_count <- if (is.na(singular_fits['1'])) 0 else singular_fits['1']
  singular_perc <- round(singular_count / total_ran * 100)

  vw_message('\nModel fit summary')
  vw_message("* {.strong Singular model fits}: {singular_count} ({.warn {singular_perc}}%)")

  .print_median_range(fitstats, 2, 'AIC', pad = 1, note = '* median [range]')
  .print_median_range(fitstats, 3, 'ICC')
  .print_median_range(fitstats, 4, 'Marginal R\u00b2')
  .print_median_range(fitstats, 5, 'Conditional R\u00b2')
  
  invisible(NULL)
}

vw_summarize_model_est <- function(coef, term_names, verbose = TRUE) {

  if (!verbose) return(invisible(NULL))

  term_name_length <- max(nchar(term_names)) + 1L

  vw_message('\nModel estimates')
  for (n in seq_along(term_names)){
    note <- if (n == 1) '* median [range]' else ''
    .print_median_range(coef, n, term_names[n], pad = term_name_length, note = note)
  }

   invisible(NULL)
}

vw_summarize_model_clusters <- function(coef, clust, term_names, verbose = TRUE) {

  if (!verbose) return(invisible(NULL))

  term_name_length <- max(nchar(term_names)) + 1L

  vw_message('\nModel estimates')
  for (n in seq_along(term_names)){

    note <- if (n == 1) '* median beta [range]' else ''
    .print_cluster_stats(clust, coef, n, term_names[n], pad = term_name_length, note = note)
  }

   invisible(NULL)
}

#' Summarise output directory: measures × hemispheres per subdirectory
#'
#' @param outp_dir Path to the top-level output directory
#' @return Invisibly returns a named list: one element per subdir, each
#'   containing a named list of measures → character vector of hemispheres
#' 
#' @export
#' 
vw_summarize_outp_dir <- function(outp_dir) {

  files <- list.files(outp_dir, recursive = TRUE)

  # Parse only files matching {subdir}/{hemi}.{measure}.* 
  parsed <- regmatches(files, 
    regexec("^(?:([^/]+)/)?([lr]h)\\.([^.]+)\\.", files, perl = TRUE))
  parsed <- Filter(function(x) length(x) == 4L, parsed)

  if (!length(parsed)) {
    vw_message("! No matching files found in: {.file {outp_dir}}")
    return(invisible(NULL))
  }

  df <- data.frame(
    subdir = vapply(parsed, function(x) if (nzchar(x[2])) x[2] else basename(outp_dir), character(1)),
    hemi = vapply(parsed, `[[`, character(1), 3L),
    measure = vapply(parsed, `[[`, character(1), 4L),
    stringsAsFactors = FALSE)

  # Build result: named list[subdir] -> named list[measure] -> hemi vector
  result <- lapply(
    split(df, df$subdir),
    function(d) lapply(
      split(d, d$measure),
      function(m) sort(unique(m$hemi))
    )
  )

  cli::cli_rule()
  vw_message("Results directory: {.file {outp_dir}}")
  cli::cli_rule()

  for (subdir in names(result)) {
    vw_message("{.strong {subdir}}")
    measures <- result[[subdir]]

    # compute max measure name width for alignment
    max_w <- max(nchar(measures))

    for (measure in names(measures)) {
      hemis <- measures[[measure]]
      hemi_str <- paste(
        ifelse(hemis == "lh", cli::col_blue("[lh]"), cli::col_red("[rh]")),
        collapse = "  "
      )
      n_space <- 12 - cli::ansi_nchar(measure, type = "width")
      filler  <- strrep("\u00a0", max(n_space, 1))
      vw_message(c("*" = "{measure}{filler}{hemi_str}"))
    }
    cat("\n")
  }

  invisible(result)
}