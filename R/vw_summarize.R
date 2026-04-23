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

vw_summarize_model_fit <- function(fitstats, verbose = TRUE){

  if (!verbose) return(invisible(NULL))

  # Row 1: singular fits -------------------------------------------------------------
  singular_fits <- table(fitstats[1, ])

  total_ran <- sum(singular_fits)

  singular_count <- singular_fits["1"]
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