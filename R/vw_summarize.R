.compute_median_range <- function(vec, digits) {
  q <- round(stats::quantile(vec, probs = c(0, 0.5, 1), na.rm = TRUE), digits)
  q <- stats::setNames(as.list(q), c('min', 'median', 'max'))
  q
}

.compute_cluster_stats <- function(ocn, coef, idx, result_path, digits = 4) {

  clust <- ocn[idx, ]
  betas <- coef[idx, ]

  n_clusters <- if (is.null(clust)) 0L else max(clust, na.rm = TRUE)

  out <- list(n_clusters = as.integer(n_clusters), 
              coef_summary = .compute_median_range(vec=betas, digits=digits))
  
  if (n_clusters > 0) {

    fs_summary_file <- list.files(path=dirname(result_path), 
      pattern = paste0(basename(result_path), '.*\\.stack', idx, '\\.cache.*\\.cluster.summary$'), 
      recursive = TRUE, full.names = TRUE)
    
    if (length(fs_summary_file) != 1L) {
      vw_message(c('!' = 'Expected exactly 1 cluster summary for stack {idx}',
                   'i' = '{length(fs_summary_file)} found: {fs_summary_file}',
                   '>' = 'Skipping summary.'))
    } else {
      fs_summary <- utils::read.table(fs_summary_file)
      my_summary <- tapply(betas, clust, summary)
      comb_summary <- cbind(fs_summary, 
        as.data.frame(do.call(rbind, my_summary[-1]))) # remove non-significant
      # fnames <- c('ClusterNo','Max','VtxMax','Size(mm^2)','MNIX','MNIY','MNIZ','CWP',
      #   'CWPLow','CWPHi','NVtxs','WghtVtx','Annot')
      summary_map <- c('cluster_id'='V1', 'cluster_size'='V11', 'cluster_area'='V4', 'cluster_p'='V8',
                       'peak_id'='V3', 'peak_ROI'='V13', 
                       'min_coef'='Min.', 'median_coef'='Median', 'max_coef'='Max.')
      out[['cluster_summary']] <- setNames(comb_summary[summary_map], names(summary_map))
    }
  }
  out
}

.print_median_range <- function(mat, idx=NULL, name='', 
  digits = 4, pad = 15, note = '', verbose = verbose) {
  
  vec <- if (is.null(idx)) mat else mat[idx, ]
  
  q <- .compute_median_range(vec=vec, digits=digits)
  
  n_space <- pad - cli::ansi_nchar(name, type = 'width')
  filler  <- strrep('\u00a0', max(n_space, 1))

  vw_message("* {.strong {name}}:{filler}{.val {cli_round(q[['median']], digits)}} 
  [{.val {cli_round(q[['min']], digits)}}, {.val {cli_round(q[['max']], digits)}}] 
  {.time {note}}")

  q
}

.print_cluster_stats <- function(ocn_mat, coef_mat, idx, name='', result_path='.', digits = 4, pad = 15, note = '',
     verbose = verbose) {
  
  clust_stats <- .compute_cluster_stats(ocn=ocn_mat, coef=coef_mat, idx=idx, result_path=result_path, digits=digits) 

  n_clusters <- clust_stats[['n_clusters']]
  q <- clust_stats[['coef_summary']]

  n_space <- pad - cli::ansi_nchar(name, type = 'width')
  filler  <- strrep('\u00a0', max(n_space, 1))
  filler2 <- strrep('\u00a0', 3 - cli::ansi_nchar(n_clusters, type = 'width'))

  msg <- '* {.strong {name}}:{filler}{n_clusters} clusters{filler2}|'

  if (n_clusters > 0) {
    msg <- paste(msg, "{.val {cli_round(q[['median']], digits)}} [{.val {cli_round(q[['min']], digits)}}, {.val {cli_round(q[['max']], digits)}}] {.time {note}}")
  }
  
  vw_message(msg, verbose = verbose)

  clust_stats
}

vw_summarize_model_fit <- function(fitstats, random_terms, verbose = TRUE){

  if (!verbose) return(invisible(NULL))

  # Row 1: singular fits -------------------------------------------------------------
  singular_fits <- table(fitstats[1, ], useNA = 'no')

  total_ran <- sum(singular_fits)

  singular_count <- if (is.na(singular_fits['1'])) 0 else singular_fits['1']
  singular_perc <- round(singular_count / total_ran * 100)

  vw_message('\nModel fit summary')
  vw_message('* {.strong Singular model fits}: {singular_count} ({.warn {singular_perc}}%)')

  aic <- .print_median_range(fitstats, 2, 'AIC', pad = 1, note = '* median [range]')
  cR2 <- .print_median_range(fitstats, 3, 'Conditional R\u00b2')
  mR2 <- .print_median_range(fitstats, 4, 'Marginal R\u00b2')
  
  icc <- lapply(seq_along(random_terms), function(i) {
    .print_median_range(fitstats, (4L + i), paste('ICC', random_terms[i]))
  })
  names(icc) <- random_terms

  invisible(
    list('Singular fits' = list(count = singular_count, percent = singular_perc),
        'AIC' = aic, 'Conditional R\u00b2' = cR2 , 'Marginal R\u00b2' = mR2, 
        'ICC' = icc)
  )
}

vw_summarize_model_est <- function(coef, term_names, verbose = TRUE) {

  term_name_length <- max(nchar(term_names)) + 1L

  out <- list()

  vw_message('\nModel estimates', verbose = verbose)
  for (n in seq_along(term_names)){
    note <- if (n == 1) '* median [range]' else ''
    term_name <- term_names[n]
    out[[term_name]] <- .print_median_range(mat=coef, idx=n, name=term_name, pad=term_name_length, 
      note = note, verbose = verbose)
  }

  invisible(out)
}

vw_summarize_model_clusters <- function(coef, clust, term_names, result_path, verbose = TRUE) {

  term_name_length <- max(nchar(term_names)) + 1L

  out <- list()

  vw_message('\nModel estimates', verbose = verbose)
  for (n in seq_along(term_names)){
    note <- if (n == 1) '* median coef [range]' else ''
    term_name <- term_names[n]
    out[[term_name]] <- .print_cluster_stats(ocn_mat=clust, coef_mat=coef, idx=n, name=term_name, result_path=result_path, 
                                             pad=term_name_length, note = note, verbose = verbose)
  }

  invisible(out)
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
    regexec('^(?:([^/]+)/)?([lr]h)\\.([^.]+)\\.', files, perl = TRUE))
  parsed <- Filter(function(x) length(x) == 4L, parsed)

  if (!length(parsed)) {
    vw_message('! No matching files found in: {.file {outp_dir}}')
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
  vw_message('Results directory: {.file {outp_dir}}')
  cli::cli_rule()

  for (subdir in names(result)) {
    vw_message('{.strong {subdir}}')
    measures <- result[[subdir]]

    # compute max measure name width for alignment
    max_w <- max(nchar(measures))

    for (measure in names(measures)) {
      hemis <- measures[[measure]]
      hemi_str <- paste(
        ifelse(hemis == 'lh', cli::col_blue('[lh]'), cli::col_red('[rh]')),
        collapse = '  '
      )
      n_space <- 12 - cli::ansi_nchar(measure, type = 'width')
      filler  <- strrep('\u00a0', max(n_space, 1))
      vw_message(c('*' = '{measure}{filler}{hemi_str}'))
    }
    cat('\n')
  }

  invisible(result)
}

#' @title Print difference map summary statistics
#' @keywords internal
.print_diff_stats <- function(all_diff, label_a, label_b, digits = 3, pad = 8) {

  qs <- stats::quantile(all_diff, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

  stat_vals <- round(
    c(Min = qs[[1]], Q1 = qs[[2]], Median = qs[[3]],
      Mean = mean(all_diff, na.rm = TRUE), Q3 = qs[[4]], Max = qs[[5]]), digits)
  
  stat_length <- max(pad, max(nchar(stat_vals)) + 2L)

  n_finite <- sum(is.finite(all_diff))

  # Pad to a fixed width with non-breaking spaces (regular spaces get
  # collapsed/trimmed by cli's bullet rendering -- see .print_median_range()
  # and the filler <- strrep('\u00a0', ...) convention used throughout this file)
  filler <- function(x, pad=stat_length) {
    n_space <- pad - cli::ansi_nchar(x, type = 'width')
    strrep('\u00a0', max(n_space, 1))
  }

  header <- paste(vapply(names(stat_vals), 
    function(nm) paste0(filler(nm), nm), character(1)), collapse = '')
  
  values <- paste(
    vapply(names(stat_vals), function(nm) {
      plain <- cli::cli_format(cli_round(stat_vals[[nm]], digits))
      paste0(filler(plain), "{.val {cli_round(stat_vals[['", nm, "']], digits)}}")
    }, character(1)),
    collapse = '')
  
  vw_message(c(
    'i' = 'Difference map: {.strong {label_a}} - {.strong {label_b}} ({.val {n_finite}} vertices)',
    ' ' = '{header}',
    ' ' = values
  ))

  invisible(list(stats = stat_vals, n_finite = n_finite))
}

#' @title Print one cutoff line (or pair of lines) for a difference map
#' @keywords internal
.print_diff_cutoff <- function(all_diff, cutoff, label_a, label_b, n_finite) {

  if (cutoff == 0) {
    n_neg <- sum(all_diff < 0, na.rm = TRUE)
    n_pos <- sum(all_diff > 0, na.rm = TRUE)
    p_neg <- round(100 * n_neg / n_finite, 1)
    p_pos <- round(100 * n_pos / n_finite, 1)

    vw_message(c(
      '*' = '{.strong {label_a} {"<"} {label_b}} in {.val {n_neg}} vertices ({.val {p_neg}}%)',
      '*' = '{.strong {label_a} {">"} {label_b}} in {.val {n_pos}} vertices ({.val {p_pos}}%)'
    ))
    return(invisible(NULL))
  }

  n_below <- sum(all_diff < -cutoff, na.rm = TRUE)
  n_above <- sum(all_diff > cutoff, na.rm = TRUE)

  p_below <- round(100 * n_below / n_finite, 1)
  p_above <- round(100 * n_above / n_finite, 1)

  vw_message(c(
    '*' = '{.strong {label_a} {"<<"} {label_b}} in {.val {n_below}} vertices ({.val {p_below}}%) [cut-off: {.val {-cutoff}}]',
    '*' = '{.strong {label_a} {">>"} {label_b}} in {.val {n_above}} vertices ({.val {p_above}}%) [cut-off: {.val {cutoff}}]'
  ))

  invisible(NULL)
}