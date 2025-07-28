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
pretty_message <- function(content, fill = "=", tot_nchar = 60, center = TRUE) {

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
#' Print message to console if verbose = TRUE
#'
#' @param ... : content to print
#' @param verbose (default = TRUE)
#'
vw_message <- function(..., verbose = TRUE) {
  if (verbose) message(...)
  invisible()
}
