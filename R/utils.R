#' @title
#' List sub-directories till depth \code{n}
#'
#' @description
#' Small recursive utility to list directories given complex folder structure.
#' This is used in \code{\link{make_supersubject()}} for listing all sub-directories
#' until the correct level.
#'
#' @param path : the path where to begin search
#' @param n : the depth of the recursive process (i.e., how many levels of sub-directories to list)
#'
#' @author Serena Defina, 2024.
#' @return A list of directory/sub-directory names.
#'
list.dirs.till <- function(path, n) {
  res <- list.dirs(path, recursive = FALSE)
  if (n > 1) {
    add <- list.dirs.till(res, n - 1)
    return(add)
  } else {
    return(res)
  }
}

#' @title
#' Load an MGH file into memory
#'
#' @param input.file : the (full) path to the mgh file.
#'
#' @author Heath Perdoe, 2013.
#' @return A mgh object with data and various header elements
#' @export
#'
load.mgh <- function(input.file) {
  to.read <- file(input.file, "rb")
  v <- readBin(to.read, integer(), endian = "big")
  ndim1 <- readBin(to.read, integer(), endian = "big")
  ndim2 <- readBin(to.read, integer(), endian = "big")
  ndim3 <- readBin(to.read, integer(), endian = "big")
  nframes <- readBin(to.read, integer(), endian = "big")
  type <- readBin(to.read, integer(), endian = "big")
  dof <- readBin(to.read, integer(), endian = "big")
  close(to.read)
  to.read <- file(input.file,"rb")
  dump <- readBin(to.read, double(), size = 4, n = 71, endian = "big")
  x <- readBin(to.read ,double(), size = 4, n = ndim1*nframes, endian = "big")
  close(to.read)
  list(x = x, v = v, ndim1 = ndim1, ndim2 = ndim2, ndim3 = ndim3, nframes =
         nframes, type = type, dof = dof)
}


# ???? =========================================================================
# qdecr_decon <- function(x, y = environment()) {
#   y <- as.list(substitute(x, env = y))
#   y[-1] <- lapply(y[-1], eval)
#   y
# }
#
# # Get function from string =====================================================
# get_function <- function(x, ...) {
#   if(is.character(x)){
#     fn <- strsplit(x, "::")[[1]] # example > "stats" "lm"
#     x <- if(length(fn) == 1) {
#       get(fn[[1]], envir = parent.frame(), mode = "function", ...)
#     } else {
#       get(fn[[2]], envir = asNamespace(fn[[1]]), mode = "function", ...)
#     }
#   }
#   x
# }
#
# do.call2 <- function(what, args, ...) {
#   what <- get_function(what)
#   do.call(what, as.list(args), ...)
# }
#
#
# # ==============================================================================\
# # Set all arguments for given function call (model) and list of user defined
# # arguments (margs)
# set_arguments <- function(margs, model){
#   m <- names(margs)
#   if(is.null(m)) m <- rep("", length(margs))
#
#   f <- methods::formalArgs(get_function(model))
#   f2 <- formals(get_function(model))
#   # Rename arguments if necessary
#   b <- which(m[-1] == "")
#   m[b+1] <- f[!f %in% m[-1]][length(b)]
#   names(margs) <- m
#   # Add default arguments...?
#   margs <- c(margs, f2[!names(f2) %in% m])
#   if(is.symbol(margs$`...`)) margs$`...` <- NULL
#   margs
# }
