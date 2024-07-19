#' @title
#' List sub-directories till depth \code{n}
#'
#' @description
#' Small recursive utility to list directories given complex folder structure.
#' This is used in \code{\link{build_supersubject}} for listing all sub-directories
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

# ============================== MGH i/o helpers ==============================
#' @title
#' Load an MGH file into memory
#'
#' @param file_path : the (full) path to the mgh file.
#'
#' @author Heath Perdoe, 2013.
#' @return A mgh object with data and various header elements
#' @export
#'
load.mgh <- function(file_path) {
  to_read <- file(file_path, "rb")
  v <- readBin(to_read, integer(), endian = "big")
  ndim1 <- readBin(to_read, integer(), endian = "big")
  ndim2 <- readBin(to_read, integer(), endian = "big")
  ndim3 <- readBin(to_read, integer(), endian = "big")
  nframes <- readBin(to_read, integer(), endian = "big")
  type <- readBin(to_read, integer(), endian = "big")
  dof <- readBin(to_read, integer(), endian = "big")
  close(to_read)

  to_read <- file(file_path, "rb")
  dump <- readBin(to_read, double(), size = 4, n = 71, endian = "big")
  x <- readBin(to_read, double(), size = 4, n = ndim1 * nframes, endian = "big")
  close(to_read)

  list(
    x = x, v = v, ndim1 = ndim1, ndim2 = ndim2, ndim3 = ndim3, nframes =
      nframes, type = type, dof = dof
  )
}

#' @title
#' Create an object structured like an MGH object

#' @param x The vertex-wise values
#' @param version Version (default = 1)
#' @param ndim1 Width / 1st dimension (default = length(x))
#' @param ndim2 Height / 2nd dimension (default = 1)
#' @param ndim3 Depth / 3rd dimension (default = 1)
#' @param nframes Number of scalar components (default = 1)
#' @param type Data type, can be UCHAR (0), SHORT (4), INT (1) or FLOAT (3) (default = 3)
#' @param dof Degrees of freedom (default = 0)
#'
#' @author Sander Lamballais, 2018.
#' @return A mgh object with data and various header elements
#'
as.mgh <- function(x, version = 1L,
                   ndim1 = as.integer(length(x)),
                   ndim2 = 1L,
                   ndim3 = 1L, nframes = 1L, type = 3L, dof = 0L) {
  if (!is.vector(x)) {
    stop("as.mgh only support objects of class `vector` right now.")
  } else {
    out <- list(
      x = x, v = version,
      ndim1 = ndim1,
      ndim2 = ndim2,
      ndim3 = ndim3,
      nframes = nframes,
      type = type,
      dof = dof
    )
    return(out)
  }
}

#' @title
#' Save an MGH file from memory
#'
#' @description
#' R translation of save_mgh.m (by Heath Pardoe, 09/12/2013)
#'
#' @param vol The MGH object (as from load.mgh)
#' @param file_name File name to be used to save out the data
#'
#' @author Sander Lamballais, 2018.
#'
save.mgh <- function(vol, file_name) {
  MRI.UCHAR <-  0
  MRI.INT <-    1
  MRI.LONG <-   2
  MRI.FLOAT <-  3
  MRI.SHORT <-  4
  MRI.BITMAP <- 5
  MRI.TENSOR <- 6
  slices <- c(1:256)
  fid <- file(file_name, open = "wb", blocking = TRUE)
  on.exit(close(fid))
  width <- vol$ndim1
  height <- vol$ndim2
  depth <- vol$ndim3
  nframes <- vol$nframes
  writeBin(as.integer(1), fid, size = 4, endian = "big")
  writeBin(as.integer(width), fid, size = 4, endian = "big")
  writeBin(as.integer(height), fid, size = 4, endian = "big")
  writeBin(as.integer(depth), fid, size = 4, endian = "big")
  writeBin(as.integer(nframes), fid, size = 4, endian = "big")
  writeBin(as.integer(MRI.FLOAT), fid, size = 4, endian = "big")
  writeBin(as.integer(1), fid, size = 4, endian = "big")
  UNUSED.SPACE.SIZE <- 256
  USED.SPACE.SIZE <- (3 * 4 + 4 * 3 * 4)
  unused.space.size <- UNUSED.SPACE.SIZE - 2
  writeBin(as.integer(0), fid, size = 2, endian = "big")
  writeBin(as.integer(rep.int(0, unused.space.size)), fid, size = 1)
  bpv <- 4
  nelts <- width * height
  writeBin(vol$x, fid, size = 4, endian = "big")
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
