#' @title
#' Save file backed matrix (FBM) to MGH file
#'
#' @author Sander Lamballais, 2018.
#' Modified version of save_mgh.m (by Heath Pardoe, 09/12/2013)
#'
#' @param fbm a file backed matrix
#' @param fname file name to be used to save out the data
#' @param filter a vector of indices of rows to include of fbm
#'
#' @export
#'
fbm2mgh <-function(fbm, fname, filter = NULL) {

    MRI.UCHAR <-  0
    MRI.INT <-    1
    MRI.LONG <-   2
    MRI.FLOAT <-  3
    MRI.SHORT <-  4
    MRI.BITMAP <- 5
    MRI.TENSOR <- 6
    slices <- c(1:256)

    # Create .mgh file
    fid <- file(fname, open = "wb", blocking = TRUE)
    on.exit(close(fid))

    # Define dimensions
    if (!is.null(filter)){
      if (max(filter) > fbm$nrow) stop("In fbm2mgh, the defined filter exceeds the bounds.")
      if (min(filter) < 1) stop("In fbm2mgh, the defined filter contains non-positive numbers.")
    } else {
      filter <- seq_len(fbm$nrow)
    }

    width <- fbm$ncol
    height <- 1
    depth <- 1
    nframes <- length(filter)

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

    for (i in filter) writeBin(fbm[i, ], fid, size = 4, endian = "big")

    NULL
  }
