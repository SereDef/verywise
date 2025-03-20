#' @title
#' Load "phenotype" file into R based on its extension
#'
#' @description
#' Takes a path to a file and reads it into memory.
#'
#' @details
#' The function checks if the file exists, identifies the file extension, and
#' loads the data accordingly. Supported file types: .csv, .rds, .sav, and .txt.
#'
#' @param pheno_str A character string specifying the file path.
#' @param verbose (default = TRUE).
#'
#' @export
#' @return A data frame containing the loaded data
#'
load_pheno_file <- function(pheno_str, verbose=TRUE) {

  if (verbose) message('Loading phenotype file...')

  # Check that specified file exists
  if (!file.exists(pheno_str)) stop(paste0(pheno_str, " file not found."))

  # Extract file extension
  get_extension <- function(file) {
    file_name_vect <- strsplit(basename(file), split="\\.")[[1]] # split by "."
    file_extension <- file_name_vect[length(file_name_vect)] # pick last element
    return(file_extension)
  }
  file_extension <- get_extension(pheno_str)

  # Load the file based on its extension
  if (file_extension == "csv") {
    data <- utils::read.csv(pheno_str)
  } else if (file_extension == "rds") {
    data <- readRDS(pheno_str)
  } else if (file_extension == "sav") {
    if (!requireNamespace("haven", quietly = TRUE)) {
      stop("You need the 'haven' package to read .sav files.")
    }
    data <- haven::read_sav(pheno_str)
  } else if (file_extension == "txt") {
    data <- utils::read.table(pheno_str, header = TRUE, sep = "\t")
  } else {
    stop("Unsupported file extension. Please either provide one of: .rds, .csv, .sav, or .txt, or read the file in yourself.")
  }

  # Return the loaded data
  return(data)
}

#' @title
#' Check that object exists in the global environment
#'
#' @details This function checks if the specified object exists in the global environment.
#' It returns the object if found, otherwise raises an error with a message.
#'
#' @param obj_name A character string specifying the name of the object
#'
#' @return The object from the global environment (if it exists).
#'
check_pheno_obj <- function(obj_name) {

  # Check if the object exists in the global environment
  if (!exists(obj_name, envir = .GlobalEnv)) {
    stop(paste0("`", obj_name, "` not found in the global environment."))
  }

  # Return the object if it exists
  return(get(obj_name, envir = .GlobalEnv))
}
