#' @title
#' Setup a Future Plan for Parallel Processing
#'
#' @description
#' This function sets up a parallel processing plan using the `future` package.
#' It allows the user to specify the number of cores to use and the type of future
#' plan (for more information see the futureverse documentation) `multicore`)
#'
#' @param n_cores Integer. The number of cores to be used. Must be less than or equal to the available cores on the system.
#' @param plan_type Character. The type of future plan. Options are `"multisession"`, `"multicore"`, `"sequential"`, or `"cluster"`.
#' @return Invisible NULL. The function sets the future plan for parallel processing.
#' @importFrom future availableCores plan
#' @examples
#' \dontrun{
#'   setup_future_plan(cores = 4, plan_type = "multisession")
#' }
#' @export
setup_parallel <- function(n_cores, plan_type) {
  # Load the future package
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("The 'future' package is required but not installed. Install it with install.packages('future').")
  }

  # Check that cores is a positive integer
  if (!is.numeric(cores) || cores <= 0 || cores != as.integer(cores)) {
    stop("'cores' must be a positive integer.")
  }

  # Get the number of available cores
  max_cores <- future::availableCores()
  if (cores > max_cores) {
    stop(sprintf("Requested %d cores, but only %d are available.", cores, max_cores))
  }

  # Validate the plan_type
  valid_plans <- c("sequential", "multisession", "multicore", "cluster")
  if (!plan_type %in% valid_plans) {
    stop(sprintf(
      "'plan_type' must be one of: %s. You provided '%s'.",
      paste(valid_plans, collapse = ", "),
      plan_type
    ))
  }

  # On non-Unix systems, restrict the use of "multicore"
  if (plan_type == "multicore" && .Platform$OS.type != "unix") {
    stop("'multicore' is only supported on Unix-based systems. Use 'multisession' instead.")
  }

  # Set the future plan
  future::plan(strategy = plan_type, workers = cores)

  message(sprintf("Future plan set to '%s' with %d workers.", plan_type, cores))
  invisible(NULL)
}
