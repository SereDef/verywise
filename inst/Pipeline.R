# library(verywise)
#
# # Simulate FreeSurfer and phenotype dataset
# subj_dir = "/Users/Serena/Desktop/Packages/verywise/inst/extdata/example_data"
#
# n_subs = 5 # * 2 cohorts * 2 timepoints * 2 hemispheres = 40 files
# n_vert = 100
#
# data_structure = list("cohort1" = list("sessions" = c("01",'02'), "n_subjects" = n_subs),
#                       "cohort2" = list("sessions" = c("01",'02'), "n_subjects" = n_subs))
#
# simulate_dataset(path = subj_dir,
#                  data_structure = data_structure,
#                  vw_resolution = n_vert,
#                  hemi = "lh")
# simulate_dataset(path = subj_dir,
#                  data_structure = data_structure,
#                  vw_resolution = n_vert,
#                  hemi = "rh")

# ==============================================================================
remove.packages("verywise")
.rs.restartR()
devtools::install()
library(verywise)
library(profvis)

# Check number of cores
# if (requireNamespace("parallelly", quietly = TRUE)){
#   if (n_cores > parallelly::availableCores()) {
#     warning("You only have access to ", parallelly::availableCores(), " cores. ",
#             parallelly::availableCores(omit = 1), " will be used.")
#     n_cores <- parallelly::availableCores(omit = 1)
#   } else if (n_cores >= parallelly::availableCores(omit = 1)) {
#     warning("You are using ", n_cores, " cores, but you only have ",
#             parallelly::availableCores(), " in total, other processes may get slower.")
#   }
# }
#

future::plan("multisession", workers = 8) # Should let the user do it instead..?
# future::plan(
#   list(
#     future::tweak(
#       future::multicore,
#       workers = 2),
#     future::tweak(
#       future::multisession,
#       workers = 2)
#   )
# ) # 8 cores in total

# pheno = read.csv(file.path(subj_dir, "phenotype.csv"))
#
# data_list = imp2list(pheno)
# m = length(data_list)
#
# out_stats <- lapply(data_list, single_lmm, y = rnorm(nrow(pheno), 6.5, 0.3),
#                     formula = vw_thickness ~ sex * age + site + (1|id),
#                     pvalues = (m == 1))
# verywise:::vw_pool(out_stats, 1)

# Run analysis
out <- run_vw_lmm(formula = vw_thickness ~ sex * age + site + (1|id),
                  subj_dir = subj_dir,
                  # pheno = pheno,
                  apply_cortical_mask = FALSE # 100 vertices, do not have mask
)

# ------------------------------------------------------------------------------
# margs <- c(quote(lme4::lmer()), # qdecr_decon(lme4::lmer()), # ??
#            list(formula = stats::as.formula(deparse(formula)),
#                 data = data_list[[1]]))
#
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
# set_arguments <- function(margs, model){
#     m <- names(margs)
#     if(is.null(m)) m <- rep("", length(margs))
#
#     f <- methods::formalArgs(get_function(model))
#     f2 <- formals(get_function(model))
#     # Rename arguments if necessary
#     b <- which(m[-1] == "")
#     m[b+1] <- f[!f %in% m[-1]][length(b)]
#     names(margs) <- m
#     # Add default arguments...?
#     margs <- c(margs, f2[!names(f2) %in% m])
#     if(is.symbol(margs$`...`)) margs$`...` <- NULL
#     margs
# }
#
# margs <- set_arguments(margs, model)
#
# # if (is.null(margs$formula)) stop("No `formula` set in margs.")
# # iii <- c("formula", "data", "method")
# # margs[iii] <- margs[!sapply(margs, is.symbol)][iii] # make sure these are not missing?
#
# mfz <- margs[match(c("formula", "data", "subset", "weights", "na.action", "offset"),
#                    names(margs),
#                    nomatch = 0L)]
# # if (!margs$method %in% 0:5) stop("The specified `method` for fastLm is not defined as a number between 0 and 5.")
# # Number of participants
# nr <- nrow(data_list[[1]])
#
# do.call2 <- function(what, args, ...) {
#   what <- get_function(what)
#   do.call(what, as.list(args), ...)
# }
#
# # Need to assign these before, <missing> does not work.. (?)
# mfz$weights = NULL
# mfz$offset = NULL
# mfz$na.action = na.fail
#
# mx <- lapply(data_list, function(x) {
#   mfz$data <- x
#   mfz$data[, paste0("qdecr_", measure)] <- 999
#   do.call2("stats::model.frame", c(mfz)) # list(na.action = na.fail)
# })
#
# nn <- length(data_list)
# mt <- attr(mx[[nn]], "terms") # why the last datasets?
#
# # Model weights
# w <- as.vector(stats::model.weights(mx[[nn]]))
# if (!is.null(w) && !is.numeric(w)) stop("'weights' must be a numeric vector")
#
# # TESTING ==========
# # What are you testing here? special terms with qdecr_ measure, how does this work?
# # mfz2 <- mfz
# # mfz2$data <- data_list[[nn]]
# # mfz2$data[, paste0("qdecr_", measure)] <- 888
# # mx_test <- do.call2("stats::model.frame", mfz2)
# # mx_test[, paste0("qdecr_", measure)] <- 999
# # if (!identical(mx[[nn]], mx_test)) stop ("Somewhere in your formula you specified a special term related to your vertex measure",
# #                                          " (interaction, polynomial, AsIs, etc); `qdecr_fastlm` currently does not support this.")
#
# y <- stats::model.response(mx[[nn]], "numeric")
# # if (nrow(mx[[nn]]) != nr) stop("The data that you are putting into the regression has missings! \n",
# #                                "QDECR can't handle that yet; we will fix this soon!")
# ys <- if(identical(unname(y), rep(999, nrow(mx[[nn]])))) "LHS" else "RHS"
# # Test if matrix is full rank (what does that mean??)
# if (ys == "LHS") {
#   mx_test2 <- mx_test
#   mx_test2b <- stats::model.matrix(mx_test2, object = mt)
#   if (Matrix::rankMatrix(mx_test2b) < ncol(mx_test2b)) stop ("The design matrix is NOT full rank. Please check if you have collinear columns in your data.")
# }
#
# # if (stats::is.empty.model(mt)) stop("The provided model (to fastLm) is empty. Check your data + formula.")
#
# mm <- NULL
#
# # ff <- "vw_fastlm_slow"
#
# if (paste0("qdecr_", measure) %in% colnames(attr(mt, "factors")) || ys == "LHS"){
#   # Create a design (or model) matrix
#   mm <- lapply(mx, stats::model.matrix, object = mt)
#
#   # ff <- "vw_fastlm"
#
#   vw <- list(mm = mm, # Design matrix
#              mf = mx[[1]], # Model frame object
#              # ff = ff, # Fast or slow implementation
#              formula = margs$formula, # formula
#              vertex = paste0("qdecr_", measure), # name of IV
#              y = y, # model.response object (with 999)
#              ys = ys, # "LHS" or "RHS"
#              w = w) # weights (NULL)
#              #method = margs$method, # NULL
#              #backing = prepvw$backing,
#              #backing_to_remove = prepvw$backing_to_remove,
#              #so = prepvw$so)
# } else {
#   warning("Your formula for `fastLm` contains complicated terms. \n",
#           "We will rely on the slower implementation of our fastLm.")
#   vw <- prepvw
# }
# #class(vw) <- c(vw$ff, "vw")
# #vw
