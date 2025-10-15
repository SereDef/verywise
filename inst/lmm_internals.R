fast_llm <- function(formula, data = NULL, REML = TRUE,
                     control = lme4::lmerControl(), start = NULL, verbose = 0L, 
                     subset, weights, na.action, offset, contrasts = NULL){
  
  # match.call() captures the call the user typed (a language object).
  # mc will be mutated and evaluated internally;
  # mcout preservs the original call to return at the end. # TODO: remove not needed
  # TODO: removed validate control, assume it is a lmerControl object for now
  mc <- mcout <- match.call()

  # Parse and pre-compute model matrices without fitting ---------------------------
  # Replace the function name in the captured call with lme4::lFormula
  # This parses the formula + data and produces a structure describing the model: 
  # (model frame, fixed and random design matrices, random-effects structure, etc.)
  mc[[1]] <- quote(lme4::lFormula)
  # Evaluate it in the caller’s environment
  lmod <- eval(mc, parent.frame(1L))

  # Store the original formula in mcout (for output) and remove it from lmod (to avoid redundancy).
  mcout$formula <- lmod$formula
  lmod$formula <- NULL

  # Extract the response from the model frame (y) and error if it’s a multi-column response 
  # (i.e., multivariate response) # TODO removed this check, not needed probably
  
  # Create deviance function for covariance parameters (theta)
  # builds a function devfun(theta) that computes the (RE)ML or ML deviance (negative log-likelihood) 
  # for a given covariance-parameter vector theta.
  # returns a closure whose environment contains precomputed matrices and data needed for fast evaluations.
  # Separates model construction from optimization:
  # The deviance function is reusable (profiling, optimization, derivative checks).
  # By building a closure that stores sparse matrices, it avoids rebuilding heavy objects repeatedly 
  # during optimization — big speed gain.
  # devfun is a function (callable in R). Its environment(devfun) contains a pp object 
  # (preprocessing/prediction helper), lower bounds, X, y, and other precomputed pieces.
  # pp is an object (list/class) with:
  # - X  : fixed-effects design matrix
  # - Zt : random-effects transpose sparse matrix
  # - y  : response vector
  # - n  : number of obs
  # - lower : numeric vector of lower bounds for theta (usually zeros)
  # - getPars(), setPars() methods / functions to expand theta -> full parameter vectors
  devfun <- do.call(mkLmerDevfun,
                    c(lmod, list(start=start, verbose=verbose, control=control)))
    
  # environment(devfun) contains:
  # $pp     # preprocessor object: X, Zt, y, methods to build/expand parameters
  # $lower  # numeric vector (lower bounds for theta)
  # plus helper functions used inside the deviance function
  
  # Optimize deviance function over covariance parameters

  # Set calc.derivs to either the user-supplied control$calc.derivs or a default boolean: 
  # compute derivative checks only if n is smaller than a threshold check.conv.nobsmax.
  # The %||% operator is a convenience: return left if not NULL, otherwise return right.
  # Derivative calculations can be expensive for large datasets; they are enabled by default 
  # only for smaller problems.
  calc.derivs <- control$calc.derivs %||% (nrow(lmod$fr) < control$checkConv$check.conv.nobsmax)

  # optimizeLmer(devfun, ...) which runs an optimizer (e.g., nlopt, bobyqa, optim, etc.) 
  # to minimize devfun. optimizeLmer handles restarts at parameter-space edges, boundary 
  # tolerances, and optimization control.
  # Provides robust default optimization behavior with various fallbacks.
  # opt structure:
  # opt$par      # numeric vector of optimized covariance params ("theta")
  # opt$fval     # numeric — final deviance value
  # opt$converged or opt$conv  # convergence indicator
  # attr(opt, "derivs") # derivative info (gradient/hessian) attached as attribute
  opt <- optimizeLmer(devfun, 
                      optimizer = control$optimizer,
                      restart_edge = control$restart_edge,
                      boundary.tol = control$boundary.tol,
                      control = control$optCtrl,
                      verbose=verbose,
                      start=start,
                      calc.derivs=calc.derivs,
                      use.last.params=control$use.last.params)
  
  # Convergence checks
  # examine derivative information and the final parameters for possible convergence problems
  # (singular fits, boundary solutions, failed gradients).
  # To detect common problems (e.g., optimizer stopped at the boundary, gradient too large, 
  # or other convergence issues) and to return diagnostics to the user.
  # It uses:
  # attr(opt, "derivs") — derivative info (gradient / possibly Hessian) returned by the optimizer
  # opt$par — solution: numeric vector of optimized covariance params ("theta")
  # environment(devfun)$lower — lower bounds for parameters (used to detect boundary solutions).
  # control$checkConv — settings controlling tolerances etc.
  # cc is a list with items indicating whether checks passed and with messages/values helpful 
  # for diagnostic printing. e.g. list(conv = TRUE/FALSE, messages = c(...), optpar = opt$par)
  cc <- checkConv(attr(opt,"derivs"), 
                  opt$par,
                  ctrl = control$checkConv,
                  lbound = environment(devfun)$lower,
                  nobs = nrow(lmod$fr))

  mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
            mc = mcout, lme4conv=cc) ## prepare output
}


#' Parse formula and data into internal model components
#' 
#' Entry point: transforms a user-specified mixed model formula into the internal
#' model components that the rest of lme4 (mkLmerDevfun, optimizeLmer, etc.) use 
#' to fit the model numerically.
#' 
#' @param formula Mixed model formula (e.g. `y ~ x + (1 | id)`)
#' @param data data.frame containing all variables in the formula
#' @param REML logical Whether to use restricted maximum likelihood.
#' @param subset,weights,na.action,offset,contrasts Standard modeling arguments.
#' @param control a list  (from lmerControl()) giving (for \code{[g]lFormula}) all options (see \code{\link{lmerControl}} for running the model;
#' (for \code{mkLmerDevfun,mkGlmerDevfun}) options for inner optimization step;
#' (for \code{optimizeLmer} and \code{optimize[Glmer}) control parameters for nonlinear optimizer (typically inherited from the \dots argument to \code{lmerControl})
#' @param ... Extra arguments (e.g. family= for GLMMs)
#' 
#' @return \bold{lFormula, glFormula}: A list containing components,
#' \item{fr}{model frame} data.frame. Full model frame (cleaned and factorized)
#' \item{X}{fixed-effect design matrix}
#' \item{reTrms}{list containing information on random effects structure: 
#' result of \code{\link{mkReTrms}}} Random effects structure (design matrices, grouping factors, bounds, etc.)
#' \item{REML}{(lFormula only): logical flag: use restricted maximum likelihood? (Copy of argument.)}
#' \item{formula} Original full formula (Copy of argument.)
#' \item{wmsgs} named character vector Warning messages from checks
#' 
# #' @export
#' 
lFormula <- function(formula, data=NULL, REML = TRUE,
                     subset, weights, na.action, offset, contrasts = NULL,
                     control=lmerControl(), ...)
{
    # Extract and simplify the control structure
    # lmerControl() returns a nested list of control components (optCtrl, checkControl, etc.).
    # Only the checkControl part is needed at this stage (for validating the formula and data).
    control <- control$checkControl ## this is all we really need

    # Capture the function call and arguments
    # Stored in mf (for model frame construction) and mc (for recursion into glFormula if needed).
    mf <- mc <- match.call()

    # Check additional arguments and GLMM handoff
    # Exclude arguments irrelevant at this stage
    dontChk <- c("start", "verbose", "devFunOnly")
    # Extract ... and check them against allowed arguments for lmer using checkArgs()
    dots <- list(...)
    do.call(checkArgs, c(list("lmer"), dots[!names(dots) %in% dontChk]))

    # If family was passed, this isn’t a linear mixed model 
    # redirect to glFormula() (for GLMMs)
    if (!is.null(dots[["family"]])) {
        ## lmer(...,family=...); warning issued within checkArgs
        mc[[1]] <- quote(lme4::glFormula)
        if (missing(control)) mc[["control"]] <- glmerControl()
        return(eval(mc, parent.frame()))
    }

    # Validate formula and environment
    # Check the left-hand side (LHS) of the formula (y).
    # Ensure the formula environment matches the dataset and variable scope
    cstr <- "check.formula.LHS"
    checkCtrlLevels(cstr, control[[cstr]])
    denv <- checkFormulaData(formula, data,
                             checkLHS = control$check.formula.LHS == "stop")
    #mc$formula <- formula <- as.formula(formula,env=denv) ## substitute evaluated call
    formula <- as.formula(formula, env=denv)
    ## as.formula ONLY sets environment if not already explicitly set.
    ## ?? environment(formula) <- denv

    # Expand and preprocess random-effects syntax
    # get rid of || terms so update() works as expected --> Convert || (which
    # means uncorrelated random effects) into an equivalent form using | but
    # zeroing covariances.
    # Prepare formula for further processing
    RHSForm(formula) <- reformulas::expandDoubleVerts(RHSForm(formula))

    mc$formula <- formula

    # Prepare model frame call
    # Construct the call to model.frame() that will extract the relevant
    # subset of the data
    ## (DRY! copied from glFormula)
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)] # Keep essential arguments
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)

    # Build the model frame (remove random parts for now)
    fr.form <- reformulas::subbars(formula) # substitute "|" by "+" (so only fixed-effect and response variables remain)
    environment(fr.form) <- environment(formula)
    # Ensure weights and offset variables are available in the formula environment 
    # if referenced
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in '...'),
    ## so they have to be put there:
    for (i in c("weights", "offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i,get(i,parent.frame()),environment(fr.form))
    }

    # Evaluate the model frame
    # Build fr: a data frame containing the variables used in the model
    # after applying subset and NA filters
    mf$formula <- fr.form
    fr <- eval(mf, parent.frame())
    if (nrow(fr) == 0L) stop("0 (non-NA) cases")
    
    ## convert character vectors to factor (defensive)
    # Ensure character columns are converted to factors (important for random effect grouping)
    fr <- factorize(fr.form, fr, char.only=TRUE)

    # Attach metadata 
    ## store full, original formula & offset
    attr(fr,"formula") <- formula
    attr(fr,"offset") <- mf$offset
    n <- nrow(fr)

    # Extract random effect structure
    ## random effects and terms modules
    # mkReTrms() analyzes the random-effect terms ((1|group), (x|id), etc.) and constructs:
    #   Zt —> transpose of random effects model matrix
    #   flist —> factor list (grouping factors)
    #   cnms —> column names for each random effect
    #   theta, lower -> variance parameters and bounds
    # Returns reTrms: a list of matrices and indexing structures that represent the random effects design
    reTrms <- reformulas::mkReTrms(reformulas::findbars(RHSForm(formula)), fr)

    # Validate random effects structure
    # Any warnings are saved to wmsgs
    wmsgNlev <- checkNlevels(reTrms$flist, n=n, control) # Enough levels per grouping factor
    wmsgZdims <- checkZdims(reTrms$Ztlist, n=n, control, allow.n=FALSE) # Matrices have correct dimensions
    if (anyNA(reTrms$Zt)) {
        stop("NA in Z (random-effects model matrix): ",
             "please use ",
             shQuote("na.action='na.omit'"),
             " or ",
             shQuote("na.action='na.exclude'"))
    } # No NAs in Zt
    wmsgZrank <- checkZrank(reTrms$Zt, n=n, control, nonSmall = 1e6) # Sufficient rank in the random effects matrix

    # Construct fixed-effects model matrix X
    fixedform <- formula
    RHSForm(fixedform) <- reformulas::nobars(RHSForm(fixedform)) # Remove random effects from formula (nobars()).
    mf$formula <- fixedform
    ## re-evaluate model frame to extract predvars component
    fixedfr <- eval(mf, parent.frame()) # Build a fixed-effects-only model frame
    # Extract predictor variable info (used later for refitting, predictions, etc.)
    attr(attr(fr,"terms"), "predvars.fixed") <-
        attr(attr(fixedfr,"terms"), "predvars")
    ## so we don't have to fart around retrieving which vars we need
    ##  in model.frame(.,fixed.only=TRUE)
    attr(attr(fr,"terms"), "varnames.fixed") <- names(fixedfr)

    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    RHSForm(ranform) <- reformulas::subbars(RHSForm(reOnly(formula)))
    mf$formula <- ranform
    ranfr <- eval(mf, parent.frame())
    attr(attr(fr,"terms"), "predvars.random") <-
        attr(terms(ranfr), "predvars")

    ## FIXME: shouldn't we have this already in the full-frame predvars?
    # Produces a dense numeric design matrix for fixed effects
    X <- model.matrix(fixedform, fr, contrasts)#, sparse = FALSE, row.names = FALSE) ## sparseX not yet

    # Check X for rank and scaling issues
    ## backward compatibility (keep no longer than ~2015):
    if(is.null(rankX.chk <- control[["check.rankX"]]))
        rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)

    if(is.null(scaleX.chk <- control[["check.scaleX"]]))
        scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    X <- checkScaleX(X, kind=scaleX.chk)

    list(fr = fr, X = X, reTrms = reTrms, REML = REML, formula = formula,
         wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))
}

##' @param fr A model frame containing the variables needed to create an
##'   \code{\link{lmerResp}} or \code{\link{glmResp}} instance
##' @param X fixed-effects design matrix
##' @param reTrms information on random effects structure (see \code{\link{mkReTrms}})
##' @param REML (logical) fit restricted maximum likelihood model?
##' @param start starting values
##' @param verbose print output?
##' @return \bold{mkLmerDevfun, mkGlmerDevfun}: A function to calculate deviance
##' (or restricted deviance) as a function of the theta (random-effect) parameters
##' (for GlmerDevfun, of beta (fixed-effect) parameters as well).  These deviance
##' functions have an environment containing objects required for their evaluation.
##' CAUTION: The output object of \code{mk(Gl|L)merDevfun} is an \code{\link{environment}}
##' containing reference class objects (see \code{\link{ReferenceClasses}}, \code{\link{merPredD-class}},
##' \code{\link{lmResp-class}}), which behave in ways that may surprise many users. For example, if the
##' output of \code{mk(Gl|L)merDevfun} is naively copied, then modifications to the original will
##' also appear in the copy (and vice versa). To avoid this behavior one must make a deep copy
##' (see \code{\link{ReferenceClasses}} for details).
##' \cr
##' \cr
# #' @export
mkLmerDevfun <- function(fr, X, reTrms, REML = TRUE, start = NULL,
                         verbose = 0, control = lmerControl(), ...)
{
    ## FIXME: make sure verbose gets handled properly
    #if (missing(fr)) {
    ## reconstitute frame
    #}
    ## pull necessary arguments for making the model frame out of ...

    # p = number of fixed effects.
    p <- ncol(X) # maybe also do rank check on X here??
    # rho = new environment that will hold all internal objects (pp, resp, etc.)
    # This will later become the environment of devfun
    # Using an environment allows devfun to access and mutate complex objects efficiently 
    # (reference semantics).
    # Reference semantics are crucial: the internal C++ objects (merPredD, lmResp) are
    # reference classes: R environments preserve their identities
    rho <- new.env(parent=parent.env(environment()))

    # Note: having response and predictor objects separately allows modular handling of residuals,
    # weights, and likelihood computations.

    # Construct the random-effects “predictor” module -----------------------------------------------
    # Call the constructor for merPredD — a Reference Class (implemented in C++)
    # that encapsulates all the linear-algebra structures needed for mixed-model computations.
    # Fields:
    # $X        : numeric matrix (n × p) - the fixed-effect design matrix
    # $Zt       : dgCMatrix (sparseMatrix, q × n) - transposed random-effects design matrix (sparse)
    # $theta    : numeric vector of length (random-effect parameters) - initial covariance parameters
    # $Lambdat  : dgCMatrix (sparse lower triangular) - Cholesky-like factor related to random effects
    # $n        : number of observations
    # $Lind     : integer index vector - indexing structure mapping theta elements to blocks of Z
    # Methods:
    # $updateDecomp(), $solve(), $setTheta(), $getTheta(), $Ptr(), etc.

    # This object is the workhorse for evaluating the mixed model’s likelihood.
    rho$pp <- do.call(merPredD$new, c(reTrms[c("Zt","theta","Lambdat","Lind")],
                                      n=nrow(X), list(X=X)))
    
    # Set up the response module ---------------------------------------------------------------
    # REMLpass value tells the C++ layer how many parameters to include in the REML correction (since REML
    # likelihood subtracts fixed-effect degrees of freedom).
    REMLpass <- if(REML) p else 0L

    # mkRespMod() returns a a Reference class object of class "lmResp" that handles 
    # response vector y, weights, offsets, and residual calculations.
    # Fields:
    # $y      : numeric vector of responses
    # $mu     : current fitted values
    # $wt     : observation weights
    # $REML   : integer flag (0 or p)
    # $sigma  : residual scale parameter
    # Methods:
    # $updateMu(), $updateWrss(), $Ptr(), etc.
    rho$resp <-
        if(missing(fr))
             mkRespMod(    REML = REMLpass, ...)
        else mkRespMod(fr, REML = REMLpass)
    ## FIXME / note: REML does double duty as rank of X and a flag for using
    ## REML maybe this should be mentioned in the help file for
    ## mkRespMod??  currently that help file says REML is logical.  a
    ## consequence of this double duty is that it is impossible to fit
    ## a model with no fixed effects using REML (MM: ==> FIXME)
    ## devfun <- mkdevfun(rho, 0L, verbose=verbose, control=control)

    # Define the deviance function ------------------------------------------------

    pp <- resp <- NULL # prevent R CMD check false pos. warnings (in this function only)
    
    # Assign lmer_Deviance (a compiled C routine symbol) into rho
    rho$lmer_Deviance <- lmer_Deviance
    
    # Call directly into C code via .Call(), passing pointers to pp and resp’s internal C++ objects.
    # Return the deviance (negative log-likelihood) for the given theta.
    # Direct C call avoids R-level overhead in every deviance evaluation
    # The closure allows the optimizer to repeatedly call devfun(theta) without rebuilding objects

    devfun <- function(theta)
        # C++ code updates pp and resp given theta, computes the penalized residual sum of squares
        # (or REML deviance), and returns a numeric scalar.
        .Call(lmer_Deviance, pp$ptr(), resp$ptr(), as.double(theta))
    
    # Attaches rho as the environment of devfun
    environment(devfun) <- rho

    # Compute starting values if not provided -----------------------------------
    # If all random effects are of the form 1|f and starting values not otherwise
    # provided (and response variable is present, i.e. not doing a simulation),
    # compute reasonable starting values for simple random-intercept models.
    # Helps the optimizer start near a plausible region in parameter space.

    if (is.null(start) &&
        all(reTrms$cnms == "(Intercept)") &&
        length(reTrms$flist) == length(reTrms$lower) &&
        !is.null(y <- model.response(fr))) {
        # Compute variance of within-group means
        v <- sapply(reTrms$flist, function(f) var(ave(y, f)))
        # Compute residual variance
        v.e <- var(y) - sum(v)

        if (!is.na(v.e) && v.e > 0) {
            # Ratio: crude variance estimates
            v.rel <- v / v.e
            # Set pp$theta to sqrt(v.rel) if they exceed lower bounds
            if (all(v.rel >= reTrms$lower^2)) rho$pp$setTheta(sqrt(v.rel))
        }
    }

    # Initial evaluation to “warm up” -------------------------------------
    # Call devfun() once at the initial theta values
    # Ensures all internal cached values (L, Lambda, etc.) are initialized before optimization
    # Some reference-class fields are lazily populated only when first used; this “touch” 
    # ensures consistency.
    if (length(rho$resp$y) > 0)  ## only if non-trivial y
        devfun(rho$pp$theta) # one evaluation to ensure all values are set
    
    # Store lower bounds and return closure --------------------------------
    # Add vector of lower bounds for theta to the environment
    # These are needed later by the optimizer and convergence checks.
    rho$lower <- reTrms$lower # to be more consistent with mkGlmerDevfun

    # Return the devfun function — whose environment now contains:
    # - pp, resp, lower, and lmer_Deviance
    devfun # this should pass the rho environment implicitly

    # Conceptual summary: 
    # mkLmerDevfun()
    #    ↓
    # build predictor module (pp: merPredD)
    #    ↓
    # build response module (resp: lmResp)
    #    ↓
    # define devfun(theta) → .Call(lmer_Deviance, pp, resp, theta)
    #    ↓
    # set starting values if needed
    #    ↓
    # warm up with one evaluation
    #    ↓
    # return devfun

    # Key ideas:
    # Use of Reference Classes (merPredD, lmResp): 
    #       needed for efficient in-place matrix updates during optimization
    # Store everything in an environment (rho): 
    #       Provides lexical scope & reference semantics for devfun
    # Return a function (devfun) instead of an object: 
    #       Optimizers in R work with functions returning scalar objectives
    # .Call() interface to C++: 
    #       Dramatic speedup by evaluating deviance in compiled code
    # Optional REMLpass = p	
    #       Corrects likelihood for the number of fixed effects
    # Built-in heuristic for start
    #       Stable optimization starts without user input
    # One initial devfun call
    #       Initializes all internal slots (avoids NA or stale values)

    # Use
    #   lmod <- lFormula(y ~ x + (1|g), data = df)
    #   devfun <- do.call(mkLmerDevfun, c(lmod, list(REML = TRUE)))

    # Content
    # str(environment(devfun))
    # List of 4
    # $ pp             :RefClass 'merPredD' [package "lme4"] - predictor / random effects machinery
    # $ resp           :RefClass 'lmResp' [package "lme4"] - response data and residuals
    # $ lmer_Deviance  : symbol - external pointer / symbol	compiled C++ function to compute deviance
    # $ lower          : num [1:n] 0 0 0 - lower bounds for theta

    # Because they are reference-based, if you modify them in place ($setTheta() etc.), 
    # those changes propagate — copies of devfun will share the same underlying state 
    # unless you deep-copy the environment.

}

#' optimizeLmer actually minimizes that devfun(theta) to find the best covariance 
#' parameters.
#' 
#' It’s essentially the “fitting engine” for random-effect covariance parameters.
#' 
#' It chooses a numeric optimizer (like bobyqa, nloptwrap, optim), runs it, checks 
#' if the fit hits parameter boundaries, possibly restarts optimization from a better
#' starting point, and returns the result.
#' 
#' @param devfun a deviance function, as generated by \code{\link{mkLmerDevfun}}
#' @param optimizer character or list. Which optimizer to use (e.g. "bobyqa", "Nelder_Mead")
#' @param restart_edge logical. Whether to restart optimization if solution is near a 
#'  parameter boundary
#' @param boundary.tol numeric. Tolerance to check if solution lies on the boundary.
#' @param start numeric vector or NULL. Starting values for θ
#' @param verbose integer. Level of printed output	
#' @param control Optimizer control list
#' @param ... Extra arguments passed through to the optimizer
#' 
##' @return \bold{optimizeLmer}: Results of an optimization.
optimizeLmer <- function(devfun,
                         optimizer=    formals(lmerControl)$optimizer,
                         restart_edge= formals(lmerControl)$restart_edge,
                         boundary.tol = formals(lmerControl)$boundary.tol,
                         start = NULL,
                         verbose = 0L,
                         control = list(),
                         ...) {
    
    verbose <- as.integer(verbose)

    # Extract the environment of the deviance function
    rho <- environment(devfun)

    # optwrap() unifies different optimizers under one consistent API.
    opt <- optwrap(optimizer,
                   devfun,
                   # Use either the user-supplied starting values or extract 
                   # the current θ from pp.
                   # If start is NULL, returns rho$pp$theta.
                   # If non-NULL, validates the provided vector.
                   getStart(start, rho$pp),
                   lower=rho$lower,
                   control=control,
                   adj=FALSE, verbose=verbose,
                   ...)
    # Retunrs a list-like structure 
    # opt$par   # optimized theta values (numeric vector)
    # opt$fval  # minimum deviance (scalar)
    # opt$conv  # convergence code (0 = success)
    # opt$message  # any optimizer messages
    # attr(opt, "derivs")  # (optional) gradient/Hessian info
    
    # If restart_edge = TRUE, check if any variance parameters θ are exactly at 
    # their lower bounds (usually 0).
    # If so, restart optimization from that slightly perturbed point.
    # Mixed models often converge to boundary values (zero variance components).
    # Some of those are genuine; others are numerical artifacts.
    # Restarting helps escape false boundaries caused by flat likelihood surfaces.

    if (restart_edge) {
        ## FIXME: should we be looking at rho$pp$theta or opt$par
        ##  at this point???  in koller example (for getData(13)) we have
        ##   rho$pp$theta=0, opt$par=0.08

        # check if any variance parameters (θ) are exactly at their lower bounds (usually 0).
        if (length(bvals <- which(rho$pp$theta==rho$lower)) > 0) {
            # bvals = vector of indices of of theta values that equal their lower bound (e.g., 0).

            ## *don't* use numDeriv -- cruder but fewer dependencies, no worries
            ##  about keeping to the interior of the allowed space
            
            # Create a 'deep' copy so modifying theta inside the test won't affect the 
            # reference object.
            theta0 <- new("numeric",rho$pp$theta)

            # Evaluate current deviance
            d0 <- devfun(theta0)

            # Define small perturbation
            btol <- 1e-5  ## FIXME: make user-settable?

            # Compute crude numerical gradients near the boundary to see if pushing
            # slightly inside lowers deviance
            bgrad <- sapply(bvals,
                            function(i) {
                                bndval <- rho$lower[i]
                                theta <- theta0
                                theta[i] <- bndval+btol
                                # Approximate the derivative (slope) of the deviance function just inside 
                                # the boundary (for each boundary parameter i)
                                (devfun(theta)-d0)/btol
                            })
            
            # Reset internal state
            devfun(theta0) # Ensure that after these gradient checks, rho$pp$theta is reset to original

            ## FIXME: allow user to specify ALWAYS restart if on boundary?
            if (any(is.na(bgrad))) {
                warning("some gradient components are NA near boundaries, skipping boundary check")
                return(opt)
            } else {
                # If gradient is negative, deviance decreases when moving inside the feasible
                # region — meaning the optimizer likely got stuck at boundary prematurely.
                if (any(bgrad < 0)) {
                    if (verbose) message("some theta parameters on the boundary, restarting")
                    # Restart optimization using current optimum opt$par as starting point.
                    opt <- optwrap(optimizer,
                                   devfun,
                                   opt$par,
                                   lower=rho$lower, control=control,
                                   adj=FALSE, verbose=verbose,
                                   ...)
                }
            } ## bgrad not NA
        }
    } ## if restart.edge

    # Boundary check (tolerance)
    # If a nonzero boundary.tol is set, call check.boundary() to verify whether the 
    # solution lies within tolerance of the boundary (and possibly issue warnings).
    # Protects against degenerate or singular fits (e.g., variance = 0 or correlation = ±1).
    # This check is purely diagnostic — doesn’t modify the fit.
    if (boundary.tol > 0)
        check.boundary(rho, opt, devfun, boundary.tol)
    else
        opt
}


#' Calls the numerical optimizer (like "bobyqa", "Nelder_Mead", "optimx", or "nloptwrap") 
#' and standardizes their behavior

#' This is a unified wrapper that:
#' Prepares optimizer-specific parameters (control tweaks),
#' Calls the chosen optimizer (via getOptfun()),
#' Collects results,
#' Standardizes convergence messages and derivative calculations,
#' Returns a consistent object for higher-level use (optimizeLmer()).
#' 
#' @param optimizer	character or function. Which optimizer to use (e.g. "bobyqa", "Nelder_Mead", "nloptwrap", "optimx")
#' @param fn function. The deviance function to minimize (usually from mkLmerDevfun())
#' @param par numeric vector. Initial parameter values (e.g. θ)
#' @param lower,upper numeric vector. Parameter bounds.
#' @param control list. Optimizer control options
#' @param adj logical. “Adjustment” flag for second optimization pass (e.g. in glmer/nlmer)
#' @param calc.derivs logical. Whether to compute gradients/Hessian after fit
#' @param use.last.params logical. Whether to reuse last internal θ values from the devfun environment
#' @param verbose integer. Verbosity level.
optwrap <- function(optimizer, fn, par, lower = -Inf, upper = Inf,
                    control = list(), adj = FALSE, calc.derivs = TRUE,
                    use.last.params = FALSE,
                    verbose = 0L)
{
    ## control must be specified if adj==TRUE;
    ##  otherwise this is a fairly simple wrapper
    
    # Retrieve the actual R function corresponding to the chosen optimizer name.
    # e.g.:
    # "bobyqa" → bobyqa() ---> not a correct mapping just an example
    # "Nelder_Mead" → Nelder_Mead()
    # "nloptwrap" → nloptwrap()
    optfun <- getOptfun(optimizer)

    # Extract a human-readable optimizer name
    # If the user passed a function directly, deparse() extracts its symbol.
    optName <- if(is.character(optimizer)) optimizer
    else ## "good try":
        deparse(substitute(optimizer))[[1L]]
    
    # Ensure bounds have correct length
    lower <- rep_len(lower, length(par))
    upper <- rep_len(upper, length(par))

    # “Adjustment”: modify control parameters for a second optimization phase (used in nlmer or glmer).
    # if (adj)
    #     ## control parameter tweaks: only for second round in nlmer, glmer
    #     switch(optName,
    #            # Tweak trust region parameters rhobeg and rhoend (controls exploration step sizes).
    #            "bobyqa" = {
    #                if(!is.numeric(control$rhobeg)) control$rhobeg <- 0.0002
    #                if(!is.numeric(control$rhoend)) control$rhoend <- 2e-7
    #            },
    #            "Nelder_Mead" = {
    #                if (is.null(control$xst))  {
    #                    thetaStep <- 0.1
    #                    nTheta <- length(environment(fn)$pp$theta)
    #                    betaSD <- sqrt(diag(environment(fn)$pp$unsc()))
    #                    control$xst <- 0.2* c(rep.int(thetaStep, nTheta),
    #                                          pmin(betaSD, 10))
    #                }
    #                if (is.null(control$xt)) control$xt <- control$xst*5e-4
    #            })
    
    # Handle verbosity per optimizer
    # Different optimizers have different control argument names for verbosity.
    # Map them on to the right control field:
    # bobyqa → iprint
    # Nelder_Mead → verbose
    # nloptwrap → print_level
    # Adds a warning for unknown optimizers that ignore verbose.
    switch(optName,
           "bobyqa" = {
               if(all(par == 0)) par[] <- 0.001  ## minor kludge
               if(!is.numeric(control$iprint)) control$iprint <- min(verbose, 3L)
           },
           "Nelder_Mead" = control$verbose <- verbose,
           "nloptwrap" = control$print_level <- min(as.numeric(verbose),3L),
           ## otherwise:
           if(verbose) warning(gettextf(
               "'verbose' not yet passed to optimizer '%s'; consider fixing optwrap()",
                                        optName), domain = NA)
           )
    
    # Build a unified list of arguments that will be passed to the actual optimizer via do.call().
    arglist <- list(fn = fn, par = par, lower = lower, upper = upper, control = control)
    
    # Special handling for optimx
    ## optimx: must pass method in control (?) because 'method' was previously
    ## used in lme4 to specify REML vs ML
    if (optName == "optimx") {
        if (is.null(method <- control$method))
            stop("must specify 'method' explicitly for optimx")
        arglist$control$method <- NULL
        arglist <- c(arglist, list(method = method))
    }

    # Catch and record warnings during optimization
    # Call the optimizer using do.call(optfun, arglist).
    # Optimizers indirectly modify internal objects of fn (deviance function
    # evaluating fn(par).
    # Any warnings raised during optimization are caught by withCallingHandlers() and stored in curWarnings.
    # so they can be later included in the output.

    ## FIXME: test!  effects of multiple warnings??
    ## may not need to catch warnings after all??
    curWarnings <- list()
    opt <- withCallingHandlers(do.call(optfun, arglist),
                               warning = function(w) {
                                   curWarnings <<- append(curWarnings,list(w$message))
                               })
    ## cat("***",unlist(tail(curWarnings,1)))
    ## FIXME: set code to warn on convergence !=0

    # Post-processing optimizer-specific quirks
    ## post-fit tweaking
    if (optName == "bobyqa") {
        # bobyqa() returns its convergence code in ierr; copy it to a common field convergence
        opt$convergence <- opt$ierr
    } else if (optName == "Nelder_Mead") {
        # Nelder_Mead reports “ran out of iterations” as convergence — code 4 is used for that.
        ## fix-up: Nelder_Mead treats running out of iterations as "convergence" (!?)
        if (opt$NM.result==4) opt$convergence <- 4
    } else if (optName == "optimx") {
        # Normalize optimx’s complex matrix-style output into a simple list.
        opt <- list(par = coef(opt)[1,],
                    fvalues = opt$value[1],
                    method = method,
                    conv = opt$convcode[1],
                    feval = opt$fevals + opt$gevals,
                    message = attr(opt,"details")[,"message"][[1]])
    }

    # Handle convergence warnings
    # Use getConv() helper to extract the convergence code (e.g., 0 = OK, nonzero = warning).
    # Builds a message via getMsg() and issues a warning (records the message in curWarnings)
    if ((optconv <- getConv(opt)) != 0) {
        wmsg <- paste("convergence code",optconv,"from",optName)
        if (!is.null(getMsg(opt))) wmsg <- paste0(wmsg,": ",getMsg(opt))
        warning(wmsg)
        curWarnings <<- append(curWarnings,list(wmsg))
    }
    ## pp_before <- environment(fn)$pp
    ## save(pp_before,file="pp_before.RData")

    # Compute derivatives (optional)
    # Needed for later convergence checks (checkConv()).
    if (calc.derivs) {
        if (use.last.params) {
            # If use.last.params, first save and restores the function’s internal parameter
            # state (since the deviance function’s environment holds reference-class objects
            # that can change in place).
            ## +0 tricks R into doing a deep copy ...
            ## otherwise element of ref class changes!
            ## FIXME:: clunky!!
            orig_pars <- opt$par
            orig_theta <- environment(fn)$pp$theta+0
            orig_pars[seq_along(orig_theta)] <- orig_theta
        }
        if (verbose > 10) cat("computing derivatives\n")
        # Compute the gradient/Hessian at the optimum
        derivs <- deriv12(fn, opt$par, fx = opt$value)

        if (use.last.params) {
            # The fn() evaluation restores internal environment consistency 
            # after derivative computation.
            ## run one more evaluation of the function at the optimized
            ##  value, to reset the internal/environment variables in devfun ...
            fn(orig_pars)
        }
    } else derivs <- NULL

    # Reset function environment if needed
    # Ensure that the internal pp and resp objects in fn’s environment
    # are synchronized with the final parameter vector.
    if (!use.last.params) {
        ## run one more evaluation of the function at the optimized
        ##  value, to reset the internal/environment variables in devfun ...
        fn(opt$par)
    }

    # Return results
    # par	Optimized parameters
    # value or fval	Objective (deviance) value
    # convergence	Convergence code
    # message	Optimizer message
    # optimizer	Name or function used
    # control	Final control parameters
    # warnings	List of warnings caught
    # derivs	Gradient/Hessian info if computed
    structure(opt, ## store all auxiliary information
              optimizer = optimizer,
              control   = control,
              warnings  = curWarnings,
              derivs    = derivs)
}

#' Wrapper around the nloptr package optimizer that’s used internally 
#' to perform parameter optimization when fitting mixed-effects models.
#' 
#' This defines the function inside a local() environment so that the 
#' internal variable defaultControl is private — it won’t pollute the 
#' global environment.
#' This is a common R trick to give a function persistent default data.
#' After local() executes, the result is a single function (the one defined
#'  inside it), but the enclosing environment keeps defaultControl around
#'  invisibly.
#' 
nloptwrap <- local({
    # Define default optimizer settings
    ## define default control values in environment of function ...
    defaultControl <- list(algorithm="NLOPT_LN_BOBYQA", # Which optimization algorithm to use, a derivative-free local optimizer.
                           xtol_abs=1e-8, # Absolute tolerance on parameter changes — when small, optimization stops.
                           ftol_abs=1e-8, # Absolute tolerance on function value changes — another stopping criterion.
                           maxeval=1e5) # Maximum number of function evaluations.
    ##
    function(par, fn, lower, upper, control=list(),...) {
        # Merge defaults with user-supplied control
        # Check each name in defaultControl (algorithm, xtol_abs, etc.)
        # If the user didn't provide a value for that element in control, 
        # it fills it in from the defaults.
        for (n in names(defaultControl))
            if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
        
        # Call the underlying optimizer
        res <- nloptr(x0=par, # starting parameter values (vector).
            eval_f=fn, # the deviance function.
            lb=lower, ub=upper, # lower and upper bounds.
            opts=control, # control list containing the options.
            ...)
        # nloptr() returns a list containing fields like:
        # Field	Meaning
        #   solution	Vector of optimized parameter values.
        #   objective	Function value at optimum.
        #   iterations	Number of iterations or evaluations.
        #   status	    Integer status code (0 = success, >0 = okay, <0 = error).
        #   message	    Human-readable status message.
        with(res, list(par   = solution,
                       fval  = objective,
                       feval = iterations,
                       ## ?nloptr: "integer value with the status of the optimization (0 is success)"
                       ## most status>0 are fine (e.g. 4 "stopped because xtol_rel was reached"
                       ## but status 5 is "ran out of evaluations"
                       conv  = if (status<0 || status==5) status else 0,
                       message = message))
    }
})

