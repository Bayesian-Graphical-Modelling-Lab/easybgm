#' @title Bayesian analysis of networks
#'
#' @description Easy estimation of a Bayesian graphical models to obtain
#'   conditional (in)dependence relations between variables in a network.
#'
#' @name easybgm
#'
#' @param data An n x p matrix or dataframe containing the variables for n
#'   independent observations on p variables.
#'
#' @param type Specifies the type of data. There are two ways to use this argument:
#'
#'   \strong{1. A single string}:
#'   \itemize{
#'     \item \code{"continuous"}: For continuous (Gaussian) data. Estimates a
#'       Gaussian Graphical Model (GGM).
#'     \item \code{"ordinal"}: For ordinal (Likert-type) data. Estimates an
#'       Ordinal Markov Random Field (OMRF).
#'     \item \code{"binary"}: For binary (0/1) data. Estimates an Ising model. 
#'     \item \code{"blume-capel"}: For Blume-Capel ordinal data. Requires a
#'       reference category via the \code{baseline_category} argument.
#'    \item \code{"mixed"}: For data with both continuous and discrete
#'       variables. Requires the \code{not_cont} argument to indicate which
#'       variables are not continuous (see below). This value should only be specified when
#'       fitting a mixed model using \code{package = "BDgraph"} or \code{package = "BGGM"}.
#'       See below on how to specify mixed models using bgms.
#'   }
#'
#'   \strong{2. A character vector of length p} (per-variable specification,
#'   bgms only): Each element specifies the type of the corresponding column in
#'   \code{data}. Valid values are \code{"ordinal"}, \code{"continuous"},
#'   \code{"blume-capel"}, and \code{"binary"}.
#'
#'   For example: \code{type = c("ordinal", "ordinal", "continuous")} specifies
#'   that the first two columns are ordinal and the third is continuous.
#'
#'   \emph{Note:} Per-variable type vectors are only supported with the
#'   \code{bgms} package.
#'
#' @param package The R-package used for fitting the network model. Optional;
#'   if not specified, \code{bgms} is used as the default for all data types.
#'   Supported options:
#'   \itemize{
#'     \item \code{"bgms"} (Default): Supports the following data types: continuous,
#'       ordinal, binary, and blume-capel or a combination of the three.
#'     \item \code{"BDgraph"}: Supports continuous (fits a GGM), mixed (fits a
#'       GCGM). For continuous data, missing values are not
#'       allowed; use \code{na.omit()} on the data first. 
#'     \item \code{"BGGM"}: Supports continuous and mixed data.
#'   }
#'
#' @param iter Number of iterations for the sampler. The default depends on the
#'   package:
#'   \itemize{
#'     \item \code{bgms}: 1e3 (1,000 iterations)
#'     \item \code{BDgraph}: 1e4 (10,000 iterations)
#'     \item \code{BGGM}: 1e4 (10,000 iterations)
#'   }
#'   The recommended number of iterations depends on the data, model complexity,
#'   and desired precision. Check the convergence diagnostics in the output to
#'   determine if more iterations are needed.
#'
#' @param save Logical. Should the posterior samples be obtained
#'   (default = \code{FALSE})? If \code{TRUE}, the output includes a
#'   \code{samples_posterior} matrix with the posterior samples for each edge
#'   weight parameter. Setting \code{centrality = TRUE} automatically sets
#'   \code{save = TRUE}.
#'
#' @param centrality Logical. Should the centrality measures be extracted
#'   (default = \code{FALSE})? Note that this will significantly increase
#'   computation time. Automatically sets \code{save = TRUE}.
#'
#' @param progress Logical. Should a progress bar be shown
#'   (default = \code{TRUE})?
#'   
#' @param not_cont A binary vector of length p, required when
#'   \code{type = "mixed"} and a single string is used. Each element indicates
#'   whether the corresponding variable is not continuous
#'   (\code{1} = not continuous/ordinal, \code{0} = continuous). This is used to
#'   construct a per-variable type vector internally. Only neccesary when 
#'   using \code{package = "BDgraph"} or \code{package = "BGGM"}. 
#'   
#' @param baseline_category Integer or vector, required if at least one variable 
#'   is of type \code{"blume-capel"}. Baseline category used in
#'   Blume--Capel variables. Can be a single integer (applied to all) or a
#'   vector of length \code{p}. 
#'
#' @param ... Additional arguments passed to the fitting functions of the
#'   underlying packages (e.g., prior specifications). See the
#'   \strong{Prior specification} section in Details for available prior options per
#'   package and the package help files for all other potential arguments.
#'
#' @return An object of class \code{easybgm} containing the following elements:
#'
#' \strong{Always returned:}
#' \itemize{
#'   \item \code{parameters}: A p x p matrix of posterior mean partial
#'     association estimates.
#'   \item \code{inc_probs}: A p x p matrix of posterior inclusion
#'     probabilities.
#'   \item \code{inc_BF}: A p x p matrix of posterior inclusion Bayes factors.
#'   \item \code{structure}: A p x p adjacency matrix of the median probability
#'     model (edges with posterior inclusion probability > 0.5).
#'   \item \code{model}: A string indicating the model type (e.g.,
#'     \code{"continuous"}, \code{"ordinal"}, \code{"mixed"}).
#'   \item \code{thresholds}: Threshold/intercept parameters (bgms only). The
#'     format depends on the model type: a matrix for ordinal/binary models,
#'     \code{NULL} for continuous models, or a list for mixed models.
#' }
#'
#' \strong{Returned for bgms and BDgraph:}
#' \itemize{
#'   \item \code{structure_probabilities}: Posterior probabilities of all
#'     visited graph structures (values between 0 and 1).
#'   \item \code{graph_weights}: Number of times each graph structure was
#'     visited.
#'   \item \code{sample_graphs}: Identifiers for each visited graph structure.
#' }
#'
#' \strong{Returned for bgms only:}
#' \itemize{
#'   \item \code{convergence_parameter}: The Gelman-Rubin (R-hat) convergence
#'     statistic for each edge weight parameter. Values close to 1 indicate good
#'     convergence.
#'   \item \code{MCSE_BF}: A matrix with the 95 percent Monte Carlo confidence
#'     interval for each inclusion Bayes factor.
#' }
#'
#' \strong{Returned when edge_prior = "Stochastic-Block" (bgms only):}
#' \itemize{
#'   \item \code{sbm}: A list containing Stochastic Block Model results, including
#'     \code{posterior_num_blocks} (posterior probabilities for each number of
#'     clusters), \code{posterior_mean_allocations} (posterior mean cluster
#'     assignments), \code{posterior_mode_allocations} (posterior mode cluster
#'     assignments), and \code{posterior_mean_coclustering_matrix} (a p x p matrix
#'     of pairwise co-clustering proportions).
#' }
#'
#' \strong{Interpretable parameter scales (bgms only):}
#'
#' In addition to the raw pairwise interaction parameters in
#' \code{parameters}, the following transformations are provided when the
#' model type supports them. They are \code{NULL} otherwise.
#' \itemize{
#'   \item \code{partial_correlations}: A matrix of posterior mean partial
#'     correlations. Available for continuous (GGM) models (full p x p matrix)
#'     and for the continuous block of mixed models. \code{NULL} for ordinal
#'     models.
#'   \item \code{precision_matrix}: The posterior mean precision
#'     (inverse covariance) matrix. Same availability as partial correlations.
#'   \item \code{log_odds}: A matrix of posterior mean log adjacent-category
#'     odds ratios. Available for ordinal/binary models (full p x p matrix)
#'     and for the discrete block of mixed models. \code{NULL} for continuous
#'     models.
#' }
#'
#' \strong{Returned when save = TRUE:}
#' \itemize{
#'   \item \code{samples_posterior}: A k x iter matrix of posterior samples for
#'     each edge weight parameter (k = p*(p-1)/2 edges).
#' }
#'
#' \strong{Returned when centrality = TRUE:}
#' \itemize{
#'   \item \code{centrality}: An iter x p matrix of centrality values for each
#'     node at each iteration.
#' }
#'
#' @details
#'
#' \strong{Data types and package support}
#'
#' The table below summarizes which data types are supported by each backend
#' package:
#'
#' \tabular{lccc}{
#'   \strong{Data type}    \tab \strong{bgms} \tab \strong{BDgraph} \tab \strong{BGGM} \cr
#'   continuous             \tab Yes (default)  \tab Yes              \tab Yes           \cr
#'   ordinal                \tab Yes (default)  \tab Yes              \tab No            \cr
#'   binary                 \tab Yes (always)   \tab No               \tab No            \cr
#'   mixed                  \tab Yes (default)  \tab Yes              \tab Yes           \cr
#'   blume-capel            \tab Yes (default)  \tab No               \tab No            \cr
#' }
#'
#' \strong{How \code{type} and \code{not_cont} relate}
#'
#' There are two ways to specify mixed-type data:
#' \enumerate{
#'   \item \strong{Using \code{type = "mixed"} with \code{not_cont}}: Set
#'     \code{type = "mixed"} and provide a binary vector \code{not_cont} of
#'     length p. Internally, this is translated to a per-variable type vector
#'     where \code{not_cont == 1} maps to \code{"ordinal"} and
#'     \code{not_cont == 0} maps to \code{"continuous"}.
#'   \item \strong{Using a per-variable type vector} (bgms only): Directly
#'     specify the type of each variable, e.g.,
#'     \code{type = c("ordinal", "continuous", "ordinal")}. This is more
#'     flexible and does not require the \code{not_cont} argument.
#' }
#'
#' \strong{Prior specification}
#'
#' Users may wish to deviate from the default (uninformative) prior
#' specifications. This can be done by passing additional arguments via
#' \code{...} to the fitting function of the chosen package. We give an
#' overview of the available prior arguments per package below.
#'
#' \emph{bgms} (>= 0.2.0.0): the preferred interface uses prior-constructor
#' objects from the \code{bgms} package. Pass them through \code{...}:
#' \itemize{
#'   \item \code{interaction_prior}: A parameter prior on pairwise interactions.
#'     Use \code{\link[bgms]{cauchy_prior}(scale)} (default
#'     \code{cauchy_prior(scale = 1)}), \code{\link[bgms]{normal_prior}(scale)},
#'     or \code{\link[bgms]{beta_prime_prior}(alpha, beta)}.
#'     For example, a cauchy prior with scale 1 would be specified with adding the 
#'     argument \code{interaction_prior = cauchy_prior(1)} to the easybgm call. 
#'   \item \code{threshold_prior}: A parameter prior on threshold (main effect)
#'     parameters. Use \code{\link[bgms]{beta_prime_prior}(alpha, beta)}
#'     (default \code{beta_prime_prior(0.5, 0.5)}),
#'     \code{\link[bgms]{cauchy_prior}(scale)}, or
#'     \code{\link[bgms]{normal_prior}(scale)}.
#'     For example, a cauchy prior with scale 1 would be specified with adding the 
#'     argument \code{threshold_prior = cauchy_prior(1)} to the easybgm call.
#'   \item \code{means_prior}: A prior on the means of continuous variables in
#'     mixed MRF models. Default \code{normal_prior(scale = 1)}.
#'   \item \code{precision_scale_prior}: A prior on the diagonal entries of the
#'     precision matrix (GGM and mixed MRF). Use
#'     \code{\link[bgms]{gamma_prior}(shape, rate)} (default) or
#'     \code{\link[bgms]{exponential_prior}(rate)}.
#'   \item \code{edge_prior}: An indicator prior on edge inclusion. Use
#'     \code{\link[bgms]{bernoulli_prior}(inclusion_probability)} (default
#'     \code{bernoulli_prior(0.5)}; \code{inclusion_probability} can also be a
#'     symmetric \eqn{p \times p} matrix of edge-specific probabilities),
#'     \code{\link[bgms]{beta_bernoulli_prior}(alpha, beta)}, or
#'     \code{\link[bgms]{sbm_prior}(alpha, beta, alpha_between, beta_between, dirichlet_alpha, lambda)}
#'     for the Stochastic Block Model prior.
#'     For example, a bernoulli prior with prior probabilit of 0.5 would be 
#'     specified with adding the argument \code{edge_prior = bernoulli_prior(0.5)} 
#'     to the easybgm call.
#' }
#'
#'
#' \emph{BDgraph}:
#' \itemize{
#'   \item \code{df.prior}: Degrees of freedom of the prior G-Wishart
#'     distribution on the precision matrix. Default is 3.
#'   \item \code{g.prior}: Prior probability of edge inclusion. Can be a
#'     scalar (same for all edges) or a matrix (edge-specific). 
#'     This can also be a symmetric pxp matrix of edge-specific inclusion probabilities.
#'     Default is 0.5.
#' }
#'
#' \emph{BGGM}:
#' \itemize{
#'   \item \code{prior_sd}: Standard deviation of the prior on interaction
#'     parameters (approximately the scale of a beta distribution). Default is
#'     0.25.
#' }
#'
#' We encourage researchers to conduct prior sensitivity checks.
#'
#' @export
#'
#' @import bgms
#' @importFrom BDgraph bdgraph bdgraph.mpl plinks
#' @importFrom BGGM explore select
#' @importFrom utils packageVersion
#'
#' @examples
#'
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(Wenchuan)[1:100, 1:5]
#'
#' # --- Continuous data (default: bgms) ---
#' fit <- easybgm(data, type = "continuous",
#'                 iter = 100  # for demonstration only; increase for real analyses
#'                 )
#' summary(fit)
#'
#' \donttest{
#' # --- Mixed data using per-variable type vector (bgms only) ---
#' dat3 <- data[, 1:3]
#' fit_vec <- easybgm(dat3,
#'                     type = c("ordinal", "ordinal", "continuous"),
#'                     iter = 100)
#'
#' # --- Extract posterior samples and centrality ---
#' fit_full <- easybgm(data, type = "continuous",
#'                      iter = 100,
#'                      centrality = TRUE, save = TRUE)
#'
#' # --- Using BDgraph for continuous data ---
#' fit_bd <- easybgm(data, type = "continuous",
#'                    package = "BDgraph",
#'                    iter = 100)
#'
#' # --- Using BGGM for continuous data ---
#' fit_bggm <- easybgm(data, type = "continuous",
#'                      package = "BGGM",
#'                      iter = 100)
#'                    
#' }



easybgm <- function(data, type, package = NULL, 
                    save = FALSE, centrality = FALSE,
                    iter = 1e3, progress = TRUE, 
                    baseline_category = NULL, not_cont = NULL,
                    ...){
  
  # --- Handle vector type (per-variable specification) ---
  # When type is a vector of length > 1 (e.g., c("ordinal", "ordinal", "continuous")),
  # it specifies the variable type for each column. This is only supported with bgms.
  is_vector_type <- length(type) > 1
  
  if(is_vector_type) {
    valid_types <- c("ordinal", "continuous", "blume-capel", "binary")
    invalid <- type[!type %in% valid_types]
    if(length(invalid) > 0) {
      warning("The following variable type(s) are not recognized: ",
              paste0("'", unique(invalid), "'", collapse = ", "), ". ",
              "Valid types are: ", paste(valid_types, collapse = ", "), ". ",
              "Please check for typos.",
              call. = FALSE)
      stop("Invalid variable types detected. See the warning message for more details.",
           call. = FALSE)
    }
    if(length(type) != ncol(data)) {
      stop("When 'type' is a vector, its length (", length(type), ") must equal ",
           "the number of columns in 'data' (", ncol(data), ").",
           call. = FALSE)
    }
    # Map "binary" to "ordinal" internally
    type[type == "binary"] <- "ordinal"
    package <- "package_bgms"
  }
  
  if(!is_vector_type && length(type) == 1 && type == "mixed" && is.null(not_cont)){
    stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
         call. = FALSE)
  }
  
  dots <- list(...)
  # reference_category still needs to be checked for backwards compatability
  has_reference <- "reference_category" %in% names(dots)
  has_baseline  <- !is.null(baseline_category)
  
  # If type contains "blume-capel", a reference category must be present
  if (any(type == "blume-capel") && !(has_reference || has_baseline)) {
    stop("For the Blume-Capel model, a reference category needs to be specified.
         If type is 'blume-capel' it specifies the reference category in the Blume-Capel model.
         Should be an integer within the range of integer scores observed for the
         'blume-capel' variable. Can be a single number specifying the reference
         category for all Blume-Capel variables at once, or a vector of length
         p. bgm ignores its elements for other variable types. 
         For bgms version smaller than 0.1.6, use the reference_category argument.
         For all package versions including and older than 0.1.6., the baseline_category argument.",
         call. = FALSE)
  }
  
  
  # Set default values for fitting if package is unspecified
  if(any(type == "blume-capel")){
    package <- "package_bgms"
  } else if(is_vector_type){
    package <- "package_bgms"
  } else if(is.null(package)){
    package <- "package_bgms"
  } else {
    if(package == "BDgraph") package <- "package_bdgraph"
    if(package == "BGGM") package <- "package_bggm"
    if(package == "bgms") package <- "package_bgms"
  }
  
  # change the default number of iterations depending on the underlying package
  if(iter == 1e3 && package == "package_bdgraph"){
    iter <- 1e4
  } else if (iter == 1e3 && package == "package_bggm"){
    iter <- 1e4
  }
  
  if(length(type) == 1 && type == "continuous" && package == "package_bdgraph" && any(is.na(data))){
    stop("The data contains missing values which cannot be handled as continuous data by BDgraph (GGM). ",
         "Please either:\n",
         "  1) Remove missing values first (e.g., data <- na.omit(data)), or\n",
         "  2) Set type = 'mixed' to estimate a GCGM, which can handle missing data.",
         call. = FALSE)
  }
  
  
  fit <- list()
  class(fit) <- c(package, "easybgm")
  
  if(!save && centrality){
    save <- TRUE
  }
  
  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, not_cont = not_cont, iter = iter,
                    save = save, centrality = centrality, progress = progress, 
                    baseline_category = baseline_category, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop(paste("Error meassage: ", e$message, "Please consult the original message for more information.") )
    })
  
  
  # Extract the results
  res <- bgm_extract(fit, type = type,
                     save = save, not_cont = not_cont,
                     data = data, centrality = centrality, 
                     iter = iter,
                     ...)

  if(any(class(res) == "package_bgms")){
    if(any(res$convergence_parameter > 1.01, na.rm = TRUE)){
      warning("One or more of the convergence statistics are larger than 1.01.
              These values are considered concerning, indicating potential lack
              of convergence for the estimates of the pairwise interactions.
              Try fitting the model with more iterations of the sampler (iter).",
              call. = FALSE)
    }
  }
  
  # Output results
  class(res) <- c(package, "easybgm")
  return(res)
}
