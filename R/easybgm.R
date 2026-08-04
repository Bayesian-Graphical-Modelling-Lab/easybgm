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
#'   \strong{1. A single string}, applied to every variable:
#'   \itemize{
#'     \item \code{"continuous"}: For continuous (Gaussian) data. Estimates a
#'       Gaussian Graphical Model (GGM).
#'     \item \code{"ordinal"}: For ordinal (Likert-type) data. Estimates an
#'       Ordinal Markov Random Field (OMRF).
#'     \item \code{"binary"}: For binary (0/1) data. Estimates an Ising model.
#'     \item \code{"blume-capel"}: For Blume-Capel ordinal data. Requires a
#'       reference category via the \code{baseline_category} argument.
#'     \item \code{"mixed"}: For data with both continuous and discrete
#'       variables. Requires the \code{not_cont} argument to indicate which
#'       variables are not continuous.
#'   }
#'
#'   \strong{2. A character vector of length p} (per-variable specification):
#'   Each element gives the type of the corresponding column of \code{data},
#'   drawn from \code{"ordinal"}, \code{"continuous"}, \code{"blume-capel"},
#'   and \code{"binary"}. For example,
#'   \code{type = c("ordinal", "ordinal", "continuous")} specifies that the
#'   first two columns are ordinal and the third is continuous.
#'
#'   Per-variable vectors are fitted by \code{bgms} and require \code{bgms}
#'   version 0.2.0.0 or later. They are the most flexible way to specify mixed
#'   data and make \code{not_cont} unnecessary.
#'
#'   Which values are legal depends on the fitting package rather than on the
#'   \code{bgms} version directly; see \strong{Data types and package support}
#'   in the Details section below.
#'
#' @param package The R-package used for fitting the network model. Optional.
#'   Supported options:
#'   \itemize{
#'     \item \code{"bgms"}: Fits ordinal, binary, and blume-capel data, and --
#'       from \code{bgms} version 0.2.0.0 onwards -- continuous data and
#'       per-variable type vectors as well.
#'     \item \code{"BDgraph"}: Fits continuous data (a GGM), mixed and ordinal
#'       data (a GCGM), and binary data (a discrete graphical model). For
#'       continuous data, missing values are not allowed; use \code{na.omit()}
#'       on the data first.
#'     \item \code{"BGGM"}: Fits continuous, mixed, ordinal, and binary data.
#'   }
#'
#'   If \code{package} is not specified, \code{bgms} is used for all data types
#'   from \code{bgms} version 0.2.0.0 onwards. With \code{bgms} 0.1.6.3, which
#'   cannot fit continuous data, \code{bgms} is used for binary, ordinal, and
#'   blume-capel data and \code{BGGM} for continuous and mixed data.
#'
#'   An explicit \code{package} is honoured for every data type that package
#'   supports. Neither \code{BGGM} nor \code{BDgraph} implements the Blume-Capel
#'   model or accepts a per-variable type vector, so requesting either package
#'   for those produces a warning and \code{easybgm} fits the model with
#'   \code{bgms} instead.
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
#' @param not_cont A binary vector of length p, required whenever
#'   \code{type = "mixed"}. Each element indicates whether the corresponding
#'   variable is not continuous (\code{1} = not continuous/ordinal,
#'   \code{0} = continuous). When the model is fitted with \code{bgms}, this is
#'   used to construct a per-variable type vector internally; \code{BGGM} and
#'   \code{BDgraph} use it directly. It is not needed when \code{type} is
#'   already given as a per-variable vector.
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
#'     association estimates. Note that the scale differs between fitting
#'     packages: \code{BGGM} and \code{BDgraph} report partial correlations,
#'     whereas \code{bgms} reports the pairwise association parameter (the
#'     coupling that enters each conditional distribution). For \code{bgms}
#'     fits of continuous or mixed data the partial correlations are returned
#'     separately in \code{partial_correlations} (see below). Edge weights from
#'     different packages are therefore not directly comparable.
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
#'     interval for each inclusion Bayes factor, derived from the Monte Carlo
#'     standard error of the Rao-Blackwellized inclusion probability. Entries
#'     are \code{NA} where the interval is not defined, which happens when the
#'     posterior inclusion probability is numerically 0 or 1 and the Bayes
#'     factor is therefore 0 or infinite. This is routine for decisive edges
#'     rather than a sign of a problem.
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
#' Which values of \code{type} are legal depends on the package that fits the
#' model. \code{BGGM} and \code{BDgraph} fit continuous, mixed, ordinal and
#' binary data; the Blume-Capel model and per-variable \code{type} vectors are
#' fitted only by \code{bgms}. Because \code{bgms} 0.1.6.3 cannot fit
#' continuous data, the \code{bgms} column is split by version:
#'
#' \tabular{lcccc}{
#'   \strong{Data type} \tab \strong{bgms >= 0.2.0.0} \tab \strong{bgms 0.1.6.3} \tab \strong{BDgraph} \tab \strong{BGGM} \cr
#'   continuous           \tab Yes (default) \tab No             \tab Yes           \tab Yes \cr
#'   ordinal              \tab Yes (default) \tab Yes (default)  \tab Yes           \tab Yes \cr
#'   binary               \tab Yes (default) \tab Yes (default)  \tab Yes           \tab Yes \cr
#'   blume-capel          \tab Yes (default) \tab Yes (default)  \tab No            \tab No  \cr
#'   mixed                \tab Yes (default) \tab No             \tab Yes           \tab Yes \cr
#'   per-variable vector  \tab Yes           \tab No             \tab No            \tab No  \cr
#' }
#'
#' "Default" marks the package used when \code{package} is left unspecified.
#' An explicit \code{package} is honoured for every data type those packages
#' support. Requesting \code{BGGM} or \code{BDgraph} for the Blume-Capel model
#' or for a per-variable \code{type} vector produces a warning and falls back to
#' \code{bgms}.
#'
#' \strong{How \code{type} and \code{not_cont} relate}
#'
#' There are two ways to specify mixed-type data:
#' \enumerate{
#'   \item \strong{Using \code{type = "mixed"} with \code{not_cont}}: Set
#'     \code{type = "mixed"} and provide a binary vector \code{not_cont} of
#'     length p. When the model is fitted with \code{bgms}, this is translated
#'     internally to a per-variable type vector, where \code{not_cont == 1}
#'     maps to \code{"ordinal"} and \code{not_cont == 0} maps to
#'     \code{"continuous"}. \code{BGGM} and \code{BDgraph} use \code{not_cont}
#'     directly.
#'   \item \strong{Using a per-variable type vector} (\code{bgms} >= 0.2.0.0
#'     only): Directly specify the type of each variable, e.g.,
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
#'     Use \code{\link[bgms]{normal_prior}(scale)} (default
#'     \code{normal_prior(scale = 1)}), \code{\link[bgms]{cauchy_prior}(scale)},
#'     or \code{\link[bgms]{beta_prime_prior}(alpha, beta)}.
#'     For example, a cauchy prior with scale 1 would be specified with adding the
#'     argument \code{interaction_prior = cauchy_prior(1)} to the easybgm call.
#'     Note that the default is a Normal slab, which has considerably lighter
#'     tails than the Cauchy used by earlier versions of \code{bgms}, and so
#'     constrains weakly identified edges more tightly.
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
#'     \code{\link[bgms]{exponential_prior}(rate)} (default
#'     \code{exponential_prior(eta = 1)}) or
#'     \code{\link[bgms]{gamma_prior}(shape, rate)}.
#'   \item \code{precision_graph_prior}: How the graph prior is applied to the
#'     precision matrix (GGM and mixed MRF). Either \code{"hierarchical"} (the
#'     default), which tracks the normalizing constant, or \code{"joint"}.
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
#' For backwards compatibility, the legacy flat arguments below are still
#' accepted via \code{...} and are translated into the constructors above, so
#' they do not raise \code{bgms} deprecation warnings:
#' \itemize{
#'   \item \code{interaction_scale}, \code{pairwise_scale}: Scale of a Cauchy
#'     prior on the pairwise interactions, translated to
#'     \code{cauchy_prior(scale)}. Note that this is a Cauchy, whereas the
#'     current default is \code{normal_prior(scale = 1)}.
#'   \item \code{threshold_alpha}, \code{threshold_beta} (or \code{main_alpha},
#'     \code{main_beta}): Shape parameters translated to
#'     \code{beta_prime_prior(alpha, beta)}.
#'   \item \code{edge_prior}: Character string \code{"Bernoulli"},
#'     \code{"Beta-Bernoulli"}, or \code{"Stochastic-Block"}, combined with the
#'     flat hyperparameters \code{inclusion_probability},
#'     \code{beta_bernoulli_alpha}, \code{beta_bernoulli_beta},
#'     \code{dirichlet_alpha}, and \code{lambda}.
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
#' # --- Continuous data (fitted by bgms >= 0.2.0.0, otherwise by BGGM) ---
#' fit <- easybgm(data, type = "continuous",
#'                 iter = 100,   # for demonstration only; increase for real analyses
#'                 warmup = 300  # bgms defaults to 2000
#'                 )
#' summary(fit)
#'
#' \donttest{
#' # --- Mixed data using per-variable type vector (requires bgms >= 0.2.0.0) ---
#' if (utils::packageVersion("bgms") >= "0.2.0.0") {
#'   dat3 <- data[, 1:3]
#'   fit_vec <- easybgm(dat3,
#'                       type = c("ordinal", "ordinal", "continuous"),
#'                       iter = 100, warmup = 300)
#' }
#'
#' # --- Extract posterior samples and centrality ---
#' fit_full <- easybgm(data, type = "continuous",
#'                      iter = 100, warmup = 300,
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
  
  # Per-variable 'type' vectors require both bgms as the fitting package and
  # bgms >= 0.2.0.0. Resolve the fitting package first, then validate 'type'
  # against it: the legal values of 'type' depend on the engine, not on the
  # bgms version directly.
  bgms_supports_vector_type <- packageVersion("bgms") > "0.1.6.3"
  is_vector_type <- length(type) > 1

  if(is_vector_type && !bgms_supports_vector_type){
    stop("Specifying 'type' as a vector with one entry per variable requires ",
         "bgms version 0.2.0.0 or later (installed: ", packageVersion("bgms"), "). ",
         "Please either update bgms, or pass a single type for all variables.",
         call. = FALSE)
  }

  # --- Resolve the fitting package ---------------------------------------
  if(is.null(package)){
    if(bgms_supports_vector_type){
      # From bgms 0.2.0.0 onwards bgms handles every level of measurement.
      if(!is_vector_type && type %in% c("continuous", "mixed")){
        warning("Note that from bgms version 0.2.0.0 onwards, the default fit
              package changed for continuous and mixed data and is also
              bgms. If you prefer another package you can specify it with the
              package argument.",
                call. = FALSE)
      }
      package <- "package_bgms"
    } else {
      # bgms 0.1.6.3 only fits discrete data, so continuous/mixed fall to BGGM.
      package <- switch(type,
                        "continuous"  = "package_bggm",
                        "mixed"       = "package_bggm",
                        "ordinal"     = "package_bgms",
                        "binary"      = "package_bgms",
                        "blume-capel" = "package_bgms",
                        stop("Unrecognized 'type': '", type, "'.", call. = FALSE))
    }
  } else {
    package <- switch(package,
                      "BDgraph" = "package_bdgraph",
                      "BGGM"    = "package_bggm",
                      "bgms"    = "package_bgms",
                      stop("Unrecognized 'package': '", package, "'. ",
                           "Valid options are 'bgms', 'BGGM', and 'BDgraph'.",
                           call. = FALSE))

    # BGGM and BDgraph fit continuous, mixed, ordinal and binary data, so an
    # explicit package choice is honoured for those. Neither implements the
    # Blume-Capel model, and neither accepts a per-variable 'type' vector;
    # those two are bgms-only and override the choice. Warn rather than
    # switching silently.
    if(package != "package_bgms"){
      override_reason <- if(is_vector_type){
        "a per-variable 'type' vector"
      } else if(type == "blume-capel"){
        "type = 'blume-capel'"
      } else NULL

      if(!is.null(override_reason)){
        warning("BGGM and BDgraph cannot fit ", override_reason, "; the ",
                "'package' argument was overridden and bgms will be used ",
                "instead.",
                call. = FALSE)
        package <- "package_bgms"
      }
    }
  }

  # --- Validate 'type' against the resolved package ------------------------
  if(package == "package_bgms"){
    if(!bgms_supports_vector_type && !type %in% c("ordinal", "binary", "blume-capel")){
      stop("bgms version ", packageVersion("bgms"), " can only fit 'ordinal', ",
           "'binary', or 'blume-capel' data, not '", type, "'. Please either ",
           "update bgms to 0.2.0.0 or later, or fit '", type, "' data with ",
           "package = 'BGGM' or package = 'BDgraph'.",
           call. = FALSE)
    }
    # bgms has no 'mixed' type; it expresses mixedness through a per-variable
    # vector. Translate 'mixed' + not_cont into that vector.
    if(!is_vector_type && type == "mixed"){
      if(is.null(not_cont)){
        stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
             call. = FALSE)
      }
      if(length(not_cont) != ncol(data)){
        stop("'not_cont' has length ", length(not_cont), " but 'data' has ",
             ncol(data), " columns.", call. = FALSE)
      }
      type <- ifelse(not_cont == 1, "ordinal", "continuous")
      is_vector_type <- TRUE
    }
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
    if(is_vector_type && length(type) != ncol(data)) {
      stop("When 'type' is a vector, its length (", length(type), ") must equal ",
           "the number of columns in 'data' (", ncol(data), ").",
           call. = FALSE)
    }
    # Map "binary" to "ordinal" internally
    type[type == "binary"] <- "ordinal"
  } else {
    # Neither BGGM nor BDgraph implements the Blume-Capel model or accepts a
    # per-variable 'type' vector. Both were redirected to bgms above, so this
    # is a safety net.
    if(is_vector_type || type == "blume-capel"){
      stop(if(package == "package_bggm") "BGGM" else "BDgraph",
           " cannot fit ",
           if(is_vector_type) "a per-variable 'type' vector" else "'blume-capel' data",
           ". Please use package = 'bgms' for these data types.",
           call. = FALSE)
    }
    bggm_bdgraph_types <- c("continuous", "mixed", "ordinal", "binary")
    if(!type %in% bggm_bdgraph_types){
      stop("Unrecognized 'type': '", type, "'. ",
           if(package == "package_bggm") "BGGM" else "BDgraph",
           " fits: ", paste(bggm_bdgraph_types, collapse = ", "), ".",
           call. = FALSE)
    }
    if(type == "mixed" && is.null(not_cont)){
      stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
           call. = FALSE)
    }
    if(type == "mixed" && length(not_cont) != ncol(data)){
      stop("'not_cont' has length ", length(not_cont), " but 'data' has ",
           ncol(data), " columns.", call. = FALSE)
    }
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
