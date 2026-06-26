#' @title Compare networks across groups using Bayesian inference
#'
#' @description Easy comparison of networks using Bayesian inference to extract
#'   differences in conditional (in)dependence relations across groups.
#'
#' @name easybgm_compare
#'
#' @param data The data can be provided in two formats:
#'
#'   \strong{1. A list of two dataframes} (two-group comparison): Each list element
#'   is an n x p matrix or dataframe for one group. The variables (columns) must
#'   be the same across both dataframes. This format is supported by both
#'   \code{bgms} and \code{BGGM}.
#'
#'   \strong{2. A single matrix or dataframe} (multi-group comparison): An n x p
#'   matrix containing responses from all groups combined. Requires the
#'   \code{group_indicator} argument to specify which rows belong to which
#'   group. This format supports two or more groups and is only available with
#'   the \code{bgms} package.
#'
#' @param type Specifies the data type. Can be used in two ways:
#'
#'   \strong{1. A single string}:
#'   \itemize{
#'     \item \code{"continuous"}: For continuous data. Default package: BGGM.
#'     \item \code{"ordinal"}: For ordinal (Likert-type) data. Default package:
#'       bgms.
#'     \item \code{"binary"}: For binary (0/1) data. Default package: bgms.
#'     \item \code{"blume-capel"}: For Blume-Capel ordinal data. Requires
#'       \code{baseline_category}. Default package: bgms.
#'     \item \code{"mixed"}: For mixed data, requires \code{not_cont}. Default 
#'        package: BGGM.
#'   }
#'
#'   \strong{2. A character vector of length p} (per-variable specification,
#'   bgms only): Each element specifies the type of the corresponding column.
#'   Valid values are \code{"ordinal"}, \code{"blume-capel"}, and
#'   \code{"binary"}.
#'
#'   For example: \code{type = c("ordinal", "ordinal", "blume-capel")}.
#'
#' @param package The R-package used for fitting the comparison model. 
#'   If not specified, \code{bgms} is used as the default for ordinal, binary,
#'   and blume-capel data and \code{BGGM} as default for continous and mixed data. 
#'
#' @param not_cont A binary vector of length p, required when
#'   \code{type = "mixed"}. Each element indicates whether the corresponding
#'   variable is not continuous (\code{1} = not continuous/ordinal,
#'   \code{0} = continuous).
#'
#' @param group_indicator An integer vector of length n specifying group
#'   membership for each row in \code{data}. Required when \code{data} is a
#'   single matrix/dataframe (multi-group comparison). Supports two or more
#'   groups (e.g., \code{rep(c(1, 2, 3), each = 50)} for three groups of 50
#'   observations each). Only available with the \code{bgms} package.
#'
#' @param iter Number of iterations for the sampler. The default is 1e4.
#'   The recommended number of iterations depends on the data and model
#'   complexity. Check convergence diagnostics in the output.
#'
#' @param save Logical. Should the posterior samples be obtained
#'   (default = \code{TRUE})? If \code{TRUE}, the output includes a
#'   \code{samples_posterior} matrix with posterior samples for each difference
#'   parameter.
#'
#' @param progress Logical. Should a progress bar be shown
#'   (default = \code{TRUE})?
#'
#' @param ... Additional arguments passed to the fitting functions of the
#'   underlying packages (e.g., prior specifications). Consult the documentation
#'   of \code{bgms} and \code{BGGM} for the specific options available.
#'
#' @return An object of class \code{easybgm_compare} with the following
#'   elements:
#'
#' \strong{Always returned:}
#' \itemize{
#'   \item \code{parameters}: A p x p matrix of posterior mean differences in
#'     partial associations across groups.
#'   \item \code{inc_probs}: A p x p matrix of posterior inclusion probabilities
#'     for group differences (i.e., the probability that an edge differs between
#'     groups).
#'   \item \code{inc_BF}: A p x p matrix of inclusion Bayes factors for group
#'     differences.
#'   \item \code{structure}: A p x p adjacency matrix of the median probability
#'     model for differences (edges with posterior inclusion probability > 0.5).
#'   \item \code{model}: A string indicating the data type used.
#' }
#'
#' \strong{Returned for bgms only:}
#' \itemize{
#'   \item \code{structure_probabilities}: Posterior probabilities of all
#'     visited graph structures.
#'   \item \code{graph_weights}: Number of times each structure was visited.
#'   \item \code{sample_graph}: Identifiers for each visited structure.
#'   \item \code{convergence_parameter}: The Gelman-Rubin (R-hat) convergence
#'     statistic for each difference parameter. Values close to 1 indicate
#'     good convergence.
#' }
#'
#' \strong{Returned when save = TRUE:}
#' \itemize{
#'   \item \code{samples_posterior}: A k x iter matrix of posterior samples for
#'     each difference parameter (k = p*(p-1)/2 edges).
#' }
#'
#' @details
#'
#' \strong{Data types and package support for group comparison}
#'
#' \tabular{lcc}{
#'   \strong{Data type}  \tab \strong{bgms} \tab \strong{BGGM} \cr
#'   ordinal              \tab Yes (default)  \tab No            \cr
#'   binary               \tab Yes (always)   \tab Yes            \cr
#'   blume-capel          \tab Yes (default)  \tab No            \cr
#'   continuous            \tab No             \tab Yes (default) \cr
#'   mixed                 \tab No             \tab Yes (default) \cr
#' }
#'
#'
#' \strong{Two-group vs. multi-group comparison}
#'
#' \itemize{
#'   \item \strong{Two-group}: Provide \code{data} as a list of two
#'     dataframes. Supported by both \code{bgms} and \code{BGGM}.
#'   \item \strong{Multi-group}: Provide \code{data} as a single
#'     matrix/dataframe and specify \code{group_indicator}. Supports two or more
#'     groups. Only available with \code{bgms}.
#' }
#'
#' \strong{Prior specification}
#'
#' Users may wish to adjust priors via the \code{...} argument. We summarize
#' the bgms options here and refer to \code{\link[bgms]{bgmCompare}} and
#' \code{\link[BGGM]{explore}} for full details.
#'
#' \emph{bgms} (>= 0.2.0.0): the preferred interface uses prior-constructor
#' objects:
#' \itemize{
#'   \item \code{interaction_prior}: Prior on the baseline pairwise
#'     interactions. Use \code{\link[bgms]{cauchy_prior}(scale)} (default
#'     \code{cauchy_prior(scale = 1)}), \code{\link[bgms]{normal_prior}(scale)},
#'     or \code{\link[bgms]{beta_prime_prior}(alpha, beta)}.
#'   \item \code{threshold_prior}: Prior on threshold parameters. Default
#'     \code{beta_prime_prior(0.5, 0.5)}.
#'   \item \code{difference_prior}: Indicator prior on group differences.
#'     Use \code{\link[bgms]{bernoulli_prior}(inclusion_probability)} (default
#'     \code{bernoulli_prior(0.5)}) or
#'     \code{\link[bgms]{beta_bernoulli_prior}(alpha, beta)}.
#'   \item \code{difference_scale}: Scale of the prior on the magnitude of
#'     pairwise differences. Default 1.
#'   \item \code{difference_selection}: Logical, whether to perform Bayesian
#'     selection on group differences. Default \code{TRUE}.
#' }
#'
#' For backwards compatibility, the legacy flat arguments below are still
#' accepted via \code{...} and are translated internally:
#' \itemize{
#'   \item \code{pairwise_scale}: Scale of the Cauchy interaction prior.
#'   \item \code{difference_prior}: Character string \code{"Bernoulli"}
#'     (default) or \code{"Beta-Bernoulli"}.
#'   \item \code{difference_probability}: Prior probability of a difference for
#'     the Bernoulli prior. Default 0.5.
#'   \item \code{beta_bernoulli_alpha}, \code{beta_bernoulli_beta}: Shape
#'     parameters of the Beta-Bernoulli difference prior. Both default to 1.
#'   \item \code{threshold_alpha}, \code{threshold_beta} (or
#'     \code{main_alpha}, \code{main_beta}): Beta-prime parameters of the
#'     threshold prior. Both default to 0.5.
#' }
#'
#' We always encourage researchers to conduct prior sensitivity checks.
#'
#' @export
#' @import bgms
#' @importFrom BGGM explore select
#' @importFrom utils packageVersion
#'
#' @examples
#'
#' \donttest{
#' library(easybgm)
#' library(bgms)
#'
#' data <- na.omit(ADHD)
#'
#' # --- Two-group comparison (list input) ---
#' group1 <- data[1:10, 1:3]
#' group2 <- data[11:20, 1:3]
#'
#' fit <- easybgm_compare(list(group1, group2),
#'                 type = "binary", save = TRUE,
#'                 iter = 50  # for demonstration only
#'                 )
#' summary(fit)
#'
#' # --- Multi-group comparison (single dataframe + group_indicator) ---
#' fit_multi <- easybgm_compare(data[1:200, 1:5],
#'                 group_indicator = rep(c(1, 2, 3, 4), each = 50),
#'                 type = "binary", save = TRUE,
#'                 iter = 100  # for demonstration only
#'                 )
#' summary(fit_multi)
#' }


easybgm_compare <- function(data, 
                            type, 
                            package = NULL, 
                            not_cont = NULL, 
                            group_indicator = NULL, 
                            iter = 1e3,
                            save = TRUE, 
                            progress = TRUE,
                            ...){
  
  if(!is.list(data) && is.null(group_indicator)){
    stop("Your data can't be read. There are two options of providing your data: 1) Provide two datasets in a list containing only the two datasets, or for ordinal data with the bgms pacakge > 0.1.6. 2) provide the data as a matrix or data.frame together with specifying the 'group_indicator' argument, which then also allows for multi-group comparison.",
         call. = FALSE)
  }
  
  if(packageVersion("bgms") > "0.1.6.3"){
    # --- Handle vector type (per-variable specification) ---
    is_vector_type <- length(type) > 1
    
    if(is_vector_type) {
      valid_types <- c("ordinal", "blume-capel", "binary")
      invalid <- type[!type %in% valid_types]
      if(length(invalid) > 0) {
        warning("The following variable type(s) are not recognized for group comparison: ",
                paste0("'", unique(invalid), "'", collapse = ", "), ". ",
                "Valid types are: ", paste(valid_types, collapse = ", "), ". ",
                "Note: 'continuous' is not supported for group comparison via bgms.",
                call. = FALSE)
        stop("Invalid variable types detected. See the warning message for more details.",
             call. = FALSE)
      }
      ncols <- if(is.list(data) && !is.data.frame(data)) ncol(data[[1]]) else ncol(data)
      if(length(type) != ncols) {
        stop("When 'type' is a vector, its length (", length(type), ") must equal ",
             "the number of columns in 'data' (", ncols, ").",
             call. = FALSE)
      }
      type[type == "binary"] <- "ordinal"
      package <- "package_bgms_compare"
    }
    
    
    if(!is_vector_type && length(type) == 1 && type == "mixed" && is.null(not_cont)){
      stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
           call. = FALSE)
    }
  } else if(packageVersion("bgms") < "0.2.0.0"){
    if(type == "mixed" & is.null(not_cont)){
      stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
           call. = FALSE)
    }
  }
  
  
  dots <- list(...)
  has_reference <- "reference_category" %in% names(dots)
  has_baseline  <- "baseline_category" %in% names(dots)
  
  if (any(type == "blume-capel") && !(has_reference || has_baseline)) {
    stop("For the Blume-Capel model, a reference category needs to be specified.
         Use the baseline_category argument to specify the reference category.",
         call. = FALSE)
  }
  
  # Set default values for fitting if package is unspecified
  if(type %in% c("blume-capel", "ordinal")){
    package <- "package_bgms_compare"
  } else if(length(type) > 1){
    package <- "package_bgms_compare"
  } else if(is.null(package)){
    if(type == "ordinal") package <- "package_bgms_compare"
    if(type == "binary") package <- "package_bgms_compare"
    if(type == "blume-capel") package <- "package_bgms_compare"
    if(type == "continuous") package <- "package_bggm_compare"
    if(type == "mixed") package <- "package_bggm_compare"
  } else {
    if(package == "BGGM") package <- "package_bggm_compare"
    if(package == "bgms") package <- "package_bgms_compare"
  }
  
  
  if((is.data.frame(data) || is.matrix(data)) && package == "package_bggm_compare"){
    stop("Your data can't be read. For continuous data fit with BGGM, 
         you can only provide two datasets in a list.",
         call. = FALSE)
  }
  
  fit <- list()
  class(fit) <- c(package, "easybgm")
  
  # Fit the model
  tryCatch(
    {fit <- bgm_fit(fit, data = data, type = type, not_cont = not_cont, 
                    group_indicator = group_indicator, 
                    iter = iter,
                    save = save, progress = progress, ...)
    },
    error = function(e){
      # If an error occurs, stop running the code
      stop(paste("Error meassage: ", e$message, "Please consult the original message for more information.") )
    })
  
  # Extract the results
  res <- bgm_extract(fit, type = type,
                     save = save, not_cont = not_cont,
                     group_indicator = group_indicator, 
                     data = data, 
                     ...)
  
  # Output results
  class(res) <- c(package, "easybgm_compare", "easybgm")
  return(res)
}
