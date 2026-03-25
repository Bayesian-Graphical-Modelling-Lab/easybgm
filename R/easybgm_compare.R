#' @title Compare networks across groups using Bayesian inference
#'
#' @description Easy comparison of networks using Bayesian inference to extract
#'   differences in conditional (in)dependence relations across groups.
#'
#' @name easybgm_compare
#'
#' @param data The data can be provided in two formats:
#'
#'   \strong{1. A list of two dataframes} (two-group comparison): Each element
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
#' @param type Specifies the data type. Currently supported types for group
#'   comparison:
#'   \itemize{
#'     \item \code{"ordinal"}: For ordinal (Likert-type) data. Default package:
#'       bgms.
#'     \item \code{"binary"}: For binary (0/1) data. Always uses bgms.
#'     \item \code{"blume-capel"}: For Blume-Capel ordinal data. Requires
#'       \code{baseline_category}. Always uses bgms.
#'     \item \code{"continuous"}: Only supported with
#'       \code{package = "BGGM"} (must be set explicitly). Not yet supported via
#'       bgms.
#'     \item \code{"mixed"}: Only supported with
#'       \code{package = "BGGM"} (must be set explicitly). Requires
#'       \code{not_cont}. Not yet supported via bgms for group comparison.
#'   }
#'
#' @param package The R-package used for fitting the comparison model. Optional;
#'   if not specified, \code{bgms} is used as the default for ordinal, binary,
#'   and blume-capel data. Supported options:
#'   \itemize{
#'     \item \code{"bgms"} (Default): Supports ordinal, binary, and blume-capel
#'       data. Supports both two-group (list input) and multi-group (single
#'       dataframe + \code{group_indicator}) comparisons.
#'     \item \code{"BGGM"}: Supports continuous and mixed data. Only supports
#'       two-group comparison (list input). Binary data is always routed to
#'       bgms regardless.
#'   }
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
#' @param iter Number of iterations for the sampler. The default depends on the
#'   package:
#'   \itemize{
#'     \item \code{bgms}: 1e4 (10,000 iterations)
#'     \item \code{BGGM}: 1e4 (10,000 iterations)
#'   }
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
#'   binary               \tab Yes (always)   \tab No            \cr
#'   blume-capel          \tab Yes (default)  \tab No            \cr
#'   continuous            \tab No             \tab Yes (explicit)\cr
#'   mixed                 \tab No             \tab Yes (explicit)\cr
#' }
#'
#' For continuous or mixed data comparison, you must explicitly set
#' \code{package = "BGGM"}. If no package is specified for these types, the
#' function will return an informative error.
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
#' Users may wish to adjust priors via the \code{...} argument. Consult the
#' documentation of \code{\link[bgms]{bgm}} and \code{\link[BGGM]{explore}} for
#' specific prior options.
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
  if(type %in% c("continuous", "mixed") && (is.null(package) || package == "bgms")){
    stop("Group comparison via bgms currently only supports ordinal and binary data types.
         For continuous or mixed data comparison, explicitly set package = 'BGGM'.",
         call. = FALSE)
  }
  if(type == "mixed" & is.null(not_cont)){
    stop("Please provide a binary vector of length p specifying the not continuous variables
         (1 = not continuous, 0 = continuous).",
         call. = FALSE)
  }


  dots <- list(...)
  has_reference <- "reference_category" %in% names(dots)
  has_baseline  <- "baseline_category" %in% names(dots)

  # Example: If type == "blume-capel", then at least one must be present
  if (type == "blume-capel" && !(has_reference || has_baseline)) {
    stop("For the Blume-Capel model, a reference category needs to be specified. 
         If type is 'blume-capel' it specifies the reference category in the Blume-Capel model.
         Should be an integer within the range of integer scores observed for the
         'blume-capel' variable. Can be a single number specifying the reference
         category for all Blume-Capel variables at once, or a vector of length
         p where the i-th element contains the reference category for
         variable i if it is Blume-Capel, and bgm ignores its elements for
         other variable types. The value of the reference category is also recoded
         when bgm recodes the corresponding observations. Only required if there is at
         least one variable of type ``blume-capel''.
         For bgms version smaller than 0.1.6, use the reference_category argument. 
         For all package versions including and older than 0.1.6., the baseline_category argument.",
         call. = FALSE)
  }
  
  # Set default values for fitting if package is unspecified
  if(is.null(package)){
    if(type == "ordinal") package <- "package_bgms_compare"
    if(type == "binary") package <- "package_bgms_compare"
    if(type == "blume-capel") package <- "package_bgms_compare"
  } else {
    if(package == "BGGM") package <- "package_bggm_compare"
    if(package == "bgms") package <- "package_bgms_compare"
    if(type == "binary") package <- "package_bgms_compare"
  }

  # change default number of iterations for BGGM fits
  if (iter == 1e3 && package == "package_bgms_compare"){
    iter <- 1e4
  }
  
  if(is.data.frame(data) &&  package == "package_bggm_compare"){
    stop("Your data can't be read. For continuous data fit with BGGM, you can only provide two datasets in a list.",
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
