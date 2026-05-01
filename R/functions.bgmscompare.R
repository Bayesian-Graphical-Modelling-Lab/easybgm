# --------------------------------------------------------------------------------------------------
# 1. Fitting function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_fit.package_bgms_compare <- function(fit, type, data, group_indicator, iter, save,
                                         not_cont, progress, ...){

  # Map binary to ordinal (vector-safe)
  type[type == "binary"] <- "ordinal"

  # Determine model label for display
  if(length(type) > 1) {
    model_label <- if(length(unique(type)) == 1) unique(type) else "mixed"
  } else {
    model_label <- type
  }

  if(is.list(data) && !is.data.frame(data)){
    group_indicator <- NULL
  }

  bgmcompare_args <- translate_bgmcompare_prior_args(list(...))
  if(is.null(group_indicator)){
    bgms_fit <- do.call(
      bgmCompare, c(list(x = data[[1]], y = data[[2]], iter = iter,
                         variable_type = type,
                         display_progress = progress),
                    bgmcompare_args)
    )

    fit$model <- model_label
    fit$packagefit <- bgms_fit
    if(is.null(colnames(data[[1]]))){
      fit$var_names <- paste0("V", 1:ncol(data[[1]]))
    } else {
      fit$var_names <- colnames(data[[1]])
    }
  } else {
    bgms_fit <- do.call(
      bgmCompare, c(list(x = data, group_indicator = group_indicator,
                         iter = iter,
                         variable_type = type,
                         display_progress = progress),
                    bgmcompare_args)
    )

    fit$model <- model_label
    fit$packagefit <- bgms_fit
    if(is.null(colnames(data))){
      fit$var_names <- paste0("V", 1:ncol(data))
    } else {
      fit$var_names <- colnames(data)
    }
  }

  class(fit) <- c("package_bgms_compare", "easybgm")
  return(fit)
}




# --------------------------------------------------------------------------------------------------
# 2. Extracting results function
# --------------------------------------------------------------------------------------------------
#' @export
bgm_extract.package_bgms_compare <- function(fit, type, save, group_indicator,
                                             not_cont, data, centrality, ...){
  if(any(class(fit) == "easybgm")){
    varnames <- fit$var_names
    fit <- fit$packagefit
  } else if (any(class(fit) == "bgmCompare")){
    varnames <- extract_arguments(fit)$data_columnnames
    if(is.null(varnames)){
      varnames <- paste0("V", 1:extract_arguments(fit)$num_variables)
    }
  }

  ######--------------------
  ## Two group estimation
  ######--------------------
  if(is.null(group_indicator)){

    args <- extract_arguments(fit)
    args$save <- TRUE
    dots <- list(...)
    if (args$difference_prior[1] == "Bernoulli") {
      if("difference_probability" %in% dots){
        edge.prior <- args$difference_probability
        args$inclusion_probability_difference <- edge.prior
      } else {
        edge.prior <- 0.5
        args$inclusion_probability_difference <- edge.prior
      }
    } else { # if BB or SBM
      edge.prior <- args$difference_selection_alpha /
        (args$difference_selection_alpha + args$difference_selection_beta)

      # otherwise it saves the wrong values (could be done more elegantly)
      args$inclusion_probability_difference <- edge.prior
    }

    bgms_res <- list()

    p <- args$num_variables
    pars <- summary(fit)$pairwise_diff$mean
    bgms_res$parameters <- vector2matrix(pars, p = p)
    colnames(bgms_res$parameters) <- varnames
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters),
                                 nrow = nrow(bgms_res$parameters))
    inc_prob_mat <- extract_posterior_inclusion_probabilities(fit)
    diag(inc_prob_mat) <- 0
    bgms_res$inc_probs <- inc_prob_mat
    bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1 - edge.prior))
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)

    #Obtain structure information
    bgms_res$group_estimates <- extract_group_params(fit)$pairwise_effects_groups
    bgms_res$parameters_g1 <- vector2matrix(extract_group_params(fit)$pairwise_effects_groups[, 1], p = p)
    bgms_res$parameters_g2 <- vector2matrix(extract_group_params(fit)$pairwise_effects_groups[, 2], p = p)

    structures <- apply(extract_indicators(fit), 1, paste0, collapse="")
    table_structures <- as.data.frame(table(structures))
    bgms_res$structure_probabilities <- table_structures[,2]/nrow(extract_indicators(fit))
    bgms_res$graph_weights <- table_structures[,2]
    bgms_res$sample_graph <- as.character(table_structures[, 1])
    bgms_res$samples_posterior <- extract_pairwise_interactions(fit)

    bgms_res$convergence_parameter <- extract_rhat(fit)$pairwise_differences
  }

  ######--------------------
  ## Multi-group estimation
  ######--------------------
  if(!is.null(group_indicator)){

    args <- extract_arguments(fit)
    args$save <- TRUE
    dots <- list(...)
    if (args$difference_prior[1] == "Bernoulli") {
      if("difference_probability" %in% dots){
        edge.prior <-  args$difference_probability
        args$inclusion_probability_difference <- edge.prior
      } else {
        edge.prior <- 0.5
        args$inclusion_probability_difference <- edge.prior
      }
    } else { # if BB or SBM
      edge.prior <- args$difference_selection_alpha /
        (args$difference_selection_alpha + args$difference_selection_beta)

      # otherwise it saves the wrong values (could be done more elegantly)
      args$inclusion_probability_difference <- edge.prior
    }

    bgms_res <- list()

    p <- args$num_variables
    # Compute average group difference
    diffs <- summary(fit)$pairwise_diff
    edge_labels <- sub(" .*", "", diffs$parameter)
    edges <- unique(edge_labels)
    average_difference <- numeric(length(edges))
    for (i in seq_along(edges)) {
      edge_name <- edges[i]
      idx <- edge_labels == edge_name
      vals <- diffs$mean[idx]
      average_difference[i] <- mean(vals, na.rm = TRUE)
    }
    bgms_res$parameters <- vector2matrix(average_difference, p = p)
    colnames(bgms_res$parameters) <- varnames
    bgms_res$overall_estimate <- vector2matrix(extract_group_params(fit)$pairwise_effects_groups[, 1], p = p)
    bgms_res$structure <- matrix(1, ncol = ncol(bgms_res$parameters),
                                 nrow = nrow(bgms_res$parameters))
    inc_prob_mat <- extract_posterior_inclusion_probabilities(fit)
    diag(inc_prob_mat) <- 0
    bgms_res$inc_probs <- inc_prob_mat
    bgms_res$inc_BF <- (bgms_res$inc_probs/(1-bgms_res$inc_probs))/(edge.prior /(1 - edge.prior))
    bgms_res$structure <- 1*(bgms_res$inc_probs > 0.5)

    structures <- apply(extract_indicators(fit), 1, paste0, collapse="")
    table_structures <- as.data.frame(table(structures))
    bgms_res$structure_probabilities <- table_structures[,2]/nrow(extract_indicators(fit))
    bgms_res$graph_weights <- table_structures[,2]
    bgms_res$sample_graph <- as.character(table_structures[, 1])
    bgms_res$samples_posterior <- extract_pairwise_interactions(fit)
    bgms_res$convergence_parameter <- extract_rhat(fit)$pairwise_baseline
    bgms_res$group_estimates <- extract_group_params(fit)$pairwise_effects_groups
    bgms_res$multi_group <- "Multi-group"

  }

  # Adapt column names of output
  colnames(bgms_res$inc_probs) <- colnames(bgms_res$parameters)
  colnames(bgms_res$inc_BF) <- colnames(bgms_res$parameters)

  bgms_res$model <- if(length(type) > 1) {
    if(length(unique(type)) == 1) unique(type) else "mixed"
  } else {
    type
  }
  bgms_res$fit_arguments <- args
  bgms_res$edge.prior <- edge.prior # otherwise it stores a whole matrix

  output <- bgms_res
  class(output) <- c("package_bgms_compare", "easybgm_compare", "easybgm")
  return(output)
}