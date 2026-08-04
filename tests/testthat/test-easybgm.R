
### how do i vary the versions of bgms with easybgm 
##### CROSS-SECTIONAL
###-------------
### Estimation checks 
###-------------

test_that("easybgm returns expected structure across valid type–package combos", {
  set.seed(123)
  
  # Subsample small data to stay fast on CRAN
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:20, 1:5]
  p <- ncol(dat)
  itr <- 10
  
  if(packageVersion("bgms") > "0.1.6.3"){
  # Test only core combinations
  combos <- list(
    ### BGGM
    list(type = "continuous", pkg = "BGGM", sv = F, cnt = F),
    list(type = "continuous", pkg = "BGGM", sv = T, cnt = T),
    list(type = "mixed", pkg = "BGGM", sv = T, cnt = T),
    ### BDGRAPH
    list(type = "mixed",      pkg = "BDgraph", sv = F, cnt = F),
    list(type = "continuous",  pkg = "BDgraph", sv = F, cnt = F),
    ### bgms
    list(type = "binary",     pkg = "bgms", sv = F, cnt = F),
    list(type = "binary",     pkg = "bgms", sv = T, cnt = T),
    list(type = "binary",     pkg = "bgms", sv = F, cnt = T),
    list(type = "blume-capel", pkg = "bgms", sv = T, cnt = T),
    list(type = "binary", pkg = "bgms", sv = T, cnt = T, sbm = "Stochastic-Block"),
    list(type = "continuous", pkg = "bgms", sv = F, cnt = F),
    list(type = "continuous", pkg = "bgms", sv = T, cnt = T),
    list(type = "mixed",      pkg = "bgms", sv = T, cnt = T),
    list(type = c("ordinal", "ordinal", "continuous", "continuous", "ordinal"), pkg = "bgms", sv = F, cnt = F)
  )
  } else if(packageVersion("bgms") < "0.2.0.0"){
    # Test only core combinations
    combos <- list(
      ### BGGM
      list(type = "continuous", pkg = "BGGM", sv = F, cnt = F),
      list(type = "continuous", pkg = "BGGM", sv = T, cnt = T),
      list(type = "mixed", pkg = "BGGM", sv = T, cnt = T),
      ### BDGRAPH
      list(type = "mixed",      pkg = "BDgraph", sv = F, cnt = F),
      list(type = "continuous",  pkg = "BDgraph", sv = F, cnt = F),
      ### bgms
      list(type = "binary",     pkg = "bgms", sv = F, cnt = F),
      list(type = "binary",     pkg = "bgms", sv = T, cnt = T),
      list(type = "binary",     pkg = "bgms", sv = F, cnt = T),
      list(type = "blume-capel", pkg = "bgms", sv = T, cnt = T),
      list(type = "binary", pkg = "bgms", sv = T, cnt = T, sbm = "Stochastic-Block")
    )
  }
  
  for (cmb in combos) {
    t <- cmb$type
    pkg <- cmb$pkg
    sv <- cmb$sv
    cnt <- cmb$cnt
    if(!is.null(cmb$sbm)) {sbm <- cmb$sbm}
    
    not_cont <- if (length(t) == 1 && t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL

    # bgms defaults to warmup = 2000, which dominates the runtime at these tiny
    # iteration counts. 300 is the shortest warmup bgms does not warn about.
    # BGGM and BDgraph have no warmup argument.
    extra <- if (identical(pkg, "bgms")) list(warmup = 300) else list()

    base_args <- list(
      data       = dat,
      type       = t,
      package    = pkg,
      iter       = itr,          # tiny for speed
      save       = sv,
      centrality = cnt,
      progress   = FALSE
    )

    if(length(t) == 1 && t == "blume-capel"){
      suppressWarnings({
        res <- do.call(easybgm, c(base_args, extra,
                                  list(not_cont = not_cont,
                                       baseline_category = 2)))
      })} else if(!is.null(cmb$sbm)){
        suppressWarnings({
          res <- do.call(easybgm, c(base_args, extra,
                                    list(edge_prior = sbm)))
        })
      } else {
        suppressWarnings({
          res <- do.call(easybgm, c(base_args, extra,
                                    list(not_cont = not_cont)))
        })
      }
    
    # --- class check ---
    expect_true(inherits(res, c("easybgm")))
    expect_true(any(grepl("package_", class(res))))  # backend tag present
    
    # --- field presence check ---
    expect_true(all(c("parameters", "inc_probs", "inc_BF", "structure", "model") %in% names(res)))
    
    # --- dimensions check ---
    expect_equal(dim(res$parameters), c(p, p))
    expect_equal(dim(res$inc_probs),  c(p, p))
    expect_equal(dim(res$inc_BF),     c(p, p))
    expect_equal(dim(res$structure),  c(p, p))
    
    # --- sanity check ---
    expect_false(all(is.na(res$parameters)))
    expect_false(all(is.na(res$inc_probs))) 
    
    expect_no_error(summary(res))
    
    if(sv == TRUE && pkg == "BGGM") {
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(itr, k))
      expect_equal(dim(res$centrality),  c(itr, p))
    } 
    if(cnt == TRUE && pkg == "bgms"){
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(4*itr, k))
      expect_equal(dim(res$centrality),  c(4*itr, p))
    }
    if(!is.null(cmb$sbm)){
      expect_equal(length(res$sbm), 4)
    }
    print(paste0("Finished easybgm: Package: ", cmb$pkg, "; Type: ", cmb$type, "; Centrality: ", cmb$cnt))
    
  }
})

###-------------
### Plotting functions test
###-------------

test_that("plotting functions work across valid type–package combos", {
  set.seed(123)
  
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:20, 1:5]
  p   <- ncol(dat)
  
  combos <- list(
    list(type = "continuous", pkg = "BGGM"),
   # list(type = "mixed",      pkg = "BDgraph"),
    list(type = "binary",     pkg = "bgms")
  )
  
  for (cmb in combos) {
    t   <- cmb$type
    pkg <- cmb$pkg
    not_cont <- if (t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL
    
    
    if(pkg == "BDgraph") {
      suppressMessages({
        res <- easybgm(
          data       = dat,
          type       = t,
          package    = pkg,
          iter       = 10,
          save       = FALSE,
          centrality = TRUE,
          progress   = FALSE,
          not_cont   = not_cont
        )
      }) 
    } else {
      # bgms defaults to warmup = 2000; BGGM has no warmup argument
      extra <- if (identical(pkg, "bgms")) list(warmup = 300) else list()
      suppressMessages({
        res <- do.call(easybgm, c(
          list(
            data       = dat,
            type       = t,
            package    = pkg,
            iter       = 10,
            save       = TRUE,
            centrality = TRUE,
            progress   = FALSE,
            not_cont   = not_cont
          ), extra))
      })
    }
    
    # --- edge evidence ---
    g1 <- invisible(plot_edgeevidence(res))
    expect_true(inherits(g1, c("ggplot", "qgraph")))
    
    # --- network ---
    g2 <- invisible(plot_network(res))
    expect_true(inherits(g2, c("ggplot", "qgraph")))
    
    # --- structure plots (skip for BGGM) ---
    if (pkg != "BGGM") {
      g3 <- invisible(plot_structure_probabilities(res))
      expect_s3_class(g3, "ggplot")
      
      g4 <- invisible(plot_complexity_probabilities(res))
      expect_s3_class(g4, "ggplot")
      
      g5 <- invisible(plot_structure(res))
      expect_true(inherits(g5, c("ggplot", "qgraph")))
    }
    
    # --- posterior parameter HDI ---
    if(pkg != "BDgraph"){
      g6 <-    suppressWarnings({invisible(plot_parameterHDI(res))})
      expect_s3_class(g6, "ggplot")
      
      # --- centrality ---
      g7 <- invisible(plot_centrality(res))
      expect_s3_class(g7, "ggplot")
    }
  }
})

# # TEst only possible to include post 0.2.0.0 version
# test_that("easybgm defaults to bgms for all data types", {
#   data("Wenchuan", package = "bgms")
#   dat <- na.omit(Wenchuan)[1:20, 1:5]
#   
#   suppressWarnings({
#     res <- easybgm(dat, type = "continuous", iter = 10, progress = FALSE)
#   })
#   expect_true("package_bgms" %in% class(res))
#   
#   suppressWarnings({
#     res2 <- easybgm(dat, type = "ordinal", iter = 10, progress = FALSE)
#   })
#   expect_true("package_bgms" %in% class(res2))
# })

test_that("easybgm errors for BDgraph continuous with missing data", {
  data("Wenchuan", package = "bgms")
  dat_with_na <- Wenchuan[1:20, 1:5]  # Wenchuan has NAs
  
  expect_error(
    easybgm(dat_with_na, type = "continuous", package = "BDgraph",
            iter = 10, progress = FALSE),
    "missing values"
  )
})


##### NETWORK COMPARISON

test_that("easybgm_compare errors for continuous/mixed without BGGM", {
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:20, 1:5]
  group_dat <- list(dat[1:10, ], dat[11:20, ])
  expect_error(
    suppressWarnings(
      easybgm_compare(group_dat, type = "continuous", package = "bgms")
    ),
    "Invalid variable types detected"
  )
  expect_error(
    suppressWarnings(
      easybgm_compare(group_dat, type = "mixed", package = "bgms")
    ),
    "Invalid variable types detected"
  )
})


test_that("easybgm_compare accepts a per-variable type vector", {
  skip_if(packageVersion("bgms") <= "0.1.6.3")
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:30, 1:3]
  grp <- rep(1:2, length.out = nrow(dat))
  fit <- suppressWarnings(
    easybgm_compare(dat, type = rep("ordinal", 3), group_indicator = grp,
                    iter = 50, warmup = 300, progress = FALSE)
  )
  expect_s3_class(fit, "package_bgms_compare")

  # a vector whose length does not match the number of columns is rejected
  expect_error(
    suppressWarnings(
      easybgm_compare(dat, type = rep("ordinal", 2), group_indicator = grp,
                      iter = 50, progress = FALSE)
    ),
    "must equal"
  )
})


test_that("easybgm_compare returns expected structure across valid type–package combos", {
  set.seed(123)
  
  # Subsample small data to stay fast on CRAN
  data("Wenchuan", package = "bgms")
  dat <- as.data.frame(na.omit(Wenchuan)[1:90, 1:5])
  p <- ncol(dat)
  itr <- 10
  
  # Test only core combinations
  combos <- list(
    ### BGGM
    list(type = "continuous", pkg = "BGGM", sv = F),
    list(type = "continuous", pkg = "BGGM", sv = T),
    list(type = "mixed", pkg = "BGGM", sv = T),
    ### bgms
    list(type = "binary",     pkg = "bgms", sv = F),
    list(type = "binary",     pkg = "bgms", sv = T),
    list(type = "binary",     pkg = "bgms", sv = T, multi_group = T)
  )
  
  for (cmb in combos) {
    t <- cmb$type
    pkg <- cmb$pkg
    sv <- cmb$sv
    
    # bgms defaults to warmup = 2000; BGGM has no warmup argument
    extra <- if (identical(pkg, "bgms")) list(warmup = 300) else list()

    if(!is.null(cmb$multi_group)){
      group <- rep(c(1, 2, 3), each = 30)

      suppressMessages({
        res <- do.call(easybgm_compare, c(
          list(
            data       = dat,
            type       = t,
            package    = pkg,
            iter       = itr,          # tiny for speed
            save       = sv,
            group_indicator = group,
            progress   = FALSE
          ), extra))
      })
    } else {
      group_dat <- list(dat[1:45, ], dat[46:90, ])
      not_cont <- if (t == "mixed") c(TRUE, TRUE, rep(FALSE, p - 2)) else NULL

      suppressWarnings({
        res <- do.call(easybgm_compare, c(
          list(
            data       = group_dat,
            type       = t,
            package    = pkg,
            iter       = itr,          # tiny for speed
            save       = sv,
            progress   = FALSE,
            not_cont   = not_cont
          ), extra))
      })
    }
    # --- class check ---
    expect_true(inherits(res, c("easybgm_compare")))
    expect_true(any(grepl("package_", class(res))))  # backend tag present
    
    # --- field presence check ---
    expect_true(all(c("parameters", "inc_probs", "inc_BF", "structure", "model") %in% names(res)))
    
    # --- dimensions check ---
    expect_equal(dim(res$parameters), c(p, p))
    expect_equal(dim(res$inc_probs),  c(p, p))
    expect_equal(dim(res$inc_BF),     c(p, p))
    expect_equal(dim(res$structure),  c(p, p))
    
    # --- sanity check ---
    expect_false(all(is.na(res$parameters)))
    expect_false(all(is.na(res$inc_probs)))
    
    if(sv == TRUE && pkg != "bgms") {
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(itr, k))
    }
    if(sv == TRUE && pkg == "bgms"){
      k <- p*(p-1)/2
      expect_equal(dim(res$samples_posterior), c(4*itr, k))
    }
    
    
    print(paste0("Finished easybgm_compare: Package: ", cmb$pkg, "; Type: ", cmb$type))
  }
})


###-------------
### Regression checks for bgms >= 0.2.0.0 extraction
###-------------

test_that("bgms centrality uses the lower-triangle edge ordering", {
  # bgms has stored pairwise interactions in lower-triangle order since at least
  # 0.1.6.3, so this holds on both supported bgms versions.
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:5]

  res <- suppressWarnings(
    easybgm(dat, type = "ordinal", iter = 50, warmup = 300,
            save = TRUE, centrality = TRUE, progress = FALSE)
  )
  p <- ncol(res$parameters)

  # bgms stores pairwise interactions in lower-triangle column order, the same
  # order res$parameters is filled from.
  expected <- t(apply(res$samples_posterior, 1, function(r)
    rowSums(abs(vector2matrix(r, p, bycolumn = FALSE)))))
  expect_equal(unname(res$centrality), unname(expected))

  # The upper-triangle fill (BGGM's order) genuinely differs, so the check above
  # would fail if the ordering regressed.
  wrong <- t(apply(res$samples_posterior, 1, function(r)
    rowSums(abs(vector2matrix(r, p, bycolumn = TRUE)))))
  expect_false(isTRUE(all.equal(unname(expected), unname(wrong))))
})

test_that("structure is the median probability model when save = FALSE", {
  # both supported bgms versions report inclusion probabilities the same way
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:4]

  res <- suppressWarnings(
    easybgm(dat, type = "ordinal", iter = 50, warmup = 300,
            save = FALSE, progress = FALSE)
  )
  expect_equal(unname(res$structure), unname(1 * (res$inc_probs > 0.5)))
  expect_true(all(diag(res$structure) == 0))
})

test_that("entry points accept raw bgms fit objects", {
  # Exercised on both supported bgms versions: an S7 object on bgms >= 0.2.0.0
  # and a plain S3 list on 0.1.6.3. The prior constructors only exist on the
  # newer version, so the SBM prior is specified in whichever form applies.
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:4]

  fit <- suppressWarnings(
    bgms::bgm(dat, iter = 50, warmup = 300, chains = 2, display_progress = FALSE)
  )
  # these two list methods used to fail on the missing `save` fit argument
  expect_no_error(suppressWarnings(plot_centrality(list(fit, fit))))
  expect_no_error(suppressWarnings(plot_prior_sensitivity(list(fit, fit))))

  # clusterBayesfactor used to call names()<- on the fit, which an S7 object
  # does not allow, and read $sbm, which a raw fit does not carry
  sbm_arg <- if (packageVersion("bgms") > "0.1.6.3") {
    bgms::sbm_prior()
  } else {
    "Stochastic-Block"
  }
  fit_sbm <- suppressWarnings(
    bgms::bgm(dat, edge_prior = sbm_arg, iter = 50, warmup = 300,
              chains = 2, display_progress = FALSE)
  )
  expect_no_error(suppressWarnings(clusterBayesfactor(fit_sbm)))
})

test_that("legacy interaction_scale does not raise a bgms deprecation warning", {
  skip_if(packageVersion("bgms") <= "0.1.6.3")
  set.seed(123)
  data("Wenchuan", package = "bgms")
  dat <- na.omit(Wenchuan)[1:40, 1:3]

  w <- character(0)
  withCallingHandlers(
    easybgm(dat, type = "ordinal", iter = 50, warmup = 300, progress = FALSE,
            interaction_scale = 2.5),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("deprecat", w, ignore.case = TRUE)))
})
