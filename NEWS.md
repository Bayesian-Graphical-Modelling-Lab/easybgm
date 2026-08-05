# easybgm 0.5.0

## Support for bgms >= 0.2.0.0

* easybgm now works with the S7 fit objects returned by bgms >= 0.2.0.0. This
  switches off the temporary backward-compatibility shim that bgms shipped for
  easybgm 0.4.0, so fits are no longer converted to S3 lists and the per-fit
  compatibility warning is gone.
* `bgms` is now the default fitting package for every data type, including
  `type = "continuous"` and `type = "mixed"`, which previously defaulted to
  BGGM. Pass `package = "BGGM"` or `package = "BDgraph"` to keep the old
  behaviour.
* `type` may now be given as a per-variable character vector, for example
  `type = c("ordinal", "ordinal", "continuous")`.
* Priors may now be given as bgms prior objects (`cauchy_prior()`,
  `normal_prior()`, `beta_prime_prior()`, `bernoulli_prior()`,
  `beta_bernoulli_prior()`, `sbm_prior()`). The older flat arguments still work
  and are translated internally.

## Blume-Capel main effects

* Fits with at least one Blume-Capel variable now return
  `blume_capel_parameters`, a data frame holding the posterior mean, posterior
  standard deviation, 95% credible interval and R-hat of the linear and
  quadratic effect of each Blume-Capel variable, together with the baseline
  category it was fitted with. Unlike the category thresholds of an ordinal
  variable, these two parameters are usually of substantive interest, so they
  are also printed by `summary()` rather than left in the fit object.
* With `save = TRUE`, the posterior draws of those effects are returned in
  `samples_blume_capel`.
* Baseline categories are reported on the scale of the input data. bgms recodes
  discrete scores to start at 0 and shifts the baseline category with them, so
  the value it stores internally can be lower than the one the user supplied.

## Bug fixes

* For Blume-Capel variables, the two columns of `thresholds` were labelled
  `cat (1)` and `cat (2)`, the same headers bgms uses for genuine category
  thresholds. They are in fact the linear and quadratic effect, and are now
  named accordingly. Where Blume-Capel and ordinal variables share one matrix
  the headers cannot describe both, so the per-row meaning is recorded in the
  matrix's `"variable_type"` attribute.
* `print()` on an unsummarised `easybgm` object printed the closing notes twice,
  once from the summary it prints internally and once from its own tail.

* `centrality` for bgms fits was computed from a mis-permuted edge matrix: the
  posterior samples were read back in BGGM's upper-triangle order rather than
  the lower-triangle order bgms uses. Per-node strengths were therefore
  permuted, and `plot_centrality()` reported them under the wrong node labels.
  The ordering is now stated explicitly at each call site.
* `structure` was returned as a complete graph (a matrix of ones) whenever
  `save = FALSE`, which is the default. `plot_structure()` consequently drew a
  fully connected network. It is now the median probability model in both
  branches, matching the documentation.
* The Monte Carlo interval in `MCSE_BF` mixed two estimators: it took the
  binomial variance of the raw indicator average but divided it by the
  effective sample size of the Rao-Blackwellized chain, and attached the result
  to a Rao-Blackwellized Bayes factor. The interval was too wide by up to about
  40%. It is now computed from the Monte Carlo standard error that bgms reports
  for the Rao-Blackwellized inclusion probability.
* `plot_centrality()` and `plot_prior_sensitivity()` failed on lists of raw
  bgms fit objects, because bgms no longer reports `save` among the fit
  arguments. Both now work, and both record the model type correctly.
* `clusterBayesfactor()` failed on a raw bgms fit object with "invalid to use
  names()<- on an S4 object". It now reads the prior and the block posterior
  through the bgms extractor functions, and gives an informative error when the
  fit was not estimated with the Stochastic Block Model prior.
* The legacy `interaction_scale` argument no longer leaks a bgms deprecation
  warning; it is translated to `cauchy_prior(scale)` like `pairwise_scale`.
* Corrected the documented defaults for `interaction_prior` and
  `precision_scale_prior`, and documented `precision_graph_prior` and
  `difference_family`.

## Results that change

Fitting with bgms >= 0.2.0.0 changes several numbers relative to easybgm 0.4.0
with bgms 0.1.6.3. None of these is a bug in either package:

* **Edge weights are about half their former size.** bgms now reports pairwise
  parameters on the association scale, the coefficient entering each
  conditional as `2 * omega * x`, where 0.1.6.3 stored `2 * omega`. This
  affects `parameters`, `samples_posterior`, `centrality`, and every plot drawn
  from them.
* **Edge weights are not on the same scale across fitting packages.** BGGM and
  BDgraph report partial correlations; bgms reports the pairwise association
  parameter. For bgms fits of continuous and mixed data the partial
  correlations and the precision matrix are returned separately, in
  `partial_correlations` and `precision_matrix`, and `summary()` now states
  which scale the reported edge weights are on.
* **Inclusion probabilities and Bayes factors are Rao-Blackwellized.** They no
  longer saturate at 0 or 1 on short chains, so inclusion Bayes factors are
  finite where they used to be `0` or `Inf`, and an edge can cross the
  median-probability threshold differently than before.
* **The default interaction prior changed** from a Cauchy to
  `normal_prior(scale = 1)`, on the new coordinate. A Normal slab has much
  lighter tails than a Cauchy and constrains weakly identified edges more
  tightly.
* **`convergence_parameter` is the classic split-R-hat.** bgms 0.1.6.3 applied
  a degrees-of-freedom adjustment that reported about 1.29 on nearly saturated
  indicators, that is, on the most decisive edges. Those now report near 1. `NA`
  and `Inf` are possible when all chains are identical or stuck.
* **Group comparisons keep every category any group observes.** bgms 0.1.6.3
  merged categories that were not observed in every group, which biased the
  affected variable's pairwise parameters. Results move most where groups have
  unequal category support.

## Other changes

* easybgm continues to support bgms 0.1.6.3, as stated in `DESCRIPTION`. The
  fixes above apply to both bgms versions, and the test suite now exercises
  them on 0.1.6.3 as well as on 0.2.0.0 rather than skipping them.
* The examples and tests now pass `warmup` explicitly. bgms defaults to
  `warmup = 2000` regardless of `iter`, which dominated the runtime of the
  examples. `warmup` is a `bgm()` argument in both supported bgms versions.
* Removed the unused `vdiffr` dependency and the `LazyData` field, and dropped
  some dead version-gating code.
