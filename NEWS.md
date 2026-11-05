# csmGmm 0.6.0

## Breaking changes

* `sameDirAlt` no longer has a default value in `symm_fit_ind_EM()` and
  `symm_fit_ind_EM_noAssumption()` (previously `FALSE`). Users must now set
  `sameDirAlt = TRUE` for replication analysis or `sameDirAlt = FALSE` for
  mediation and pleiotropy analysis. This prevents replication analyses from
  accidentally using the larger (non-replication) alternative space.
* `symm_fit_cor_EM()`, `symm_fit_cor_EM_noAssumption()`,
  `symm_fit_cor_EM_rho()`, and `symm_fit_cor_EM_fulllik()` gain a required
  `sameDirAlt` argument, placed before `eps`. Code that passed `eps` or
  `checkpoint` by position should now pass them by name.
* The correlated-model functions (`symm_fit_cor_EM()`,
  `symm_fit_cor_EM_noAssumption()`, `symm_fit_cor_EM_rho()`, and
  `symm_fit_cor_EM_fulllik()`) now check that `testStats` is a matrix with
  exactly two columns and stop with an informative error otherwise.

## New features

* New `plot_fit_hist()` for checking model fit. For each dimension, it plots a
  histogram of the observed test statistics with the fitted marginal mixture
  density, the null N(0, 1) component, and the combined alternative component
  overlaid, plus the observed counts minus the counts expected under the fitted
  null. Optional `nBins`, `zLim`, and `yLim` arguments control the display.
* New `get_marginal_mixture()` collapses a fitted csmGmm into the
  one-dimensional marginal Gaussian mixture (means and weights) for each
  dimension. It is used internally by `plot_fit_hist()`.

## Documentation and minor changes

* Examples for all model-fitting functions now simulate signals in both
  directions and pass `sameDirAlt` explicitly.
* The vignette now sets `sameDirAlt` explicitly in every model-fitting call.

