# Estimate variance using Lanczos method

Estimate squared standard error (variance) for a parameter or derived
quantity without forming the Hessian matrix directly

## Usage

``` r
lanczos_variance(
  Hv,
  v,
  k = c(25, 30),
  min_spectral_ratio = 1e-10,
  orthogonalize = FALSE
)
```

## Arguments

- Hv:

  function that calculates the product `Hv`

- v:

  vector to use when calculating variance, either an indicator for a
  single parameter, or a gradient evaluated at the MLE for a derived
  quantity

- k:

  can be a vector

- min_spectral_ratio:

  is the ratio of minimum to maximum ratio, where values \< minimum are
  truncated to minimum (min_spectral_ratio=0 disables truncation)

- orthogonalize:

  Whether to do two-pass Gram-Schmidt re-normalization (much slower)
