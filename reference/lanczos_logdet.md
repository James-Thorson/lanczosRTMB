# Estimate log-determinant likelihood using Lanczos method

Estimate log-determinant of inner Hessian matrix with a stochastic
approximation

## Usage

``` r
lanczos_logdet(Hv, k, m, n, seed = NULL, orthogonalize = TRUE)
```

## Arguments

- Hv:

  function that calculates the product `Hv`

- k:

  dimension for Kyrlov subspace

- m:

  number of probe-vectors to use for approximating average and standard
  deviation of log-determinant

- n:

  length of parameters (and necessary probe-vector)

- seed:

  if not NULL, then sets the seed. This is helfpul given that the
  Hutchinson probe vectors are randomly sampled, and comparisons have
  lower variance using a fixed seed.

- orthogonalize:

  Whether to do two-pass Gram-Schmidt re-normalization (much slower)

## Details

For a model with independent random effects, the variance of stochastic
trace estimation should be zero
