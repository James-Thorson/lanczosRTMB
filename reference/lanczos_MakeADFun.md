# Approximate log-marginal likelihood using Lanczos method (EXPERIMENTAL)

Make a function that returns the Lanczos-Laplace approximation for a
user-supplied joint likelihood and designated random effects

## Usage

``` r
lanczos_MakeADFun(
  func,
  parameters,
  random,
  k,
  profile = NULL,
  m = 3,
  seed = 123
)
```

## Arguments

- func:

  Function taking a parameter list (or parameter vector) as input.

- parameters:

  Parameter list (or parameter vector) used by `func`.

- random:

  As [MakeADFun](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

- k:

  dimension for Kyrlov subspace

- profile:

  As [MakeADFun](https://rdrr.io/pkg/TMB/man/MakeADFun.html).

- m:

  number of probe-vectors to use for approximating average and standard
  deviation of log-determinant

- seed:

  if not NULL, then sets the seed. This is helfpul given that the
  Hutchinson probe vectors are randomly sampled, and comparisons have
  lower variance using a fixed seed.

## Value

An object (list) of class `tinyVAST`. Elements include:

- par:

  parameter-vector of fixed effects

- fn:

  function that returns the Lanczos-Laplace approximation given fixed
  effects

- env:

  environment of local variables

- Hv:

  function that returns the Hessian-vector product (for use in
  debugging)
