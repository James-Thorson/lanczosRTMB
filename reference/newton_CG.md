# Nonlinear minimizer using line search with approximated Newton solution

Nonlinear minimizer designed for cheap Hessian-vector products using
[make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md),
involving iterating a linear search along the truncated conjugate
gradient for a Newton.

## Usage

``` r
newton_CG(
  par,
  fn,
  gr,
  Hq,
  gr_tol = 1e-08,
  e_ratio = 0.01,
  maxit_newton = 100,
  maxit_CG = min(100, length(par)),
  c1 = 0.01,
  beta = 0.5,
  silent = FALSE
)
```

## Arguments

- par:

  initial parameter vector

- fn:

  function to evaluate negative log-likelihood

- gr:

  function to evaluate gradient of negative log-likelihood

- Hq:

  efficient Hessian-vector product function, e.g.,
  [make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md)

- gr_tol:

  early stopping condition for gradient of Newton solver

- e_ratio:

  early stopping condition for error of CG, relative to initial gradient

- maxit_newton:

  maximum iterations for Newton solver

- maxit_CG:

  maximum iterations for CG solution for each Newton iteration

- c1:

  stopping condition for line search given CG solution in each Newton
  iteration

- beta:

  updates in line search stepsize alpha when Armijo sufficient decrease
  condition fails

- silent:

  Be silent or print progress?
