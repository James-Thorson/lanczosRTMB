# Package index

## Core function

User-level function that combines computational and optimization tools

- [`lanczos_MakeADFun()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_MakeADFun.md)
  : Approximate log-marginal likelihood using Lanczos method

## Approximate properties

Copmutational tools, using Lanczos to approximate model properties

- [`lanczos_logdet()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_logdet.md)
  : Estimate log-determinant likelihood using Lanczos method
- [`lanczos_nll()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_nll.md)
  : Estimate log-marginal likelihood using Lanczos method
- [`lanczos_sample()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_sample.md)
  : Sample from Lanczos method
- [`lanczos_variance()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_variance.md)
  : Estimate variance using Lanczos method
- [`dgmrf_lanczos()`](https://james-thorson.github.io/lanczosRTMB/reference/dgmrf_lanczos.md)
  : Approximate the log-density of a GMRF using Lanczos

## Optimization tools

Optimization tools using Krylov subspaces

- [`CG()`](https://james-thorson.github.io/lanczosRTMB/reference/CG.md)
  : Truncated conjugate gradient using Hessian-vector products
- [`newton_CG()`](https://james-thorson.github.io/lanczosRTMB/reference/newton_CG.md)
  : Nonlinear minimizer using line search with approximated Newton
  solution
- [`dgmrf_icpl()`](https://james-thorson.github.io/lanczosRTMB/reference/dgmrf_icpl.md)
  : Approximate the log-density of a GMRF using incomplete Cholesky

## Helper functions

Helper tools, used for components of an analysis

- [`lanczos()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos.md)
  : Calculate Lanczos approximation
- [`make_Hq()`](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md)
  : Make function to calculate H %\*% q
- [`estimate_diag_hutchinson()`](https://james-thorson.github.io/lanczosRTMB/reference/estimate_diag_hutchinson.md)
  : Estimate diag(H) using Hutchinson probes

## Experimental and development

Experimental tools, under development.

- [`lanczos_fixedQ()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_fixedQ.md)
  : Calculate Lanczos approximation using fixed Q
