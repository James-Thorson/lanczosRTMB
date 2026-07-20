# Package index

## Approximate properties

Core tools, using Lanczos to approximate model properties

- [`lanczos_logdet()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_logdet.md)
  : Estimate log-determinant likelihood using Lanczos method
- [`lanczos_nll()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_nll.md)
  : Estimate log-marginal likelihood using Lanczos method
- [`lanczos_sample()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_sample.md)
  : Sample from Lanczos method
- [`lanczos_variance()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_variance.md)
  : Estimate variance using Lanczos method
- [`lanczos_dgmrf()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_dgmrf.md)
  : Approximate the log-density of a GMRF using Lanczos

## Helper functions

Helper tools, used for components of an analysis

- [`lanczos()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos.md)
  : Calculate Lanczos approximation
- [`make_Hq()`](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md)
  : Make function to calculate H %\*% q

## Experimental and development

Experimental tools, under development.

- [`lanczos_MakeADFun()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_MakeADFun.md)
  : Approximate log-marginal likelihood using Lanczos method
- [`lanczos_fixedQ()`](https://james-thorson.github.io/lanczosRTMB/reference/lanczos_fixedQ.md)
  : Calculate Lanczos approximation using fixed Q
- [`CG()`](https://james-thorson.github.io/lanczosRTMB/reference/CG.md)
  : Truncated conjugate gradient using Hessian-vector products
- [`newton_CG()`](https://james-thorson.github.io/lanczosRTMB/reference/newton_CG.md)
  : Nonlinear minimizer using line search with approximated Newton
  solution
