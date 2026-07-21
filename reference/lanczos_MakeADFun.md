# Approximate log-marginal likelihood using Lanczos method

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
  method = "newton_CG",
  seed = 123,
  make_gr = TRUE,
  pu_update = c("FD", "exact"),
  silent = TRUE
)
```

## Arguments

- func:

  Function taking a parameter list (or parameter vector) as input.

- parameters:

  Parameter list (or parameter vector) used by `func`.

- random:

  Character vector defining the random effect parameters. See also
  `regexp`.

- k:

  dimension for Kyrlov subspace

- profile:

  Parameters to profile out of the likelihood (this subset will be
  appended to `random` with Laplace approximation disabled).

- m:

  number of probe-vectors to use for approximating average and standard
  deviation of log-determinant

- method:

  whether to use
  [newton_CG](https://james-thorson.github.io/lanczosRTMB/reference/newton_CG.md)
  or a gradient-based low-memory option specifically "L-BFGS-B" in
  [optim](https://rdrr.io/r/stats/optim.html) to optimize the inner
  problem

- seed:

  if not NULL, then sets the seed. This is helfpul given that the
  Hutchinson probe vectors are randomly sampled, and comparisons have
  lower variance using a fixed seed.

- make_gr:

  whether to make approximated gradient using fixed probes (slow for
  large models)

- pu_update:

  when make_gr=TRUE, whether to use a finite-difference based on an AD
  tape or re-optimize the joint likelihood to get the update on random
  effects when calculating the FD for fixed effects in the
  log-determinant calculation. pu_update="FD" can be very slow for dense
  inner-Hessians.

- silent:

  Disable all tracing information?

## Value

An object (list), where elements include:

- par:

  parameter-vector of fixed effects

- fn:

  function that returns the Lanczos-Laplace approximation given fixed
  effects.

- gr:

  a function that returns the approximated gradient of `fn` with respect
  to fixed effects.

- env:

  environment of local variables

- Hq:

  function that returns the Hessian-vector product (for use in
  debugging)

## Details

The gradient uses a finite-difference applied to a fixed set of probes,
inspired by Dong et al. (2017) and the `stochasticLQ` option in
GPyTorch. Exploration suggests that a useful approximation to the
gradient of the log-marginal likelihood with respect to fixed effects
can be calculated using the gradient of the joint likelihood with
respect to the fixed effects, and a finite-difference approximation to
the log-determinant. The latter approximation is only performs well when
recalculating the Lanczos matrix \\Q\\.

## References

Dong, K., Eriksson, D., Nickisch, H., Bindel, D., & Wilson, A. G.
(2017). Scalable log determinants for Gaussian process kernel learning.
Advances in Neural Information Processing Systems, 30.
<https://proceedings.neurips.cc/paper/2017/hash/976abf49974d4686f87192efa0513ae0-Abstract.html>

## Examples

``` r
# Simulate lognormal-gamma process
set.seed(123)
library(RTMB)
n = 30
n_sum = 3
u = 0 + rnorm(n)
y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )

# Fit as GLMM
what = "jnll"
nll = function(p){
  sumexpu = sum(exp(p$u[seq_len(n_sum)]))
  nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
  nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
  jnll = -1 * ( sum(nll1) + sum(nll2) )
  return(jnll)
}

# Build
obj = lanczos_MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", k = 10 )
opt = nlminb( obj$par, obj$fn )

# Compare with RTMB
obj2 = MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
opt2 = nlminb( obj2$par, obj2$fn, obj2$gr )
opt$par - opt2$par
#>         mu      logsd      logcv 
#> -0.1585838  0.1860306 -6.1881942 

# Fit again using FD gradient for Lanczos method using fixed probe-recursion
opt3 = optim( obj$par, obj$fn, obj$gr, method = "BFGS" )
opt3$par - opt2$par
#>            mu         logsd         logcv 
#>  1.673518e-05 -2.399284e-05  5.883783e-05 
```
