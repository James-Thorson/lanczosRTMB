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
  seed = 123,
  do_grad = FALSE
)
```

## Arguments

- func:

  Function taking a parameter list (or parameter vector) as input.

- parameters:

  Parameter list (or parameter vector) used by `func`.

- k:

  dimension for Kyrlov subspace

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

- Hq:

  function that returns the Hessian-vector product (for use in
  debugging)

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
obj2 = MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u" )
opt2 = nlminb( obj2$par, obj2$fn, obj2$gr )
#> iter: 1  value: 58.37851 mgc: 2.856111 ustep: 1 
#> iter: 2  value: 57.98311 mgc: 0.6945181 ustep: 1 
#> iter: 3  value: 57.98169 mgc: 0.04932098 ustep: 1 
#> iter: 4  value: 57.98169 mgc: 0.0002842405 ustep: 1 
#> iter: 5  mgc: 9.465661e-09 
#> iter: 1  mgc: 9.465661e-09 
#> Matching hessian patterns... Done
#> outer mgc:  11.72122 
#> iter: 1  value: 39.77521 mgc: 3.895884 ustep: 1 
#> iter: 2  value: 39.74061 mgc: 0.3920337 ustep: 1 
#> iter: 3  value: 39.74061 mgc: 0.005037983 ustep: 1 
#> iter: 4  value: 39.74061 mgc: 8.569536e-07 ustep: 1 
#> iter: 5  mgc: 2.620126e-14 
#> iter: 1  value: 46.39835 mgc: 0.7570035 ustep: 1 
#> iter: 2  value: 46.39798 mgc: 0.03489372 ustep: 1 
#> iter: 3  value: 46.39798 mgc: 8.200935e-05 ustep: 1 
#> iter: 4  mgc: 4.552982e-10 
#> iter: 1  mgc: 4.552982e-10 
#> outer mgc:  1.84157 
#> iter: 1  value: 40.59271 mgc: 0.7974252 ustep: 1 
#> iter: 2  value: 40.5909 mgc: 0.0294835 ustep: 1 
#> iter: 3  value: 40.5909 mgc: 3.901207e-05 ustep: 1 
#> iter: 4  mgc: 7.772716e-11 
#> iter: 1  value: 44.68845 mgc: 0.1246689 ustep: 1 
#> iter: 2  value: 44.68844 mgc: 0.0009686411 ustep: 1 
#> iter: 3  value: 44.68844 mgc: 5.984028e-08 ustep: 1 
#> iter: 4  mgc: 1.44329e-15 
#> iter: 1  mgc: 1.44329e-15 
#> outer mgc:  0.5755699 
#> iter: 1  value: 43.92098 mgc: 0.380674 ustep: 1 
#> iter: 2  value: 43.92087 mgc: 0.01093109 ustep: 1 
#> iter: 3  value: 43.92087 mgc: 8.978169e-06 ustep: 1 
#> iter: 4  mgc: 6.059375e-12 
#> iter: 1  mgc: 6.059375e-12 
#> outer mgc:  0.71546 
#> iter: 1  value: 43.01585 mgc: 0.5486935 ustep: 1 
#> iter: 2  value: 43.01565 mgc: 0.01743663 ustep: 1 
#> iter: 3  value: 43.01565 mgc: 2.527355e-05 ustep: 1 
#> iter: 4  mgc: 5.775047e-11 
#> iter: 1  mgc: 5.775047e-11 
#> outer mgc:  0.359416 
#> iter: 1  value: 42.25048 mgc: 0.1961147 ustep: 1 
#> iter: 2  value: 42.25047 mgc: 0.003648652 ustep: 1 
#> iter: 3  value: 42.25047 mgc: 1.232595e-06 ustep: 1 
#> iter: 4  mgc: 1.403322e-13 
#> iter: 1  mgc: 1.403322e-13 
#> outer mgc:  0.1645561 
#> iter: 1  value: 42.03535 mgc: 0.2020453 ustep: 1 
#> iter: 2  value: 42.03534 mgc: 0.004116722 ustep: 1 
#> iter: 3  value: 42.03534 mgc: 1.65679e-06 ustep: 1 
#> iter: 4  mgc: 2.68896e-13 
#> iter: 1  value: 42.37252 mgc: 0.0298795 ustep: 1 
#> iter: 2  value: 42.37252 mgc: 6.610381e-05 ustep: 1 
#> iter: 3  mgc: 4.124083e-10 
#> iter: 1  mgc: 4.124083e-10 
#> outer mgc:  0.1417209 
#> iter: 1  value: 42.27696 mgc: 0.003953647 ustep: 1 
#> iter: 2  value: 42.27696 mgc: 9.898509e-07 ustep: 1 
#> iter: 3  mgc: 9.192647e-14 
#> iter: 1  mgc: 9.192647e-14 
#> outer mgc:  0.02801003 
#> iter: 1  value: 42.11956 mgc: 0.0445184 ustep: 1 
#> iter: 2  value: 42.11956 mgc: 0.0001876571 ustep: 1 
#> iter: 3  mgc: 3.314731e-09 
#> iter: 1  mgc: 3.314731e-09 
#> outer mgc:  0.04623801 
#> iter: 1  value: 41.96495 mgc: 0.04784547 ustep: 1 
#> iter: 2  value: 41.96495 mgc: 0.0002183498 ustep: 1 
#> iter: 3  mgc: 4.517449e-09 
#> iter: 1  mgc: 4.517449e-09 
#> outer mgc:  0.02447321 
#> iter: 1  value: 41.81629 mgc: 0.04627014 ustep: 1 
#> iter: 2  value: 41.81629 mgc: 0.0002053923 ustep: 1 
#> iter: 3  mgc: 4.020242e-09 
#> iter: 1  mgc: 4.020242e-09 
#> outer mgc:  0.002854085 
#> iter: 1  value: 41.79253 mgc: 0.006424874 ustep: 1 
#> iter: 2  value: 41.79253 mgc: 3.939604e-06 ustep: 1 
#> iter: 3  mgc: 1.480149e-12 
#> iter: 1  mgc: 1.480149e-12 
#> outer mgc:  0.0002425434 
#> iter: 1  value: 41.78946 mgc: 0.0008453203 ustep: 1 
#> iter: 2  value: 41.78946 mgc: 6.814674e-08 ustep: 1 
#> iter: 3  mgc: 1.998401e-15 
#> iter: 1  mgc: 1.998401e-15 
#> outer mgc:  3.242996e-05 
#> iter: 1  value: 41.78943 mgc: 3.211222e-06 ustep: 1 
#> iter: 2  mgc: 9.829915e-13 
#> iter: 1  mgc: 9.829915e-13 
#> outer mgc:  1.636354e-06 
#> iter: 1  mgc: 9.829915e-13 
opt$par - opt2$par
#>            mu         logsd         logcv 
#>  3.931741e-05 -5.483657e-05  1.143543e-04 
```
