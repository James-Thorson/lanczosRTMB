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
obj2 = MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
opt2 = nlminb( obj2$par, obj2$fn, obj2$gr )
opt$par - opt2$par
#>            mu         logsd         logcv 
#> -2.380471e-05  4.309889e-05 -9.203910e-05 
```
