# Estimate log-determinant likelihood using Lanczos method

Estimate log-determinant of inner Hessian matrix with a stochastic
approximation

## Usage

``` r
lanczos_logdet(
  Hq,
  k,
  m,
  x = attr(Hq, "env")$x0,
  which_random = attr(Hq, "env")$which_random,
  seed = NULL,
  orthogonalize = TRUE,
  return_extra = FALSE
)
```

## Arguments

- Hq:

  function that calculates the product `H %*% q` given probe `q` and
  parameters `x`

- k:

  dimension for Kyrlov subspace

- m:

  number of probe-vectors to use for approximating average and standard
  deviation of log-determinant

- x:

  parameter vector used when calculating the Hessian matrix

- which_random:

  integer-vector indicating which elements of `x` correspond to random
  effects, where the probe `q` then has length `length(which_random)`

- seed:

  if not NULL, then sets the seed. This is helfpul given that the
  Hutchinson probe vectors are randomly sampled, and comparisons have
  lower variance using a fixed seed.

- orthogonalize:

  Whether to do two-pass Gram-Schmidt re-normalization (much slower)

- return_extra:

  whether to return probes and other internal constructions.

## Details

For a model with independent random effects, the variance of stochastic
trace estimation should be zero

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
  nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
  nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
  jnll = -1 * ( sum(nll1) + sum(nll2) )
  return(jnll)
}
params = list(u=u, mu = 0, logsd = 0, logcv = 0)

# Make RTMB object
obj = RTMB::MakeADFun( nll, params, random = "u", silent = TRUE )

# Make Lanczos object
tape = MakeTape( nll, params )
Hq = make_Hq( tape, unlist(params), which_random = 1:30 )

# Compare determinant at start values
lanczos_logdet( Hq, k = 10, m = 3 )
#> [1] 18.43706 18.43706 18.43706
H = obj$env$spHess(par = obj$env$par, random = TRUE)
sum(log(eigen(H)$values))
#> [1] 18.43706

# Compare determinant at new values
x_new = unlist(params)
  x_new['logsd'] = 1
lanczos_logdet( Hq, x = x_new, k = 10, m = 3 )
#> [1] -1.454869 -1.454869 -1.454869
H = obj$env$spHess(par = x_new, random = TRUE)
sum(log(eigen(H)$values))
#> [1] -1.454869
```
