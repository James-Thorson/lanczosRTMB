# Estimate log-determinant likelihood using Lanczos method

Estimate log-determinant of inner Hessian matrix with a stochastic
approximation

## Usage

``` r
lanczos_logdet(
  Hq,
  k,
  m,
  n,
  seed = NULL,
  orthogonalize = TRUE,
  return_extra = FALSE
)
```

## Arguments

- Hq:

  function that calculates the product `H %*% q`

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
  sumexpu = sum(exp(p$u[seq_len(n_sum)]))
  ADREPORT( sumexpu )
  REPORT( sumexpu )
  nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
  nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
  jnll = -1 * ( sum(nll1) + sum(nll2) )
  if(what == "jnll") return(jnll)
  if(what == "sumexpu") return(sumexpu)
}
obj = RTMB::MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u", silent = TRUE )
opt = nlminb( obj$par, obj$fn, obj$gr )
sdrep = sdreport(obj, bias.correct = TRUE )
H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)

# Re-do as penalized likelihood
newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
pen = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap, silent = TRUE )
opt_pen = nlminb( pen$par, pen$fn, pen$gr )

# Compare determinant
Hq = make_Hq( GetTape(pen), opt_pen$par )
lanczos_logdet( Hq, k = 10, m = 3, n = length(pen$par) )
#> [1] 44.53951 44.53951 44.53951
Matrix::determinant( H )$modulus
#> [1] 44.4956
#> attr(,"logarithm")
#> [1] TRUE
```
