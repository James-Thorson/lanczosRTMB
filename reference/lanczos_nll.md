# Estimate log-marginal likelihood using Lanczos method

Estimate log-marginal likelihood using Laplace approximation, but
replacing exact calculation of log-determinant with a stochastic
approximation

## Usage

``` r
lanczos_nll(obj, k, m, Hq, seed = NULL)
```

## Arguments

- obj:

  TMB object (output from
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html))

- k:

  dimension for Kyrlov subspace

- m:

  number of probe-vectors to use for approximating average and standard
  deviation of log-determinant

- Hq:

  function that calculates the product `H %*% q`

- seed:

  if not NULL, then sets the seed. This is helfpul given that the
  Hutchinson probe vectors are randomly sampled, and comparisons have
  lower variance using a fixed seed.

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
Hq = make_Hq(pen)
lanczos_logdet( Hq, k = 10, m = 3, n = length(pen$par) )
#> [1] 44.53951 44.53951 44.53951
Matrix::determinant( H )
#> $modulus
#> [1] 44.4956
#> attr(,"logarithm")
#> [1] TRUE
#> 
#> $sign
#> [1] 1
#> 
#> attr(,"class")
#> [1] "det"

# Compare marginal likelihoods
lanczos_nll( pen, k = 10, m = 10 )
#>      nll   sd_nll 
#> 36.46908  0.00000 
opt$obj
#> [1] 36.46907
```
