# Sample from Lanczos method

Sample from a penalized likelihood model without directly forming the
Hessian matrix

## Usage

``` r
lanczos_sample(
  Hq,
  q,
  k,
  x = attr(Hq, "env")$x0,
  nsamp = 1,
  orthogonalize = FALSE
)
```

## Arguments

- Hq:

  function that calculates the product `H %*% q` given probe `q` and
  parameters `x`

- q:

  initial vector used for defining a Kyrlov subspace

- k:

  dimension for Kyrlov subspace

- x:

  parameter vector used when calculating the Hessian matrix

- nsamp:

  number of samples

- orthogonalize:

  Whether to do two-pass Gram-Schmidt re-normalization (much slower)

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
Hq = make_Hq( GetTape(pen), opt_pen$par )

# Gradient-based Lanczos sampling
what = "sumexpu"
grad = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap, silent = TRUE )$gr( opt_pen$par )[1,]
what = "jnll"
sample_x = function(n){ x = rnorm(n); return( x / sqrt(sum(x^2))) }
samples = lanczos_sample(
  Hq = Hq,
  q = grad,
  k = 30,
  n = 1000,
  orthogonalize = TRUE
)

# Samples from Lanczos for parameters that contribute to derived quantity
samples = sweep( samples, MARGIN = 1, FUN = "+", STATS = opt_pen$par )
apply( samples, MARGIN = 1, FUN = sd )
#>  [1] 0.4785608 0.4843596 0.4114760 0.0000000 0.0000000 0.0000000 0.0000000
#>  [8] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [15] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [22] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [29] 0.0000000 0.0000000

# Samples from full Hessian should approximately match
samp2 = t(mvtnorm::rmvnorm( n = 1000, sigma = as.matrix(solve(H)) ))
apply( samp2, MARGIN = 1, FUN = sd )
#>  [1] 0.4909493 0.4696421 0.4120659 0.4423042 0.4643884 0.3977407 0.4334791
#>  [8] 0.5479412 0.5325241 0.4876392 0.4092560 0.4274928 0.4818361 0.4363509
#> [15] 0.5039307 0.3924580 0.4678129 0.6049901 0.4369197 0.5336944 0.5921870
#> [22] 0.4993766 0.5171501 0.4881900 0.5038213 0.5131423 0.4787957 0.4439222
#> [29] 0.5402735 0.4108968

# Compare bias-correction
sumexpu_z = apply( samples, MARGIN = 2, FUN = \(x) pen$report(x)$sumexpu )
mean(sumexpu_z)
#> [1] 4.473128
summary(sdrep)['sumexpu',]
#>            Estimate          Std. Error Est. (bias.correct) Std. (bias.correct) 
#>            4.064121            1.625808            4.734025                  NA 
```
