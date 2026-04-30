# Sample from Lanczos method

Sample from a penalized likelihood model without directly forming the
Hessian matrix

## Usage

``` r
lanczos_sample(Hq, q, k, nsamp = 1, orthogonalize = FALSE)
```

## Arguments

- Hq:

  function that calculates the product `H %*% q`

- k:

  dimension for Kyrlov subspace

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
Hq = make_Hq( pen )

# Gradient-based Lanczos sampling
what = "sumexpu"
grad = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap, silent = TRUE )$gr( opt_pen$par )[1,]
what = "jnll"
sample_x = function(n){ x = rnorm(n); return( x / sqrt(sum(x^2))) }
samples = lanczos_sample(
  Hq = Hq,
  v = grad,
  k = 30,
  n = 1000,
  orthogonalize = TRUE
)
#> Error in lanczos_sample(Hq = Hq, v = grad, k = 30, n = 1000, orthogonalize = TRUE): unused argument (v = grad)

# Samples from Lanczos for parameters that contribute to derived quantity
samples = sweep( samples, MARGIN = 1, FUN = "+", STAT = opt$par )
#> Error: object 'samples' not found
apply( samples, MARGIN = 1, FUN = sd )
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'X' in selecting a method for function 'apply': object 'samples' not found

# Samples from full Hessian should approximately match
samp2 = t(mvtnorm::rmvnorm( n = 1000, sigma = as.matrix(solve(H)) ))
#> Error in loadNamespace(x): there is no package called ‘mvtnorm’
apply( samp2, MARGIN = 1, FUN = sd )
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'X' in selecting a method for function 'apply': object 'samp2' not found

# Compare bias-correction
sumexpu_z = apply( samples, MARGIN = 2, FUN = \(x) pen$report(x)$sumexpu )
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'X' in selecting a method for function 'apply': object 'samples' not found
mean(sumexpu_z)
#> Error: object 'sumexpu_z' not found
summary(sdrep)['sumexpu',]
#>            Estimate          Std. Error Est. (bias.correct) Std. (bias.correct) 
#>            4.064121            1.625808            4.734025                  NA 
```
