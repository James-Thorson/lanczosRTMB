# Sample from Lanczos method

Sample from a penalized likelihood model without directly forming the
Hessian matrix

## Usage

``` r
lanczos_sample(Hv, v, k, nsamp = 1, orthogonalize = FALSE)
```

## Arguments

- Hv:

  function that calculates the product `Hv`

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
obj = RTMB::MakeADFun( nll, list(u=u, mu = 0, logsd = 0, logcv = 0), random = "u" )
opt = nlminb( obj$par, obj$fn, obj$gr )
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
#> iter: 4  mgc: 1.407763e-13 
#> iter: 1  mgc: 1.407763e-13 
#> outer mgc:  0.1645561 
#> iter: 1  value: 42.03535 mgc: 0.2020453 ustep: 1 
#> iter: 2  value: 42.03534 mgc: 0.004116722 ustep: 1 
#> iter: 3  value: 42.03534 mgc: 1.65679e-06 ustep: 1 
#> iter: 4  mgc: 2.68674e-13 
#> iter: 1  value: 42.37252 mgc: 0.0298795 ustep: 1 
#> iter: 2  value: 42.37252 mgc: 6.610381e-05 ustep: 1 
#> iter: 3  mgc: 4.124088e-10 
#> iter: 1  mgc: 4.124088e-10 
#> outer mgc:  0.1417209 
#> iter: 1  value: 42.27696 mgc: 0.003953647 ustep: 1 
#> iter: 2  value: 42.27696 mgc: 9.898509e-07 ustep: 1 
#> iter: 3  mgc: 9.148238e-14 
#> iter: 1  mgc: 9.148238e-14 
#> outer mgc:  0.02801003 
#> iter: 1  value: 42.11956 mgc: 0.0445184 ustep: 1 
#> iter: 2  value: 42.11956 mgc: 0.0001876571 ustep: 1 
#> iter: 3  mgc: 3.31473e-09 
#> iter: 1  mgc: 3.31473e-09 
#> outer mgc:  0.04623801 
#> iter: 1  value: 41.96495 mgc: 0.04784547 ustep: 1 
#> iter: 2  value: 41.96495 mgc: 0.0002183498 ustep: 1 
#> iter: 3  mgc: 4.517448e-09 
#> iter: 1  mgc: 4.517448e-09 
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
#> iter: 3  mgc: 9.992007e-16 
#> iter: 1  mgc: 9.992007e-16 
#> outer mgc:  3.242996e-05 
#> iter: 1  value: 41.78943 mgc: 3.211221e-06 ustep: 1 
#> iter: 2  mgc: 9.825474e-13 
#> iter: 1  mgc: 9.825474e-13 
#> outer mgc:  1.636354e-06 
#> iter: 1  mgc: 9.825474e-13 
sdrep = sdreport(obj, bias.correct = TRUE )
#> iter: 1  mgc: 9.825474e-13 
#> outer mgc:  1.636354e-06 
#> iter: 1  value: 41.79255 mgc: 0.001303796 ustep: 1 
#> iter: 2  value: 41.79255 mgc: 1.620694e-07 ustep: 1 
#> iter: 3  mgc: 2.88658e-15 
#> outer mgc:  0.0269305 
#> iter: 1  value: 41.78634 mgc: 0.001303796 ustep: 1 
#> iter: 2  value: 41.78634 mgc: 1.621195e-07 ustep: 1 
#> iter: 3  mgc: 3.552714e-15 
#> outer mgc:  0.0269297 
#> iter: 1  value: 41.79749 mgc: 0.003904606 ustep: 1 
#> iter: 2  value: 41.79749 mgc: 1.23556e-06 ustep: 1 
#> iter: 3  mgc: 1.452172e-13 
#> outer mgc:  0.0311446 
#> iter: 1  value: 41.78139 mgc: 0.003912423 ustep: 1 
#> iter: 2  value: 41.78139 mgc: 1.23487e-06 ustep: 1 
#> iter: 3  mgc: 1.449951e-13 
#> outer mgc:  0.03114147 
#> iter: 1  value: 41.81138 mgc: 0.003904606 ustep: 1 
#> iter: 2  value: 41.81138 mgc: 1.232402e-06 ustep: 1 
#> iter: 3  mgc: 1.44329e-13 
#> outer mgc:  0.01275106 
#> iter: 1  value: 41.76748 mgc: 0.003912423 ustep: 1 
#> iter: 2  value: 41.76748 mgc: 1.238033e-06 ustep: 1 
#> iter: 3  mgc: 1.458833e-13 
#> outer mgc:  0.01276952 
#> outer mgc:  2.704224 
#> Re-using symbolic Cholesky
#> iter: 1  mgc: 9.825474e-13 
#> Matching hessian patterns... Done
#> outer mgc:  4.734025 
H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)

# Re-do as penalized likelihood
newmap = list(mu = factor(NA), logsd = factor(NA), logcv = factor(NA))
pen = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap )
opt_pen = nlminb( pen$par, pen$fn, pen$gr )
#> outer mgc:  1.458833e-13 
Hv = make_Hv( pen )

# Gradient-based Lanczos sampling
what = "sumexpu"
grad = RTMB::MakeADFun( nll, obj$env$parList(), map = newmap )$gr( opt$par )[1,]
#> Error in EvalADFunObject(ADFun, theta, order = order, hessiancols = cols,     hessianrows = rows, sparsitypattern = sparsitypattern, rangecomponent = rangecomponent,     rangeweight = rangeweight, dumpstack = dumpstack, doforward = doforward,     set_tail = set_tail, data_changed = data_changed): Wrong parameter length.
what = "jnll"
sample_x = function(n){ x = rnorm(n); return( x / sqrt(sum(x^2))) }
samples = lanczos_sample(
  Hv = Hv,
  v = grad,
  k = 30,
  n = 1000,
  orthogonalize = TRUE
)
#> Error: object 'grad' not found

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
