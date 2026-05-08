# Lanczos for penalized likelihood models

This vignette demonstrates how Lanczos methods can be used to
approximate standard errors and epsilon bias-correction for penalized
likelihood models, without directly forming or inverting the Hessian
matrix representing associations among random effects.

``` r

library(lanczosRTMB)
#> Loading required package: RTMB
library(RTMB)
```

## Simulate and fit a generalized linear mixed model

We first simulate data from a compound lognormal-Gamma process:

``` r

set.seed(123)
n = 30
n_sum = 3

# log-linked normal deviates
u = 0 + rnorm(n)

# Conditional gamma process
y = rgamma( n, shape = 1/0.5^2, scale = exp(u) * 0.5^2 )
```

We then define a function to compute the joint likelihood, while also
allowing us to calculate other outputs for later use:

``` r

nll = function(p){
  sumexpu = sum(exp(p$u[seq_len(n_sum)]))
  ADREPORT( sumexpu )
  REPORT( sumexpu )
  nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
  nll2 = dgamma( 
    y, 
    shape = 1/exp(2*p$logcv), 
    scale = exp(p$u) * exp(2*p$logcv), 
    log=TRUE
  )
  jnll = -1 * ( sum(nll1) + sum(nll2) )
  if(what == "jnll") return(jnll)
  if(what == "sumexpu") return(sumexpu)
  if(what == "biascorr") return(jnll + p$eps * sumexpu)
}

# Starting list of parameters
parlist = list(u=u*0, mu = 0, logsd = 0, logcv = 0, eps = 0)
```

Finally, we fit this model as a log-linked generalized linear mixed
model (GLMM):

``` r

# Build RTMB object
what = "jnll"
obj = RTMB::MakeADFun( 
  nll, 
  parlist,
  map = list(eps = factor(NA)),
  random = "u",
  silent = TRUE
)

# optimize 
opt = nlminb( obj$par, obj$fn, obj$gr )
opt
#> $par
#>         mu      logsd      logcv 
#> -0.1052225 -0.1326400 -0.5984998 
#> 
#> $objective
#> [1] 36.46907
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 13
#> 
#> $evaluations
#> function gradient 
#>       17       14 
#> 
#> $message
#> [1] "relative convergence (4)"
```

## Fit as marginal likelihood without Cholesky decomposition

We first show that we can re-fit the model using the marginal
likelihood, evaluated using Hutchinson-Lanczos instead of the Cholesky
factorization that is default in TMB:

``` r

parlist$eps = numeric()
pen0 = lanczos_MakeADFun(
  nll, 
  parlist,
  random = "u",
  method = "optim",
  k = 10,
  silent = TRUE
)

# optimize 
opt_pen0 = nlminb( pen0$par, pen0$fn, control = list(trace = 1) )
#>   0:     39.910646:  0.00000  0.00000  0.00000
#>   1:     36.513207: -0.0610391 -0.219325 -0.433051
#>   2:     36.499397: -0.0358829 -0.212092 -0.474385
#>   3:     36.483784: -0.0669697 -0.176820 -0.487917
#>   4:     36.474964: -0.0886213 -0.172405 -0.531567
#>   5:     36.473660: -0.0802839 -0.138859 -0.566191
#>   6:     36.469332: -0.0998972 -0.139569 -0.579935
#>   7:     36.469146: -0.101806 -0.137720 -0.588132
#>   8:     36.469110: -0.102884 -0.135692 -0.590837
#>   9:     36.469108: -0.102718 -0.135189 -0.591376
#>  10:     36.469102: -0.103304 -0.135192 -0.591852
#>  11:     36.469098: -0.103080 -0.135718 -0.592347
#>  12:     36.469097: -0.103158 -0.135528 -0.592121
#>  13:     36.469096: -0.103277 -0.135326 -0.592318
#>  14:     36.469093: -0.103770 -0.134572 -0.593141
#>  15:     36.469085: -0.103722 -0.134483 -0.593928
#>  16:     36.469080: -0.103931 -0.134178 -0.594779
#>  17:     36.469080: -0.104000 -0.134230 -0.594880
#>  18:     36.469079: -0.104071 -0.134293 -0.595092
#>  19:     36.469079: -0.104086 -0.134299 -0.595157
#>  20:     36.469078: -0.104125 -0.134036 -0.595166
#>  21:     36.469078: -0.104272 -0.133969 -0.595378
#>  22:     36.469077: -0.104121 -0.133968 -0.595598
#>  23:     36.469077: -0.104259 -0.133848 -0.595621
#>  24:     36.469076: -0.104399 -0.133814 -0.595735
#>  25:     36.469076: -0.104436 -0.133763 -0.595907
#>  26:     36.469074: -0.104640 -0.133447 -0.596539
#>  27:     36.469072: -0.105440 -0.132216 -0.599087
#>  28:     36.469072: -0.105494 -0.132188 -0.599348
#>  29:     36.469072: -0.105489 -0.132320 -0.599373
#>  30:     36.469072: -0.105498 -0.132278 -0.599335
#>  31:     36.469072: -0.105500 -0.132274 -0.599278
#>  32:     36.469072: -0.105491 -0.132299 -0.599289
#>  33:     36.469072: -0.105467 -0.132301 -0.599273
#>  34:     36.469072: -0.105474 -0.132291 -0.599265
#>  35:     36.469072: -0.105470 -0.132297 -0.599253
#>  36:     36.469072: -0.105477 -0.132350 -0.599233
#>  37:     36.469072: -0.105464 -0.132339 -0.599178
#>  38:     36.469072: -0.105462 -0.132335 -0.599181
#>  39:     36.469072: -0.105457 -0.132335 -0.599178
#>  40:     36.469072: -0.105447 -0.132337 -0.599171
#>  41:     36.469072: -0.105447 -0.132339 -0.599171
#>  42:     36.469072: -0.105447 -0.132339 -0.599171
#>  43:     36.469072: -0.105447 -0.132339 -0.599171
opt_pen0
#> $par
#>         mu      logsd      logcv 
#> -0.1054474 -0.1323386 -0.5991714 
#> 
#> $objective
#> [1] 36.46907
#> 
#> $convergence
#> [1] 1
#> 
#> $iterations
#> [1] 43
#> 
#> $evaluations
#> function gradient 
#>       73      228 
#> 
#> $message
#> [1] "false convergence (8)"
```

## Fit as a penalized likelihood model

Next, we show that the model can be refitted using penalized likelihood,
conditional upon fixed values for variance parameters:

``` r

# Define RTMB object for penalized likelihood
newmap = list(
  mu = factor(NA), 
  logsd = factor(NA), 
  logcv = factor(NA),
  eps = factor(NA)
)
pen = RTMB::MakeADFun( 
  nll, 
  obj$env$parList(), 
  map = newmap,
  silent = TRUE
)

# Re-optimize
opt_pen = nlminb( pen$par, pen$fn, pen$gr )
```

Alternatively, we provide a Newton optimizer using a truncated conjugate
gradient, which operates using Hessian-vector products without ever
forming the full Hessian matrix:

``` r

# Extract tape and make Hessian-vector-product function
tape = GetTape(pen)
Hq = make_Hq( tape, pen$par )

# run Newton optimizer
opt_pen2 = newton_CG(
  par = pen$par,
  fn = tape,
  gr = tape$jacfun(),
  Hq = Hq
)
#> value: 41.78943 mgc: 9.836576e-13
```

We can then use stochastic trace estimation and the Lanczos method to
approximate the log-determinant of the Hessian for random effects:

``` r

# log-determinant using Lanczos
Hq = make_Hq( GetTape(pen), opt_pen$par )
lanczos_logdet( Hq, k = 10, m = 3 )
#> [1] 44.4956 44.4956 44.4956

# log-determinant for marginal likelihood
H = obj$env$spHess(par = obj$env$last.par.best, random = TRUE)
Matrix::determinant( H )$modulus
#> [1] 44.4956
#> attr(,"logarithm")
#> [1] TRUE
```

Using this log-determinant, we can also compare the log-marginal
likelihood using Lanczos with the full calculation:

``` r

# Marginal likelihoods using Lanczos
lanczos_nll( pen, k = 10, m = 10 )
#>      nll   sd_nll 
#> 36.46907  0.00000

# Marginal likelihood using Laplace
opt$obj
#> [1] 36.46907
```

## Delta methods

Alternatively, we can apply Lanczos methods to standard errors for
parameters. Here, we will sample the first three random effects

``` r

samples = lanczos_sample(
   Hq = Hq,
   q = c( rep(1,3), rep(0,n-3) ),
   k = 30,
   n = 1000,
   orthogonalize = TRUE
)

#
sdrep = sdreport( obj, ignore.parm.uncertainty = TRUE )
as.list(sdrep, what = "Std. Error")$u[1:3]
#> [1] 0.4839019 0.4864751 0.4066186
apply( samples, MARGIN = 1, FUN = sd)[1:3]
#> [1] 0.4782002 0.4887318 0.4051458
```

Alternatively, we can apply Lanczos methods to approximate the Hessian
with respect to a derived quantity. This involves calculating the
gradient for the derived quantity with respect to coefficients:

``` r

what = "sumexpu"
grad = RTMB::MakeADFun( 
  nll, 
  obj$env$parList(), 
  map = newmap,
  silent = TRUE
)$gr( pen$par )[1,]
```

We can then use this gradient to approximate samples from the Hessian
matrix and use that to calculate a Monte Carlo estimator for standard
errors:

``` r

samples = lanczos_sample(
   Hq = Hq,
   q = grad,
   k = 30,
   n = 1000,
   orthogonalize = TRUE
)

# Compare SD of samples with the delta method
sumexpu1 = apply( samples, MARGIN = 2, FUN = \(v)pen$report(v)$sumexpu )
c( pen$report()$sumexpu, sd(sumexpu1) )
#> [1] 4.0641214 0.9435151
summary(sdrep)['sumexpu',]
#>   Estimate Std. Error 
#>   4.064121   1.194484
```

Alternatively, we can calculate the standard error for a derived
quantity using Lanczos within the delta method:

``` r

# 
Var = lanczos_variance(
  Hq = Hq,
  q = grad,
  k = 3
)

# Compare with the delta method
c( pen$report()$sumexpu, sqrt(Var) )
#> [1] 4.064121 1.194484
summary(sdrep)['sumexpu',]
#>   Estimate Std. Error 
#>   4.064121   1.194484
```

## Epsilon methods

Finally, we can use the gradient of the log-likelihood with respect to a
dummy variable epsilon to correct for retransformation bias:

``` r

# One-sided finite difference
what = "biascorr"
phat = obj$env$parList()
phat$eps = 0.0001
pen_hi = RTMB::MakeADFun( 
  nll, 
  parameters = phat, 
  map = newmap,
  silent = TRUE
)

# Re-optimize
opt_hi = optim( 
  pen_hi$par, pen_hi$fn, pen_hi$gr, 
  method = "L-BFGS-B", control = list(factr = 1e-2) 
)

# Calculate marginal nll for finite-difference
nll_hi = lanczos_nll( pen_hi, k = 10, m = 10, seed = 123 )
nll_mid = lanczos_nll( pen, k = 10, m = 10, seed = 123 )

# Epsilon bias-correction estimator
sdrep = sdreport( obj, bias.correct = TRUE )

#
(nll_hi['nll'] - nll_mid['nll']) / (phat$eps)
#>      nll 
#> 4.733918
summary(sdrep)['sumexpu',]
#>            Estimate          Std. Error Est. (bias.correct) Std. (bias.correct) 
#>            4.064121            1.625808            4.734025                  NA
```

Runtime for this vignette: 9.7 secs

## Works cited
