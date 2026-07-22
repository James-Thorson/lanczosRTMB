# Lanczos methods for penalized and marginal likelihood models

This vignette demonstrates how Lanczos methods can be used to
approximate standard errors and epsilon bias-correction for penalized
likelihood models, without directly forming or inverting the Hessian
matrix representing associations among random effects.

``` r

library(Matrix)
library(RTMB)
library(lanczosRTMB)
library(numDeriv)
library(memprof)
```

## Generalized linear mixed model

We first demonstrate Lanczos methods using a simple generalized linear
mixed model.

### Simulate from a GLMM

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

### Fit as marginal likelihood without Cholesky decomposition

We first show that we can re-fit the model using the marginal
likelihood, evaluated using Hutchinson-Lanczos instead of the Cholesky
factorization that is default in TMB:

``` r

parlist$eps = numeric()
pen0 = lanczos_MakeADFun(
  nll, 
  parlist,
  random = "u",
  k = 10,
  silent = TRUE
)

# optimize 
opt_pen0 = nlminb( pen0$par, pen0$fn )
opt_pen0
#> $par
#>         mu      logsd      logcv 
#> -0.1084795 -0.1284993 -0.6074894 
#> 
#> $objective
#> [1] 36.46912
#> 
#> $convergence
#> [1] 1
#> 
#> $iterations
#> [1] 11
#> 
#> $evaluations
#> function gradient 
#>       31       35 
#> 
#> $message
#> [1] "false convergence (8)"
```

### Fit as a penalized likelihood model

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
#> value: 41.78943 mgc: 9.836576e-13 ustep: 0.9999
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

### Delta methods

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

### Epsilon methods

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

## Spatial change of support

Next, we fit a model that involves spatial change of support, which
breaks sparsity in the inner Hessian. We again simulate data from the
change-in-support process:

``` r

# Settings
set.seed(123)
nx = 25
ny = 100
rho = 0.9

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
Px = bandSparse( n = nx, k = c(-1,1), diagonals = list(rep(0.5,nx),rep(0.5,nx)) )
Py = bandSparse( n = ny, k = c(-1,1), diagonals = list(rep(0.5,ny),rep(0.5,ny)) )
Ix = Diagonal(n=nx)
Iy = Diagonal(n=ny)
Qx = (Ix - rho * t(Px)) %*% (10 * Ix) %*% (Ix - rho * Px)
Qy = (Iy - rho * t(Py)) %*% (10 * Iy) %*% (Iy - rho * Py)
Q = kronecker( Qy, Qx )
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
y = rpois( nx*ny, lambda = exp(1 + x) )
which_seen = sample( seq_len(nx*ny), size = nx*ny/100, replace = FALSE)
sumy = sum(y)
#y[-which_seen] = NA
```

We then define the joint likelihood that includes a term for the total
across samples:

``` r

nll = function(p){
  Qx = (Ix - plogis(p$invlogis_rho) * t(Px)) %*% (exp(2*p$logtau) * Ix) %*% (Ix - plogis(p$invlogis_rho) * Px)
  Qy = (Iy - plogis(p$invlogis_rho) * t(Py)) %*% (exp(2*p$logtau) * Iy) %*% (Iy - plogis(p$invlogis_rho) * Py)
  Q = kronecker( Qy, Qx )
  loglik1 = dgmrf(p$x, Q = Q, log = TRUE)
  loglik2 = sum(dpois(y, exp(p$mu + p$x), log=TRUE), na.rm=TRUE)
  loglik3 = dnorm(log(sumy), log(sum(exp(p$mu + p$x))), sd = 0.01, log = TRUE)
  -1 * ( loglik1 + loglik2 + loglik3 )
}
parlist = list( x=rnorm(nx*ny), logtau = 0, invlogis_rho = 0, mu = 0 )
```

We then fit this as normal in RTMB:

``` r

start_time = Sys.time()
fit_RTMB = with_monitor({
  obj_RTMB = MakeADFun(
    nll,
    parlist,
    random = "x",
    silent = TRUE
  )
  nlminb(
    obj_RTMB$par,
    obj_RTMB$fn,
    obj_RTMB$gr,
    control = list(trace = 1)
  )
})
#>   0:     6695.9923:  0.00000  0.00000  0.00000
#>   1:     6553.9980: 0.0994045 0.00420131 0.0100547
#>   2:     6480.2008: 0.192113 0.0314836 0.0357614
#>   3:     6090.0149: 0.384639 0.664367 0.485641
#>   4:     6050.5174: 0.567800 0.677494 0.515123
#>   5:     5976.8525: 0.532551 0.847680 0.581331
#>   6:     5893.1189: 0.743109  1.11859 0.724962
#>   7:     5874.6455: 0.809255  1.61679 0.176454
#>   8:     5843.9260: 0.756846  1.35265 0.869935
#>   9:     5837.0556: 0.828967  1.36482 0.884478
#>  10:     5828.5730: 0.812307  1.43648 0.896698
#>  11:     5814.9003: 0.898049  1.55851 0.895360
#>  12:     5795.1786: 0.998974  1.83899 0.906668
#>  13:     5785.5759:  1.08507  2.02378  1.12443
#>  14:     5783.3742:  1.15799  2.24198  1.08329
#>  15:     5783.1478:  1.15832  2.24106  1.14629
#>  16:     5783.1369:  1.15702  2.24079  1.15827
#>  17:     5783.1365:  1.15674  2.24118  1.15868
#>  18:     5783.1365:  1.15671  2.24123  1.15851
#>  19:     5783.1365:  1.15671  2.24123  1.15850
fit_RTMB$run_time = Sys.time() - start_time
```

We can see that the Hessian is dense:

``` r

Matrix::image( obj_RTMB$env$spHess(random=TRUE) )
```

![](lanczos_files/figure-html/unnamed-chunk-17-1.png)

Alternatively, we can fit this using Lanczos methods, which are less
sensitive to the lack-of-sparsity. Here, we avoid constructing the
gradient of the Laplace-Lanczos approximation to the log-marginal
likelihood, pending improvements:

``` r

start_time = Sys.time()
fit_lanczos = with_monitor({
  obj = lanczos_MakeADFun(
    nll,
    parlist,
    random = "x",
    k = 40,
    make_gr = TRUE,
    silent = TRUE
  )
  nlminb(
    obj$par,
    obj$fn,
    \(x) obj$gr(x, method = "simple"),
    control = list(trace = 1)
  )
})
#>   0:     6700.4179:  0.00000  0.00000  0.00000
#>   1:     6558.9301: 0.0994048 0.00408596 0.0100989
#>   2:     6485.6350: 0.192093 0.0312778 0.0359720
#>   3:     6102.0546: 0.379578 0.663808 0.488471
#>   4:     6064.7262: 0.572814 0.673877 0.517953
#>   5:     5986.5229: 0.524891 0.807029 0.653174
#>   6:     5914.1339: 0.743371  1.13174 0.644615
#>   7:     5906.5239: 0.793666  1.65720 0.0663996
#>   8:     5863.2962: 0.762242  1.34634 0.784280
#>   9:     5857.2336: 0.835181  1.37681 0.813692
#>  10:     5845.8018: 0.821016  1.43142 0.876393
#>  11:     5821.7346: 0.986174  1.73194 0.968901
#>  12:     5820.3569:  1.12752  2.20733 0.652802
#>  13:     5809.9740:  1.13352  2.18266  1.24038
#>  14:     5809.2796:  1.10236  2.11329  1.16778
#>  15:     5809.2561:  1.09407  2.08951  1.16689
#>  16:     5809.2553:  1.09517  2.09205  1.16495
#>  17:     5809.2553:  1.09529  2.09214  1.16503
#>  18:     5809.2553:  1.09529  2.09210  1.16502
fit_lanczos$run_time = Sys.time() - start_time
```

Where we can compare the runtime for the two models:

``` r

runtime = c(
  RTMB = fit_RTMB$run_time,
  Lanczos = fit_lanczos$run_time
)
knitr::kable( runtime, digits=2, caption="Run-times" )
```

|         | x           |
|:--------|:------------|
| RTMB    | 49.29 secs  |
| Lanczos | 134.90 secs |

Run-times {.table}

And we can also compare memory use:

``` r

par(mfrow = c(2,1) )
mem1 = memprof:::used_memory_total_by_time(fit_RTMB$memory_use)
mem1[,'used'] = mem1[,'used'] - min(mem1[,'used'])
plot( mem1, type = "l", main = "RTMB")
mem2 = memprof:::used_memory_total_by_time(fit_lanczos$memory_use)
mem2[,'used'] = mem2[,'used'] - min(mem2[,'used'])
plot( mem2, type = "l", main = "Lanczos")
```

![](lanczos_files/figure-html/unnamed-chunk-20-1.png)

which shows that Lanczos uses similar memory:

``` r

max_memory = c(
  RTMB = max(mem1[,'used']) / 1e9,
  Lanczos = max(mem2[,'used']) / 1e9
)
knitr::kable( max_memory, digits=2, caption="Maximum memory use (GB)" )
```

|         |    x |
|:--------|-----:|
| RTMB    | 1.49 |
| Lanczos | 1.58 |

Maximum memory use (GB) {.table}

Runtime for this vignette: 3.74 mins

### Works cited
