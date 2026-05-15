# Nonlinear minimizer using line search with approximated Newton solution

Nonlinear minimizer designed for cheap Hessian-vector products using
[make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md),
involving iterating a linear search along the truncated conjugate
gradient for a Newton.

## Usage

``` r
newton_CG(
  par,
  fn,
  gr,
  Hq,
  gr_tol = 1e-08,
  e_ratio = 0.01,
  maxit_newton = 100,
  maxit_CG = min(100, length(par)),
  c1 = 0.01,
  beta = 0.5,
  silent = FALSE
)
```

## Arguments

- par:

  initial parameter vector

- fn:

  function to evaluate negative log-likelihood

- gr:

  function to evaluate gradient of negative log-likelihood

- Hq:

  efficient Hessian-vector product function, e.g.,
  [make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md)

- gr_tol:

  early stopping condition for gradient of Newton solver

- e_ratio:

  early stopping condition for error of CG, relative to initial gradient

- maxit_newton:

  maximum iterations for Newton solver

- maxit_CG:

  maximum iterations for CG solution for each Newton iteration

- c1:

  stopping condition for line search given CG solution in each Newton
  iteration

- beta:

  updates in line search stepsize alpha when Armijo sufficient decrease
  condition fails

- silent:

  Be silent or print progress?

## Examples

``` r
library(Matrix)
library(RTMB)

# Settings
n = 10^4
rho = 0.99

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
y = x + 0.1 * rnorm(n)
which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
y[-which_seen] = NA

nll = function(p){
  -dgmrf(p$x, Q = Q, log = TRUE) - sum(dnorm(y, p$x, sd = 0.1, log=TRUE), na.rm=TRUE)
}
parlist = list( x=rnorm(n) )

tape = MakeTape(nll, parlist)
gr = tape$jacfun()
Hq = make_Hq( tape, x = unlist(parlist) )
H = gr$jacfun(sparse = TRUE)

start_time = Sys.time()
opt1 = optim(
  par = unlist(parlist),
  fn = tape,
  gr = gr,
  method = "L-BFGS-B",
  control = list(
    maxit = 1e4
  )
)
opt1$runtime = Sys.time() - start_time

opt2 = newton_CG(
  par = unlist(parlist),
  fn = tape,
  gr = gr,
  Hq = Hq
)
#> value: 49346.74 mgc: 35.80117 
#> value: 10369.95 mgc: 2.873467 
#> value: 8326.889 mgc: 0.2209165 
#> value: 8285.138 mgc: 0.02620617 
#> value: 8284.761 mgc: 0.002207074 
#> value: 8284.757 mgc: 0.0004615199 
#> value: 8284.757 mgc: 8.618456e-05 
#> value: 8284.757 mgc: 7.827714e-06 
#> value: 8284.757 mgc: 1.042457e-05 
#> value: 8284.757 mgc: 5.630421e-06 
#> value: 8284.757 mgc: 5.406651e-06 
#> value: 8284.757 mgc: 5.397407e-06 
#> value: 8284.757 mgc: 5.360276e-06 
#> value: 8284.757 mgc: 5.348572e-06 
#> value: 8284.757 mgc: 5.344506e-06 
#> Not improving much - will try early exit...PD hess?: TRUE 

# Compare the estimates and speed
matplot( cbind(opt1$par, opt2$par), type = "l", col = c("black","blue","red"), lty = "solid" )

c(opt1$runtime, opt2$runtime)
#> Time differences in secs
#> [1] 3.213687 1.103188
```
