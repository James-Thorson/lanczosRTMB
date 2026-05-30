# Make function to calculate H %\*% q

Given a TMB object, make a function that efficiently calculates a
Hessian-vector product (HVP) `H %*% q`.

## Usage

``` r
make_Hq(
  tape,
  x0,
  which_random = seq_along(x0),
  method = c("reverse-on-reverse", "sparse", "FD-on-reverse")
)
```

## Arguments

- tape:

  Alternative to specifying `obj`

- x0:

  parameter vector `x` (or list coersed to vector) used as default when
  evaluating `H`

- which_random:

  integer-vector indicating which elements of `x` correspond to random
  effects, where the probe `q` then has length `length(which_random)`

- method:

  See details

## Value

A function with two arguments:

- q a probe vector with length `length(which_random)`

- a vector with length `length(x0)`, where the Hessian is calculated

## Details

This can then be used e.g. in Lanczos methods when H is too large to
construct explicitly.

When `method = "reverse-on-reverse"`, `make_Hq` calculates a HVP without
constructing H itself, and instead using `grad_u( grad_u(f) %** q)`
given function f(x) that returns the negative log-likelihood given
`x = u` with fixed `v`

When `method = "sparse"`, `make_Hq` instead calculates a HVP by
calculating and storing the sparse Hessian in a local environment. The
resulting function can be used with `update_H = FALSE` to use the
pre-calculated Hessian as-is, or `update_H = TRUE` to recalculate the
sparse Hessian, store the update in the local environment and then
calculate the HVP. `update_H = FALSE` is then useful when repeadly using
the same Hessian in a HVP.

When `method = "FD-on-verse"`, `make_Hq` instead calculates a two-sided
finite-difference approximation to a forward-on-reverse autodiff
calculation, using `delta = 1e-6` forward-on-reverse (and FD of autodiff
gradients) is efficient given that `grad_u(f) %** q` has length of one.

`qprime` is defined internally where `qprime[which_random] = q` and
`qprime[!which_random] = 0`, where `length(qprime)` is equal to
`length(x)`

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
nll = function(p){
  nll1 = dnorm(p$u, mean=p$mu, sd=exp(p$logsd), log=TRUE)
  nll2 = dgamma(y, shape = 1/exp(2*p$logcv), scale = exp(p$u) * exp(2*p$logcv), log=TRUE)
  jnll = -1 * ( sum(nll1) + sum(nll2) )
  return(jnll)
}

# Build with RTMB
params = list(u=u, mu = 0, logsd = 0, logcv = 0)
obj = MakeADFun( nll, params, random = "u", silent = TRUE )

# Build with with bespoke function
tape = MakeTape( nll, params )
which_random = seq_len(n)
Hq = make_Hq(
  tape,
  x = params,
  which_random = which_random
)

# Compare them
q = rnorm(length(which_random))
all.equal(
  Hq(q),
  (obj$env$spHess(par = unlist(params), random=TRUE) %*% q)[,1]
)
#> [1] TRUE

# Compare them when passing new value
x_new = unlist(params)
  x_new['logsd'] = 1
all.equal(
  Hq(q, x_new),
  (obj$env$spHess(par = x_new, random=TRUE) %*% q)[,1]
)
#> [1] TRUE

# Compare speed with explicit sparse H
Hq2 = make_Hq(
  tape,
  x = params,
  which_random = which_random,
  method = "sparse"
)
x_new[which_random] = rnorm(length(which_random))

# Compare speed with finite-difference approximation to forward-on-reverse
Hq3 = make_Hq(
  tape,
  x = params,
  which_random = which_random,
  method = "FD-on-reverse"
)
x_new[which_random] = rnorm(length(which_random))

system.time(Hq(q, x_new))
#>    user  system elapsed 
#>       0       0       0 
system.time(Hq2(q, x_new))
#>    user  system elapsed 
#>   0.006   0.000   0.006 
system.time(Hq2(q, x_new, update_H = FALSE))
#>    user  system elapsed 
#>       0       0       0 
```
