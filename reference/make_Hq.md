# Make function to calculate H %\*% q

Given a TMB object, make a function that efficiently calculates
`H %*% q` without constructing H itself, and instead using
`grad_u( grad_u(f) %** q)` given function f(x) that returns the negative
log-likelihood given `x = u` with fixed `v`. This can then be used e.g.
in Lanczos methods when H is too large to construct explicitly

## Usage

``` r
make_Hq(tape, x0, which_random = seq_along(x0))
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

## Value

A function with two arguments:

- q a probe vector with length `length(which_random)`

- a vector with length `length(x0)`, where the Hessian is calculated

## Details

The output `Hq = make_Hq( tape, x )` takes as argument a probe
\\\mathbf{q}\\ and outputs \\\mathbf{Hq}\\. To change the point at which
\\\mathbf{Hq}\\ is evaluated, assign a new value to `attr(Hq,"env")$x`.
RTMB then does a `force.update()` to update the tape based on that new
value.

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

# Build with Lanczos
tape = MakeTape( nll, params )
which_random = 1:30
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
```
