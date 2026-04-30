# Make function to calculate H %\*% q

Given a TMB object, make a function that efficiently calculates
`H %*% q` without constructing H itself, and instead using
`grad_u( grad_u(f) %** q)` given function f(x) that returns the negative
log-likelihood given `x = u` with fixed `v`. This can then be used e.g.
in Lanczos methods when H is too large to construct explicitly

## Usage

``` r
make_Hq(tape, x)
```

## Arguments

- tape:

  Alternative to specifying `obj`

- x:

  parameter vector `x` used when evaluating `H`

## Details

The output `Hq = make_Hq( tape, x )` takes as argument a probe
\\\mathbf{q}\\ and outputs \\\mathbf{Hq}\\. To change the point at which
\\\mathbf{Hq}\\ is evaluated, assign a new value to `attr(Hq,"env")$x`.
RTMB then does a `force.update()` to update the tape based on that new
value.

## Examples

``` r
u = rnorm(100)
y = rpois(length(u), exp(u))
nll = function(p) -1 * ( sum(dnorm(p$u,log=TRUE)) + sum(dpois(y,exp(p$u),log=TRUE)) )
obj = RTMB::MakeADFun( nll, list(u=u), silent = TRUE )
Hq = make_Hq( GetTape(obj), obj$par )
# Confirm
q = rnorm( length(obj$par) )
all.equal( Hq(q)[1,], (obj$he()%*%q)[,1] )
#> [1] TRUE
```
