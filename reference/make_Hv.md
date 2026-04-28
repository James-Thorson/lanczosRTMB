# Make function to calculate H %\*% v

Given a TMB object, make a function that efficiently calculates
`H %*% v` without constructing H itself, and instead using
`grad_x( grad_x(f) %** v)` given function f(x) that returns the negative
log-likelihood. This can then be used e.g. in Lanczos methods when H is
too large to construct explicitly

## Usage

``` r
make_Hv(obj, par = obj$env$last.par.best, tape)
```

## Arguments

- obj:

  TMB object (output from
  [`TMB::MakeADFun`](https://rdrr.io/pkg/TMB/man/MakeADFun.html))

- par:

  parameter vector `p` used when evaluating `H`

## Examples

``` r
u = rnorm(100)
y = rpois(length(u), exp(u))
nll = function(p) -1 * ( sum(dnorm(p$u,log=TRUE)) + sum(dpois(y,exp(p$u),log=TRUE)) )
obj = RTMB::MakeADFun( nll, list(u=u) )
Hv = make_Hv( obj )
# Confirm
v = rnorm( length(obj$par) )
all.equal( Hv(v)[1,], (obj$he()%*%v)[,1] )
#> [1] TRUE
```
