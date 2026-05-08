# Truncated conjugate gradient using Hessian-vector products

Approximate the soution to a linear system \\Hq = b\\ using truncated
conjugate gradient, without forming the Hessian matrix directly and
instead using Hessian-vector products e.g., using forwards-on-reverse
autodiff Adapted from `mcmcsae::CG` under GPL-3 licence to use
Hessian-vector products

## Usage

``` r
CG(
  b,
  Hq,
  x = 0 * b,
  Minv = Diagonal(n = length(b)),
  max.it = length(b),
  e = 1e-10,
  silent = TRUE
)
```

## Arguments

- b:

  vector to solve for

- Hq:

  efficient Hessian-vector product function, e.g.,
  [make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md)

- x:

  initial guess for solution

- Minv:

  preconditioner matrix (uses Identity by default)

- max.it:

  maximum iterations (often should be less than `length(b)`)

- e:

  error criterion

- silent:

  Be silent or print progress?

## Examples

``` r
# Settings
n = 10^4
rho = 0.99

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
#> Error in bandSparse(n = n, k = c(-1), diagonals = list(rep(1, n))): could not find function "bandSparse"
Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
#> Error in Diagonal(n): could not find function "Diagonal"
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'A' in selecting a method for function 'Cholesky': object 'Q' not found
y = x + 0.1 * rnorm(n)
#> Error: object 'x' not found
which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
y[-which_seen] = NA
#> Error: object 'y' not found

nll = function(p){
  -dgmrf(p$x, Q = Q, log = TRUE) - sum(dnorm(y, p$x, sd = 0.1, log=TRUE), na.rm=TRUE)
}
parlist = list( x=rnorm(n) )

tape = MakeTape(nll, parlist)
#> Error in f(x): object 'Q' not found
gr = tape$jacfun()
#> Error: object 'tape' not found
Hq = make_Hq( tape, x = unlist(parlist) )
#> Error: object 'tape' not found
H = gr$jacfun(sparse = TRUE)
#> Error: object 'gr' not found

#
b = gr(unlist(parlist))[1,]
#> Error in gr(unlist(parlist)): could not find function "gr"
x = solve( H(unlist(parlist)), b)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'a' in selecting a method for function 'solve': could not find function "H"
out = CG(
  b = b,
  Hq = Hq
)
#> Error: object 'b' not found

#
plot( x, out$x )
#> Error: object 'x' not found
```
