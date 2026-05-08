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
