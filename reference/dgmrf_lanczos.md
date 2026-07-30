# Approximate the log-density of a GMRF using Lanczos

Calculate log-density of a GMRF using Lanczos recursion

## Usage

``` r
dgmrf_lanczos(
  x,
  Q,
  k,
  mu = rep(0, length(x)),
  m = 3,
  log = FALSE,
  seed = NULL,
  orthogonalize = TRUE
)
```

## Arguments

- x:

  parameter vector used when calculating the Hessian matrix

- Q:

  Sparse precision matrix

- k:

  dimension for Kyrlov subspace

- mu:

  Mean parameter vector

- m:

  number of probe-vectors to use for approximating average and standard
  deviation of log-determinant

- log:

  Logical; Return log density?

- seed:

  if not NULL, then sets the seed. This is helfpul given that the
  Hutchinson probe vectors are randomly sampled, and comparisons have
  lower variance using a fixed seed.

- orthogonalize:

  Whether to do two-pass Gram-Schmidt re-normalization (much slower)

## Examples

``` r
library(RTMB)
library(Matrix)
library(lanczosRTMB)

# Create a precision for a 1D AR1 process
nx = 1e5
P = 0.5 * bandSparse( n = nx, k = c(-1,1) )
Q = (Diagonal(nx) - 0.5*t(P) ) %*% (Diagonal(nx) - 0.5*P)
x = RTMB:::rgmrf0(n = 1, Q )[,1]

# Exact log-density using RTMB
dgmrf(x, Q = Q, log = TRUE)
#> [1] -149139.8

# log-density using stochastic Lanczos quadrature
dgmrf_lanczos(x, Q = Q, k = 20, log = TRUE)
#> [1] -149074.3 -148934.0 -149245.2

# Log-density using incomplete Cholesky
# dgmrf_icpl(x, Q = Q, log = TRUE)
```
