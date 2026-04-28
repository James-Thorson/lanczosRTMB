# Assemble tri-diagonal matrix

Assemble tri-diagonal matrix from alpha and beta from Lanczos method

## Usage

``` r
lanczos(Hv, q1, k, orthogonalize = FALSE, tol = 1e-12)
```

## Arguments

- Hv:

  function that calculates the product `Hv`

- q1:

  initial vector used for defining a Kyrlov subspace

- k:

  dimension for Kyrlov subspace

- orthogonalize:

  Whether to do two-pass Gram-Schmidt re-normalization (much slower)

- tol:

  numerical tolerance for stopping algorithm given that no more terms
  are identifiable

## Examples

``` r
H = diag(exp(rnorm(5)))
v = rep(1,5)
Hv = function(v) (H %*% v)[,1]

L = lanczos(Hv, v, k = 5, ortho = TRUE)
T = tridiag(L$alpha, L$beta)

# Should match H if and only if L$m = nrow(H)
range(L$Q %*% T %*% t(L$Q) - H)
#> [1] -4.440892e-16  4.440892e-16
```
