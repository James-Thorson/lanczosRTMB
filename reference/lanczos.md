# Calculate Lanczos approximation

Calculate Q, alpha and beta for Lanzos approximation

## Usage

``` r
lanczos(Hq, q, k, x = attr(Hq, "env")$x0, orthogonalize = FALSE, tol = 1e-12)
```

## Arguments

- Hq:

  function that calculates the product `H %*% q` given probe `q` and
  parameters `x`

- q:

  initial vector used for defining a Kyrlov subspace

- k:

  dimension for Kyrlov subspace

- x:

  parameter vector used when calculating the Hessian matrix

- orthogonalize:

  Whether to do two-pass Gram-Schmidt re-normalization (much slower)

- tol:

  numerical tolerance for stopping algorithm given that no more terms
  are identifiable

## Examples

``` r
H = diag(exp(rnorm(5)))
q = rep(1,5)
Hq = function(q, x) (H %*% q)[,1]

L = lanczos(Hq, x = NULL, q, k = 5, ortho = TRUE)
T = tridiag(L$alpha, L$beta)

# Should match H if and only if L$m = nrow(H)
range(L$Q %*% T %*% t(L$Q) - H)
#> [1] -3.560819e-16  2.972869e-16
```
