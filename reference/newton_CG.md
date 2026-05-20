# Nonlinear minimizer using line search with approximated Newton solution

Nonlinear minimizer designed to use cheap Hessian-vector products using
[make_Hq](https://james-thorson.github.io/lanczosRTMB/reference/make_Hq.md).
The minimizer approximates a Newton update using truncated conjugate
gradient, while adapting a regularization parameter designed and using a
line-search for each step.

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
  line_steps = 20,
  smartsearch = TRUE,
  ustep = 1,
  power = 0.5,
  u0 = 1e-04,
  stop_if_nonPD = smartsearch,
  tol10 = 0.001,
  diagnostics = FALSE,
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
  iteration, where \\1\>c1\>0\\ seeks to prevent overshoot and resulting
  zigzagging.

- beta:

  updates in line search stepsize alpha when Armijo sufficient decrease
  condition fails

- line_steps:

  number of steps to explore for linear-search given Newton update

- smartsearch:

  Turn on adaptive stepsize algorithm for non-convex problems?

- ustep:

  Adaptive stepsize initial guess between 0 and 1.

- power:

  Parameter controlling adaptive stepsize.

- u0:

  Parameter controlling adaptive stepsize.

- stop_if_nonPD:

  whether to stop CG if recusion is not positive definite

- tol10:

  Try to exit if last 10 iterations not improved more than this.

- diagnostics:

  whether to provide extra diagnostics for each Newton iteration

- silent:

  Be silent or print progress?

## Details

This minimizer approximates Newton steps \\x\_{i+1} = x\_{i} -
H(x_i)^{-1} g(x_i)\\, where \\H(x_i) = \nabla^2 f(x_i)\\ is the Hessian
matrix and \\g(x_i) = \nabla f(x_i)\\ is the gradient of the negative
log-likelihood \\f(x_i)\\. It then uses several strategies for numerical
and computational efficiency.

1.  It approximates each Newton step \\H(x_i)^{-1} g(x_i)\\ using a
    truncated conjugate gradient algorithm that involves Hessian-vector
    products \\H(x_i) v\\. It then recursively improves the solution
    until a desired accuracy is reached, controlled by `gr_tol`, often
    using many fewer steps than the full dimension `length(x_i)`;

2.  It approximates each Newton solution \\H(x_i)^{-1} g(x_i)\\ without
    ever constructing or storing the Hessian \\H(x_i)\\, and instead
    computing Hessian-vector product \\H(x_i) v = \nabla( v^T \nabla
    f(x_i) )\\ using autodifferentiation in RTMB;

3.  After approximating each step direction \\H(x_i)^{-1} g(x_i)\\, it
    uses a line-search algorithm while decreasing the step-size if
    needed to find a suitable decrease in the objective function,
    controlled `line_steps`, `beta`, and `c1`;

4.  The CG is monitored to identify whether the Hessian matrix is
    positive definite (PD), and if it is not PD then the CG may be
    terminated, controlled by `stop_if_nonPD`.

5.  When using `smartsearch = TRUE`, if the Hessian is not PD, then a
    regularization `ustep` is decreased, corresponding to an increase in
    \\t=\frac{1}{u}-1\\ where the Newton step is actually using \\(H +
    tI)^{-1} g\\ where \\g\\ is the gradient, and \\1\>u\>0\\
    corresponds to \\t\>0\\. This increase or decrease in `ustep` (and
    associated decrease/increase in \\t\\) is copied from
    [`TMB::newton()`](https://rdrr.io/pkg/TMB/man/newton.html). If
    `smartsearch = FALSE` then \\ustep=1\\ and \\t=0\\ such that CG uses
    the Hessian corresponding to unregularized Newton steps. This
    smartsearch behavior controlled by `u0`, `ustep`, `power`, and
    `tol10`, and it corresponds to an adaptive "trust-region".

## Examples

``` r
library(Matrix)
library(RTMB)

# Settings
n = 10^3
rho = 0.99

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
y = rpois( n = n, lambda = exp(1 + x) )
which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
y[-which_seen] = NA

nll = function(p){
  -dgmrf(p$x, Q = Q, log = TRUE) - sum(dpois(y, lambda = exp(1 + p$x), log=TRUE), na.rm=TRUE)
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
  Hq = Hq,
  maxit_newton = 1e4,
  gr_tol = 1e-4
)
#> value: 89015524 mgc: 3219507 ustep: 1 
#> value: 64354055 mgc: 3213277 ustep: 1 
#> value: 31364467 mgc: 3174704 ustep: 1 
#> value: 12591638 mgc: 4594553 ustep: 1 
#> value: 5706788 mgc: 3706529 ustep: 1 
#> value: 2460777 mgc: 870265.4 ustep: 1 
#> value: 1459134 mgc: 866867.3 ustep: 1 
#> value: 856958.8 mgc: 655739.4 ustep: 1 
#> value: 573292.7 mgc: 341193.9 ustep: 1 
#> value: 358530 mgc: 330547.4 ustep: 1 
#> value: 250931.7 mgc: 35433.75 ustep: 1 
#> value: 98193.49 mgc: 34882.58 ustep: 1 
#> value: 47547.1 mgc: 30533.52 ustep: 1 
#> value: 17061.07 mgc: 15394.74 ustep: 1 
#> value: 13326.02 mgc: 716.7227 ustep: 1 
#> value: 7566.164 mgc: 666.5817 ustep: 1 
#> value: 4940.035 mgc: 604.0929 ustep: 1 
#> value: 2462.476 mgc: 195.9688 ustep: 1 
#> value: 1777.685 mgc: 732.34 ustep: 1 
#> value: 1741.084 mgc: 51.4548 ustep: 1 
#> value: 1588.798 mgc: 148.5927 ustep: 1 
#> value: 1507.604 mgc: 218.6208 ustep: 1 
#> value: 1495.131 mgc: 17.09987 ustep: 1 
#> value: 1419.398 mgc: 267.6423 ustep: 1 
#> value: 1410.093 mgc: 15.29279 ustep: 1 
#> value: 1377.204 mgc: 40.54135 ustep: 1 
#> value: 1354.708 mgc: 45.6425 ustep: 1 
#> value: 1338.097 mgc: 15.41863 ustep: 1 
#> value: 1326.49 mgc: 15.84385 ustep: 1 
#> value: 1317.276 mgc: 51.0363 ustep: 1 
#> value: 1309.429 mgc: 100.0655 ustep: 1 
#> value: 1308.314 mgc: 6.033478 ustep: 1 
#> value: 1301.54 mgc: 68.74179 ustep: 1 
#> value: 1300.979 mgc: 6.270615 ustep: 1 
#> value: 1294.885 mgc: 143.9958 ustep: 1 
#> value: 1294.701 mgc: 6.891263 ustep: 1 
#> value: 1290.016 mgc: 6.236139 ustep: 1 
#> value: 1285.92 mgc: 39.89157 ustep: 1 
#> value: 1282.313 mgc: 114.6718 ustep: 1 
#> value: 1282.283 mgc: 3.245747 ustep: 1 
#> value: 1279.039 mgc: 44.59548 ustep: 1 
#> value: 1276.891 mgc: 1.750034 ustep: 1 
#> value: 1270.842 mgc: 305.859 ustep: 1 
#> value: 1270.523 mgc: 22.47509 ustep: 1 
#> value: 1266.172 mgc: 1.307028 ustep: 1 
#> value: 1259.925 mgc: 352.0682 ustep: 1 
#> value: 1259.332 mgc: 20.07487 ustep: 1 
#> value: 1255.662 mgc: 1.1348 ustep: 1 
#> value: 1250.742 mgc: 950.9405 ustep: 1 
#> value: 1250.208 mgc: 14.05714 ustep: 1 
#> value: 1247.616 mgc: 0.7125134 ustep: 1 
#> value: 1246.375 mgc: 11.31505 ustep: 1 
#> value: 1245.379 mgc: 0.5172237 ustep: 1 
#> value: 1245.014 mgc: 24.60542 ustep: 1 
#> value: 1244.938 mgc: 0.923513 ustep: 1 
#> value: 1244.582 mgc: 31.97525 ustep: 1 
#> value: 1244.552 mgc: 2.48671 ustep: 1 
#> value: 1244.264 mgc: 3.13874 ustep: 1 
#> value: 1244.034 mgc: 2.000323 ustep: 1 
#> value: 1243.829 mgc: 46.0525 ustep: 1 
#> value: 1243.826 mgc: 1.040043 ustep: 1 
#> value: 1243.629 mgc: 4.291303 ustep: 1 
#> value: 1243.44 mgc: 3.019725 ustep: 1 
#> value: 1243.26 mgc: 33.94242 ustep: 1 
#> value: 1243.257 mgc: 1.326453 ustep: 1 
#> value: 1243.082 mgc: 7.285105 ustep: 1 
#> value: 1242.885 mgc: 142.7729 ustep: 1 
#> value: 1242.879 mgc: 13.64304 ustep: 1 
#> value: 1242.832 mgc: 0.8453799 ustep: 1 
#> value: 1241.523 mgc: 145.1961 ustep: 1 
#> value: 1241.039 mgc: 9.766808 ustep: 1 
#> value: 1240.703 mgc: 0.5281847 ustep: 1 
#> value: 1240.482 mgc: 6.910277 ustep: 1 
#> value: 1240.327 mgc: 0.310594 ustep: 1 
#> value: 1240.215 mgc: 51.70395 ustep: 1 
#> value: 1240.214 mgc: 1.47059 ustep: 1 
#> value: 1240.11 mgc: 0.3141597 ustep: 1 
#> value: 1239.898 mgc: 328.9395 ustep: 1 
#> value: 1239.881 mgc: 12.98487 ustep: 1 
#> value: 1239.805 mgc: 0.9094438 ustep: 1 
#> value: 1239.706 mgc: 4.895884 ustep: 1 
#> value: 1239.705 mgc: 0.2508401 ustep: 1 
#> value: 1239.578 mgc: 16.72867 ustep: 1 
#> value: 1239.539 mgc: 0.5207994 ustep: 1 
#> value: 1239.458 mgc: 9.958309 ustep: 1 
#> value: 1239.458 mgc: 0.1957865 ustep: 1 
#> value: 1239.383 mgc: 7.984085 ustep: 1 
#> value: 1239.318 mgc: 0.3262382 ustep: 1 
#> value: 1239.206 mgc: 10.48147 ustep: 1 
#> value: 1239.121 mgc: 0.500219 ustep: 1 
#> value: 1239.062 mgc: 2.822262 ustep: 1 
#> value: 1239.006 mgc: 8.921089 ustep: 1 
#> value: 1239.003 mgc: 0.2869496 ustep: 1 
#> value: 1238.923 mgc: 3.586376 ustep: 1 
#> value: 1238.856 mgc: 0.1792182 ustep: 1 
#> value: 1238.745 mgc: 14.6319 ustep: 1 
#> value: 1238.678 mgc: 0.7672271 ustep: 1 
#> value: 1238.628 mgc: 4.857017 ustep: 1 
#> value: 1238.614 mgc: 0.1401276 ustep: 1 
#> value: 1238.564 mgc: 1.445128 ustep: 1 
#> value: 1238.517 mgc: 0.8330959 ustep: 1 
#> value: 1238.481 mgc: 2.788303 ustep: 1 
#> value: 1238.447 mgc: 1.811714 ustep: 1 
#> value: 1238.419 mgc: 1.207904 ustep: 1 
#> value: 1238.392 mgc: 0.2471798 ustep: 1 
#> value: 1238.368 mgc: 20.3161 ustep: 1 
#> value: 1238.368 mgc: 1.455452 ustep: 1 
#> value: 1238.346 mgc: 25.38645 ustep: 1 
#> value: 1238.346 mgc: 1.135461 ustep: 1 
#> value: 1238.323 mgc: 1.478063 ustep: 1 
#> value: 1238.301 mgc: 0.5931894 ustep: 1 
#> value: 1238.28 mgc: 0.7219052 ustep: 1 
#> value: 1238.26 mgc: 9.597618 ustep: 1 
#> value: 1238.26 mgc: 0.1714888 ustep: 1 
#> value: 1238.24 mgc: 1.242776 ustep: 1 
#> value: 1238.221 mgc: 1.213026 ustep: 1 
#> value: 1238.203 mgc: 10.05505 ustep: 1 
#> value: 1238.203 mgc: 0.4123513 ustep: 1 
#> value: 1238.185 mgc: 2.133903 ustep: 1 
#> value: 1238.168 mgc: 6.480116 ustep: 1 
#> value: 1238.163 mgc: 0.33231 ustep: 1 
#> value: 1238.144 mgc: 2.524236 ustep: 1 
#> value: 1238.137 mgc: 0.1805898 ustep: 1 
#> value: 1237.934 mgc: 482.7996 ustep: 1 
#> value: 1237.887 mgc: 16.71825 ustep: 1 
#> value: 1237.839 mgc: 1.319747 ustep: 1 
#> value: 1237.804 mgc: 0.09301146 ustep: 1 
#> value: 1237.759 mgc: 137.6901 ustep: 1 
#> value: 1237.755 mgc: 5.073573 ustep: 1 
#> value: 1237.727 mgc: 0.2376452 ustep: 1 
#> value: 1237.719 mgc: 0.1747258 ustep: 1 
#> value: 1237.714 mgc: 0.2740638 ustep: 1 
#> value: 1237.71 mgc: 0.3948168 ustep: 1 
#> value: 1237.705 mgc: 0.4126858 ustep: 1 
#> value: 1237.7 mgc: 0.06912574 ustep: 1 
#> value: 1237.695 mgc: 15.6304 ustep: 1 
#> value: 1237.695 mgc: 0.1828836 ustep: 1 
#> value: 1237.69 mgc: 0.07442206 ustep: 1 
#> value: 1237.687 mgc: 0.5966552 ustep: 1 
#> value: 1237.683 mgc: 0.327073 ustep: 1 
#> value: 1237.679 mgc: 0.3303495 ustep: 1 
#> value: 1237.676 mgc: 4.061036 ustep: 1 
#> value: 1237.676 mgc: 0.3351669 ustep: 1 
#> value: 1237.672 mgc: 0.1755043 ustep: 1 
#> value: 1237.669 mgc: 0.1130465 ustep: 1 
#> value: 1237.666 mgc: 3.287955 ustep: 1 
#> value: 1237.665 mgc: 0.1522697 ustep: 1 
#> value: 1237.662 mgc: 0.1945496 ustep: 1 
#> value: 1237.659 mgc: 0.6089228 ustep: 1 
#> value: 1237.655 mgc: 2.455662 ustep: 1 
#> value: 1237.655 mgc: 0.1033828 ustep: 1 
#> value: 1237.652 mgc: 9.749426 ustep: 1 
#> value: 1237.652 mgc: 0.2001524 ustep: 1 
#> value: 1237.649 mgc: 0.212439 ustep: 1 
#> value: 1237.647 mgc: 10.75644 ustep: 1 
#> value: 1237.647 mgc: 0.3631673 ustep: 1 
#> value: 1237.644 mgc: 1.056668 ustep: 1 
#> value: 1237.644 mgc: 0.07192999 ustep: 1 
#> value: 1237.639 mgc: 1.014281 ustep: 1 
#> value: 1237.635 mgc: 0.5501944 ustep: 1 
#> value: 1237.631 mgc: 0.1705965 ustep: 1 
#> value: 1237.628 mgc: 2.816796 ustep: 1 
#> value: 1237.628 mgc: 0.07473023 ustep: 1 
#> value: 1237.624 mgc: 1.390641 ustep: 1 
#> value: 1237.621 mgc: 0.8808542 ustep: 1 
#> value: 1237.611 mgc: 16.40895 ustep: 1 
#> value: 1237.605 mgc: 0.2735365 ustep: 1 
#> value: 1237.602 mgc: 0.6491629 ustep: 1 
#> value: 1237.6 mgc: 0.5375266 ustep: 1 
#> value: 1237.597 mgc: 11.64255 ustep: 1 
#> value: 1237.597 mgc: 0.678111 ustep: 1 
#> value: 1237.595 mgc: 8.842523 ustep: 1 
#> value: 1237.595 mgc: 0.196632 ustep: 1 
#> value: 1237.593 mgc: 3.900048 ustep: 1 
#> value: 1237.593 mgc: 0.1633387 ustep: 1 
#> value: 1237.591 mgc: 0.3772736 ustep: 1 
#> value: 1237.589 mgc: 3.722599 ustep: 1 
#> value: 1237.589 mgc: 0.1990965 ustep: 1 
#> value: 1237.586 mgc: 19.105 ustep: 1 
#> value: 1237.586 mgc: 0.3272902 ustep: 1 
#> value: 1237.584 mgc: 0.5761234 ustep: 1 
#> value: 1237.582 mgc: 0.442342 ustep: 1 
#> value: 1237.58 mgc: 2.643903 ustep: 1 
#> value: 1237.58 mgc: 0.05364849 ustep: 1 
#> value: 1237.578 mgc: 0.5379474 ustep: 1 
#> value: 1237.576 mgc: 0.1266905 ustep: 1 
#> value: 1237.575 mgc: 0.5440149 ustep: 1 
#> value: 1237.573 mgc: 0.03105513 ustep: 1 
#> value: 1237.571 mgc: 0.2348746 ustep: 1 
#> value: 1237.569 mgc: 2.363607 ustep: 1 
#> value: 1237.569 mgc: 0.06348945 ustep: 1 
#> value: 1237.567 mgc: 1.378478 ustep: 1 
#> value: 1237.566 mgc: 0.05403464 ustep: 1 
#> value: 1237.564 mgc: 0.2433888 ustep: 1 
#> value: 1237.561 mgc: 0.6692679 ustep: 1 
#> value: 1237.56 mgc: 2.66106 ustep: 1 
#> value: 1237.56 mgc: 0.1612971 ustep: 1 
#> value: 1237.558 mgc: 4.977451 ustep: 1 
#> value: 1237.558 mgc: 0.1345789 ustep: 1 
#> value: 1237.557 mgc: 0.5436482 ustep: 1 
#> value: 1237.556 mgc: 5.608521 ustep: 1 
#> value: 1237.556 mgc: 0.1813347 ustep: 1 
#> value: 1237.555 mgc: 0.2526272 ustep: 1 
#> value: 1237.553 mgc: 0.1751241 ustep: 1 
#> value: 1237.552 mgc: 0.3394202 ustep: 1 
#> value: 1237.551 mgc: 1.345103 ustep: 1 
#> value: 1237.551 mgc: 0.05595323 ustep: 1 
#> value: 1237.55 mgc: 4.03254 ustep: 1 
#> value: 1237.55 mgc: 0.09771109 ustep: 1 
#> value: 1237.548 mgc: 9.295019 ustep: 1 
#> value: 1237.548 mgc: 0.3343939 ustep: 1 
#> value: 1237.545 mgc: 1.316269 ustep: 1 
#> value: 1237.545 mgc: 0.07793788 ustep: 1 
#> value: 1237.544 mgc: 22.87973 ustep: 1 
#> value: 1237.544 mgc: 1.0467 ustep: 1 
#> value: 1237.543 mgc: 0.04498305 ustep: 1 
#> value: 1237.511 mgc: 8.379994 ustep: 1 
#> value: 1237.504 mgc: 0.5018816 ustep: 1 
#> value: 1237.494 mgc: 0.03038277 ustep: 1 
#> value: 1237.492 mgc: 1.379571 ustep: 1 
#> value: 1237.49 mgc: 0.03577944 ustep: 1 
#> value: 1237.489 mgc: 0.1596464 ustep: 1 
#> value: 1237.488 mgc: 0.658749 ustep: 1 
#> value: 1237.488 mgc: 0.03858126 ustep: 1 
#> value: 1237.487 mgc: 6.332127 ustep: 1 
#> value: 1237.487 mgc: 0.4437963 ustep: 1 
#> value: 1237.487 mgc: 0.01956423 ustep: 1 
#> value: 1237.486 mgc: 0.7696525 ustep: 1 
#> value: 1237.485 mgc: 0.04056552 ustep: 1 
#> value: 1237.484 mgc: 1.738066 ustep: 1 
#> value: 1237.484 mgc: 0.1443213 ustep: 1 
#> value: 1237.484 mgc: 0.1646361 ustep: 1 
#> value: 1237.484 mgc: 0.06157323 ustep: 1 
#> value: 1237.484 mgc: 2.765803 ustep: 1 
#> value: 1237.484 mgc: 0.1022818 ustep: 1 
#> value: 1237.483 mgc: 0.6241182 ustep: 1 
#> value: 1237.483 mgc: 0.02595786 ustep: 1 
#> value: 1237.483 mgc: 18.56712 ustep: 1 
#> value: 1237.483 mgc: 0.5292831 ustep: 1 
#> value: 1237.483 mgc: 0.02103776 ustep: 1 
#> value: 1237.481 mgc: 0.6827625 ustep: 1 
#> value: 1237.48 mgc: 0.01685679 ustep: 1 
#> value: 1237.479 mgc: 1.252195 ustep: 1 
#> value: 1237.479 mgc: 0.0406483 ustep: 1 
#> value: 1237.479 mgc: 2.899194 ustep: 1 
#> value: 1237.479 mgc: 0.1344043 ustep: 1 
#> value: 1237.479 mgc: 1.433738 ustep: 1 
#> value: 1237.479 mgc: 0.06481301 ustep: 1 
#> value: 1237.479 mgc: 0.06631906 ustep: 1 
#> value: 1237.478 mgc: 1.178681 ustep: 1 
#> value: 1237.478 mgc: 0.07965458 ustep: 1 
#> value: 1237.478 mgc: 0.03013213 ustep: 1 
#> value: 1237.478 mgc: 0.0499719 ustep: 1 
#> value: 1237.478 mgc: 0.1256377 ustep: 1 
#> value: 1237.477 mgc: 3.212127 ustep: 1 
#> value: 1237.477 mgc: 0.02588252 ustep: 1 
#> value: 1237.477 mgc: 1.659809 ustep: 1 
#> value: 1237.477 mgc: 0.1404657 ustep: 1 
#> value: 1237.477 mgc: 0.1001358 ustep: 1 
#> value: 1237.477 mgc: 0.426311 ustep: 1 
#> value: 1237.477 mgc: 0.01659373 ustep: 1 
#> value: 1237.476 mgc: 4.971518 ustep: 1 
#> value: 1237.476 mgc: 0.3143253 ustep: 1 
#> value: 1237.476 mgc: 0.3314531 ustep: 1 
#> value: 1237.476 mgc: 0.02165891 ustep: 1 
#> value: 1237.476 mgc: 0.3892353 ustep: 1 
#> value: 1237.476 mgc: 0.02270085 ustep: 1 
#> value: 1237.475 mgc: 0.2337763 ustep: 1 
#> value: 1237.475 mgc: 0.07697692 ustep: 1 
#> value: 1237.474 mgc: 0.6950607 ustep: 1 
#> value: 1237.474 mgc: 0.06229591 ustep: 1 
#> value: 1237.474 mgc: 0.1647303 ustep: 1 
#> value: 1237.473 mgc: 0.7798811 ustep: 1 
#> value: 1237.473 mgc: 0.04839756 ustep: 1 
#> value: 1237.473 mgc: 2.633969 ustep: 1 
#> value: 1237.473 mgc: 0.2321419 ustep: 1 
#> value: 1237.473 mgc: 0.8728963 ustep: 1 
#> value: 1237.473 mgc: 0.04051442 ustep: 1 
#> value: 1237.472 mgc: 1.131023 ustep: 1 
#> value: 1237.472 mgc: 0.02170852 ustep: 1 
#> value: 1237.472 mgc: 0.2231234 ustep: 1 
#> value: 1237.472 mgc: 0.08656786 ustep: 1 
#> value: 1237.472 mgc: 0.2611268 ustep: 1 
#> value: 1237.471 mgc: 0.01607679 ustep: 1 
#> value: 1237.47 mgc: 0.1662974 ustep: 1 
#> value: 1237.47 mgc: 0.01324078 ustep: 1 
#> value: 1237.469 mgc: 7.279448 ustep: 1 
#> value: 1237.469 mgc: 0.07121008 ustep: 1 
#> value: 1237.468 mgc: 0.01839549 ustep: 1 
#> value: 1237.468 mgc: 5.311373 ustep: 1 
#> value: 1237.468 mgc: 0.1423479 ustep: 1 
#> value: 1237.467 mgc: 0.21836 ustep: 1 
#> value: 1237.466 mgc: 0.3738214 ustep: 1 
#> value: 1237.466 mgc: 0.009535663 ustep: 1 
#> value: 1237.465 mgc: 0.2464702 ustep: 1 
#> value: 1237.464 mgc: 0.01377742 ustep: 1 
#> value: 1237.464 mgc: 0.1923752 ustep: 1 
#> value: 1237.464 mgc: 1.243238 ustep: 1 
#> value: 1237.464 mgc: 0.02498514 ustep: 1 
#> value: 1237.464 mgc: 0.1066386 ustep: 1 
#> value: 1237.463 mgc: 0.06840805 ustep: 1 
#> value: 1237.463 mgc: 2.947593 ustep: 1 
#> value: 1237.463 mgc: 0.06942906 ustep: 1 
#> value: 1237.463 mgc: 0.2452939 ustep: 1 
#> value: 1237.463 mgc: 0.01134828 ustep: 1 
#> Not improving much - will try early exit...PD hess?: TRUE 

# Compare the estimates and speed
matplot( cbind(opt1$par, opt2$par), type = "l", col = c("black","blue","red"), lty = "solid" )

c(opt1$runtime, opt2$runtime)
#> Time differences in secs
#> [1]  3.488409 19.111351
```
