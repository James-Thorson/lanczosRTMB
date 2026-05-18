
library(Matrix)
library(RTMB)
library(lanczosRTMB)

# Settings
set.seed(123)
nx = 10^2
ny = 10^2
rho = 0.9

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
Px = bandSparse( n = nx, k = c(-1,1), diagonals = list(rep(0.5,nx),rep(0.5,nx)) )
Py = bandSparse( n = ny, k = c(-1,1), diagonals = list(rep(0.5,ny),rep(0.5,ny)) )
Ix = Diagonal(n=nx)
Iy = Diagonal(n=ny)
Qx = (Ix - rho * t(Px)) %*% (10 * Ix) %*% (Ix - rho * Px)
Qy = (Iy - rho * t(Py)) %*% (10 * Iy) %*% (Iy - rho * Py)
Q = kronecker( Qy, Qx )
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
y = rpois( nx*ny, lambda = exp(1 + x) )
which_seen = sample( seq_len(nx*ny), size = nx*ny/100, replace = FALSE)
sumy = sum(y)
y[-which_seen] = NA

nll = function(p){
  Qx = (Ix - plogis(p$invlogis_rho) * t(Px)) %*% (exp(2*p$logtau) * Ix) %*% (Ix - plogis(p$invlogis_rho) * Px)
  Qy = (Iy - plogis(p$invlogis_rho) * t(Py)) %*% (exp(2*p$logtau) * Iy) %*% (Iy - plogis(p$invlogis_rho) * Py)
  Q = kronecker( Qy, Qx )
  -dgmrf(p$x, Q = Q, log = TRUE) -
  sum(dpois(y, exp(p$mu + p$x), log=TRUE), na.rm=TRUE) -
  dnorm(log(sumy), log(sum(exp(p$mu + p$x))), sd = 0.01, log = TRUE)
}
parlist = list( x=rnorm(nx*ny), logtau = 0, invlogis_rho = 0, mu = 0 )

start_time = Sys.time()
obj = lanczos_MakeADFun(
  nll,
  parlist,
  random = "x",
  k = 100,
  #method = "L-BFGS-B",
  make_gr = FALSE,
  silent = TRUE
)
tmp = obj$fn( obj$par, gr_tol = 1e-8, smartsearch = TRUE, what = "all" )
opt = optim(
  obj$par,
  \(x) obj$fn(x, gr_tol = 1e-8, smartsearch = TRUE ),
  #\(par) obj$fn(par, gr_tol = 1e-06, e_ratio = 0.1, maxit_newton = 200 , maxit_CG = 200 ),
  control = list(trace = 6),
  method = "L-BFGS-B"
)
opt$run_time = Sys.time() - start_time

start_time = Sys.time()
obj2 = MakeADFun(
  nll,
  parlist,
  random = "x",
  silent = TRUE
)
opt2 = nlminb( obj2$par, obj2$fn, control = list(trace = 1) )
opt2$run_time = Sys.time() - start_time

