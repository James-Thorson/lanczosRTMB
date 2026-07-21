
library(Matrix)
library(RTMB)
library(lanczosRTMB)
library(numDeriv)

# Settings
set.seed(123)
nx = 10^1
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
#sumy = NA
sumy = sum(y)
#y[-which_seen] = NA

nll = function(p){
  Qx = (Ix - plogis(p$invlogis_rho) * t(Px)) %*% (exp(2*p$logtau) * Ix) %*% (Ix - plogis(p$invlogis_rho) * Px)
  Qy = (Iy - plogis(p$invlogis_rho) * t(Py)) %*% (exp(2*p$logtau) * Iy) %*% (Iy - plogis(p$invlogis_rho) * Py)
  Q = kronecker( Qy, Qx )
  loglik1 = dgmrf(p$x, Q = Q, log = TRUE)
  loglik2 = sum(dpois(y, exp(p$mu + p$x), log=TRUE), na.rm=TRUE)
  loglik3 = sum(dnorm(log(sumy), log(sum(exp(p$mu + p$x))), sd = 0.01, log = TRUE), na.rm=TRUE)
  -1 * ( loglik1 + loglik2 + loglik3 )
}
parlist = list( x=rnorm(nx*ny), logtau = 0, invlogis_rho = 0, mu = 0 )

################
# Compare initial optimizers
################

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
x = tmp$inner_opt$par

obj2 = MakeADFun(
  nll,
  parlist,
  random = "x",
  silent = TRUE
)
tmp2 = obj2$fn( obj2$par )
x2 = obj2$env$last.par[obj2$env$lrandom()]

cor(x, x2)
c(tmp$nll, tmp2)

################
# Check gradients
################

gr = grad( obj$fn, obj$par )
gr2 = obj2$gr( obj$par )
cbind( obj$par, gr, gr2 )

################
# Compare MLE and speeds
################

start_time = Sys.time()
obj = lanczos_MakeADFun(
  nll,
  parlist,
  random = "x",
  k = 40,
  #method = "L-BFGS-B",
  make_gr = FALSE,
  silent = TRUE
)
#opt = optim(
#  obj$par,
#  \(x) obj$fn(x, gr_tol = 1e-8, smartsearch = TRUE ),
#  #\(par) obj$fn(par, gr_tol = 1e-06, e_ratio = 0.1, maxit_newton = 200 , maxit_CG = 200 ),
#  control = list(trace = 6),
#  method = "L-BFGS-B"
#)
opt = nlminb(
  obj$par,
  \(x) obj$fn(x), # , gr_tol = 1e-8, smartsearch = TRUE ),
  #\(par) obj$fn(par, gr_tol = 1e-06, e_ratio = 0.1, maxit_newton = 200 , maxit_CG = 200 ),
  control = list(trace = 1)
)
opt$run_time = Sys.time() - start_time

# Try with gradient
start_time = Sys.time()
obj = lanczos_MakeADFun(
  nll,
  parlist,
  random = "x",
  k = 40,
  make_gr = TRUE,
  silent = TRUE,
  pu_update = "implicit"
)
opt = nlminb(
  obj$par,
  \(x) obj$fn(x, gr_tol = 1e-8, smartsearch = TRUE ),
  \(x) obj$gr(x, method = "simple", method.args = list(eps = 1e-8)),  # fixed_Q = FALSE,
  control = list(trace = 1)
)
opt$run_time = Sys.time() - start_time

# Try with FD gradient
start_time = Sys.time()
obj = lanczos_MakeADFun(
  nll,
  parlist,
  random = "x",
  k = 40,
  make_gr = FALSE,
  silent = TRUE
)
opt = nlminb(
  obj$par,
  \(x) obj$fn(x, gr_tol = 1e-8, smartsearch = TRUE ),
  \(x) grad(obj$fn, x), # , method = "Richardson", method.args = list(eps = 1e-4) ),
  control = list(trace = 1)
)
opt$run_time = Sys.time() - start_time

# Check with TMB
start_time = Sys.time()
obj2 = MakeADFun(
  nll,
  parlist,
  random = "x",
  silent = TRUE
)
opt2 = nlminb( obj2$par, obj2$fn, control = list(trace = 1) )
opt2$run_time = Sys.time() - start_time

################
# Debugging
################

# FD and AD gradient should match ... confirmed :)
obj2$fn(obj2$par)
obj2$gr(obj2$par)
grad( obj2$fn, obj2$par )

# FD and SLM should match ... they do NOT
obj$fn(obj$par)
obj$gr(obj$par, fixed_Q = FALSE)
grad( obj$fn, obj$par )

################
# Debug gradient components
################

# FD and SLM should match given envelope theorem ... they DO
obj$gr(obj$par, what = "all")$grad_jnll
jacobian(
  func = \(x) obj$fn(x, what = "all")$jnll,
  x = obj$par
)

#
obj$gr(obj$par, what = "all", fixed_Q = FALSE)$grad_logdet
jacobian(
  func = \(x) mean(obj$fn(x, what = "all")$logdet),
  x = obj$par
)


################
# Check
################

H = obj2$env$spHess(par = obj2$env$par, random = TRUE)
ones = rep(1, sum(obj2$env$lrandom()))
t(ones) %*% H %*% ones

################
# Check implicit solution for duhat_dtheta
################

func = nll
parameters = parlist
random = "x"
k = 40
profile = NULL
m = 3
method = "newton_CG"
seed = 123
make_gr = TRUE
pu_update = "FD"
silent = TRUE
