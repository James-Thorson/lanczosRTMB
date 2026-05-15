
library(Matrix)
library(RTMB)
library(lanczosRTMB)

## Settings
#n = 10^5
#rho = 0.999
#
##
#P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
#Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
#x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
#y = x + 0.1 * rnorm(n)
#which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
#y[-which_seen] = NA
#
#nll = function(p){
#  -dgmrf(p$x, Q = Q, log = TRUE) - sum(dnorm(y, p$x, sd = 0.1, log=TRUE), na.rm=TRUE)
#}
#parlist = list( x=rnorm(n) )

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

# Break sparsity with a change-of-support contribution
nll = function(p){
  -dgmrf(p$x, Q = Q, log = TRUE) - sum(dpois(y, exp(1 + p$x), log=TRUE), na.rm=TRUE) - dnorm(log(sumy), log(sum(exp(1 + p$x))), sd = 0.01, log = TRUE)
}
parlist = list( x=rnorm(nx*ny) )

tape = MakeTape(nll, parlist)
gr = tape$jacfun()
Hq = make_Hq( tape, x = unlist(parlist) )
#H = gr$jacfun(sparse = TRUE)
H = gr$jacfun(sparse = FALSE)

#source( R'(C:\Users\James.Thorson\Desktop\Git\lanczosRTMB\scratch\CG.R)')
start_time = Sys.time()
opt1 = optim(
  par = unlist(parlist),
  fn = tape,
  gr = gr,
  method = "L-BFGS-B",
  control = list(
    maxit = 1e5,
    factr = 1e4  # must be lower than default
  )
)
opt1$runtime = Sys.time() - start_time
opt1$grad = sqrt(sum(gr(opt1$par)^2))

opt2 = newton_CG(
  par = unlist(parlist),
  fn = tape,
  gr = gr,
  Hq = Hq,
  gr_tol = 0.001,
  #tol10 = 0, # 0 = turn off early stop
  smartsearch = FALSE,
  stop_if_nonPD = TRUE
)

opt3 = newton_CG(
  par = unlist(parlist),
  fn = tape,
  gr = gr,
  Hq = Hq,
  gr_tol = 0.001,
  smartsearch = TRUE,
  stop_if_nonPD = TRUE
)

matplot( cbind(opt1$par, opt2$par), type = "l", col = c("black","blue","red"), lty = "solid" )
cbind(opt1$runtime, opt2$runtime)
cbind(opt1$grad, opt2$grad)
cor( opt1$par, opt2$par )
cor( opt2$par, opt3$par )

