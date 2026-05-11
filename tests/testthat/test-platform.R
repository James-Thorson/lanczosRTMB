
test_that("newton_CG is working ", {
  library(Matrix)
  library(RTMB)

  # Settings
  n = 10^4
  rho = 0.99

  # Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
  P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
  Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
  x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
  y = x + 0.1 * rnorm(n)
  which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
  y[-which_seen] = NA

  nll = function(p){
    -dgmrf(p$x, Q = Q, log = TRUE) - sum(dnorm(y, p$x, sd = 0.1, log=TRUE), na.rm=TRUE)
  }
  parlist = list( x=rnorm(n) )

  tape = MakeTape(nll, parlist)
  gr = tape$jacfun()
  Hq = make_Hq( tape, x = unlist(parlist) )
  H = gr$jacfun(sparse = TRUE)

  opt1 = optim(
    par = unlist(parlist),
    fn = tape,
    gr = gr,
    method = "L-BFGS-B",
    control = list(
      maxit = 1e4
    )
  )

  opt2 = newton_CG(
    par = unlist(parlist),
    fn = tape,
    gr = gr,
    Hq = Hq
  )

  # Compare the estimates and speed
  expect_equal( opt1$par, opt2$par, tol = 1e-2 )
})

