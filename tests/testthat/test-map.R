
test_that("map is working ", {

  set.seed(123)
  x = rnorm(100)
  y = 1 + 0.5 * x + rnorm(100)

  get_nll = function( p ){
    yhat = p$par[1] + p$par[2] * x + p$par[3] * x
    -sum(dnorm(y, yhat, exp(p$par[4]), log = TRUE))
  }

  # Map off
  obj = lanczos_MakeADFun(
    func = get_nll,
    parameters = list( par = 1:4 ),
    map = list( par = factor(c(1,2,NA,4)) ),
    k = 10,
    random = NULL,
    inner_optimizer = "optim",
    make_gr = FALSE
  )
  obj0 = MakeADFun(
    func = get_nll,
    parameters = list( par = 1:4 ),
    map = list( par = factor(c(1,2,NA,4)) )
  )
  # Compare the estimates and speed
  expect_equal( obj$fn(obj$par), obj0$fn(obj0$par), tol = 1e-2 )

  # Map off
  obj = lanczos_MakeADFun(
    func = get_nll,
    parameters = list( par = 1:4 ),
    map = list( par = factor(c(1,2,2,4)) ),
    k = 10,
    random = NULL,
    inner_optimizer = "optim",
    make_gr = FALSE
  )
  obj0 = MakeADFun(
    func = get_nll,
    parameters = list( par = 1:4 ),
    map = list( par = factor(c(1,2,2,4)) )
  )
  # Compare the estimates and speed
  expect_equal( obj$fn(obj$par), obj0$fn(obj0$par), tol = 1e-2 )

})

