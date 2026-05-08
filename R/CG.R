
#' @title
#' Truncated conjugate gradient using Hessian-vector products
#'
#' @description
#' Approximate the soution to a linear system \eqn{Hq = b} using truncated conjugate gradient,
#'    without forming the Hessian matrix directly and instead using Hessian-vector products
#'    e.g., using forwards-on-reverse autodiff
#' Adapted from \code{mcmcsae::CG} under GPL-3 licence to use Hessian-vector products
#'
#' @param b vector to solve for
#' @param Hq efficient Hessian-vector product function, e.g., [make_Hq]
#' @param x initial guess for solution
#' @param Minv preconditioner matrix (uses Identity by default)
#' @param max.it maximum iterations (often should be less than \code{length(b)})
#' @param e error criterion
#' @param silent Be silent or print progress?
#'
#' @examples
#' library(Matrix)
#'
#' # Settings
#' n = 10^4
#' rho = 0.99
#'
#' # Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
#' P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
#' Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
#' x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
#' y = x + 0.1 * rnorm(n)
#' which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
#' y[-which_seen] = NA
#'
#' nll = function(p){
#'   -dgmrf(p$x, Q = Q, log = TRUE) - sum(dnorm(y, p$x, sd = 0.1, log=TRUE), na.rm=TRUE)
#' }
#' parlist = list( x=rnorm(n) )
#'
#' tape = MakeTape(nll, parlist)
#' gr = tape$jacfun()
#' Hq = make_Hq( tape, x = unlist(parlist) )
#' H = gr$jacfun(sparse = TRUE)
#'
#' #
#' b = gr(unlist(parlist))[1,]
#' Hess = H(unlist(parlist))
#' x = solve( Hess, b)
#' out = CG(
#'   b = b,
#'   Hq = Hq
#' )
#'
#' #
#' plot( x, out$x )
#'
#' @export
CG <-
function( b,
          Hq,
          x = 0 * b,
          Minv = Diagonal(n = length(b)),
          max.it = length(b),
          e = 1e-10,
          silent = TRUE ){
  # NOTES: could use Sturm with bisection to detect minimum eigenvalue for Lanczos, rather than min(eigen(H)$values)

  r <- b - Hq(x)
  z = as.vector(Minv %*% r)
  p <- z
  rz_old <- sum(r * z)
  alpha = beta = rep(NA, max.it)
  for(k in seq_len(max.it)){
    Ap <- Hq(p)
    alpha[k] <- rz_old / sum(p * Ap)
    x <- x + alpha[k] * p
    r <- r - alpha[k] * Ap
    z = as.vector(Minv %*% r)
    rz_new <- sum(r * z)
    discr <- sum(z * z)
    if (discr < e){
      break
    }
    beta[k] <- rz_new / rz_old
    p <- z + beta[k] * p
    rz_old <- rz_new
    if( isFALSE(silent) ){
      cat("iter:", k,"discr:",discr,"\n")
    }
  }
  list(
    x = x,
    k = k,
    alpha = alpha[seq_len(k)],
    beta = beta[seq_len(k)]
  )
}

#' @title
#' Nonlinear minimizer using line search with approximated Newton solution
#'
#' @description
#' Nonlinear minimizer designed for cheap Hessian-vector products using
#' [make_Hq], involving iterating a linear search along the truncated conjugate gradient
#' for a Newton.
#'
#' @inheritParams CG
#' @param par initial parameter vector
#' @param fn function to evaluate negative log-likelihood
#' @param gr function to evaluate gradient of negative log-likelihood
#' @param gr_tol early stopping condition for gradient of Newton solver
#' @param e_ratio early stopping condition for error of CG, relative to initial gradient
#' @param maxit_newton maximum iterations for Newton solver
#' @param maxit_CG maximum iterations for CG solution for each Newton iteration
#' @param c1 stopping condition for line search given CG solution in each Newton iteration
#' @param beta updates in line search stepsize alpha when Armijo sufficient decrease condition fails
#' @param silent Be silent or print progress?
#'
#' @examples
#' library(Matrix)
#' library(RTMB)
#'
#' # Settings
#' n = 10^4
#' rho = 0.99
#'
#' # Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
#' P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
#' Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
#' x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
#' y = x + 0.1 * rnorm(n)
#' which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
#' y[-which_seen] = NA
#'
#' nll = function(p){
#'   -dgmrf(p$x, Q = Q, log = TRUE) - sum(dnorm(y, p$x, sd = 0.1, log=TRUE), na.rm=TRUE)
#' }
#' parlist = list( x=rnorm(n) )
#'
#' tape = MakeTape(nll, parlist)
#' gr = tape$jacfun()
#' Hq = make_Hq( tape, x = unlist(parlist) )
#' H = gr$jacfun(sparse = TRUE)
#'
#' start_time = Sys.time()
#' opt1 = optim(
#'   par = unlist(parlist),
#'   fn = tape,
#'   gr = gr,
#'   method = "L-BFGS-B",
#'   control = list(
#'     maxit = 1e4
#'   )
#' )
#' opt1$runtime = Sys.time() - start_time
#'
#' opt2 = newton_CG(
#'   par = unlist(parlist),
#'   fn = tape,
#'   gr = gr,
#'   Hq = Hq
#' )
#'
#' # Compare the estimates and speed
#' matplot( cbind(opt1$par, opt2$par), type = "l", col = c("black","blue","red"), lty = "solid" )
#' c(opt1$runtime, opt2$runtime)
#'
#' @export
newton_CG <-
function( par,
          fn,
          gr,
          Hq,
          gr_tol = 1e-8,
          e_ratio = 0.01,
          maxit_newton = 100,
          maxit_CG = min(100,length(par)),
          c1 = 0.01,
          beta = 0.5,
          silent = FALSE ){

  start_time = Sys.time()
  x = par
  grad = gr(x)[,1]
  nll = fn(x)
  alpha_iter = CG_iter = rep(NA, maxit_newton)
  for( newton_iter in seq_len(maxit_newton) ){
    # CG for H^-1 grad
    step = CG(
      b = grad,
      Hq = \(q) Hq(q,x),
      max.it = maxit_CG,
      e = e_ratio * sum(grad^2)
    )
    CG_iter[newton_iter] = step$k

    # Update trust-region
    Tri = tridiag( step$alpha, step$beta[seq_len(length(step$alpha)-1)] )
    eigval = eigen(Tri, only.values = TRUE, symmetric = TRUE)$values
    t_new = max(0, -min(eigval) + 1e-3)

    ## Line search using c1 and beta with Armijo sufficient decrease condition
    alpha = 1
    grad_step = sum(grad * step$x)
    for( i in 1:20 ){
      x_test = x - alpha * step$x
      nll_test = fn(x_test)
      if( !is.nan(nll_test) ){
        if( nll_test < (nll - c1 * alpha * grad_step) ){
          break
        }
      }
      alpha = alpha * beta
    }
    alpha_iter[newton_iter] = alpha

    # First checks
    #if(nll_test > nll){
    #  grad = as.vector(gr(x_test))
    #  max_abs_grad = max(abs(grad))
    #  break
    #  #browser()
    #  #stop("Problem with newton_CG line search")
    #}

    # Update stuff
    nll = nll_test
    x = x_test
    grad = as.vector(gr(x))
    if( isFALSE(silent) ){
      cat("value:", nll,"mgc:",max(abs(grad)),"\n")
    }
    max_abs_grad = max(abs(grad))
    if( max_abs_grad < gr_tol ){
      break
    }
  }
  out = list(
    value = nll,
    par = x,
    max_abs_grad = max_abs_grad,
    runtime = Sys.time() - start_time,
    newton_iter = newton_iter,
    alpha_iter = alpha_iter[seq_len(newton_iter)],
    CG_iter = CG_iter[seq_len(newton_iter)]
  )
  out
}
