
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
#' @param stop_if_nonPD whether to stop CG if recusion is not positive definite
#' @param silent Be silent or print progress?
#'
#' @importFrom utils tail
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
#' # Test unregularized solution
#' b = gr(unlist(parlist))[1,]
#' Hess = H(unlist(parlist))
#' x = solve( Hess, b)
#' out = CG(
#'   b = b,
#'   Hq = Hq
#' )
#' plot( x, out$x )
#'
#' # Test trust-region regularization
#' Hess2 = Hess + 2 * Diagonal(n=n)
#' x2 = solve( Hess2, b)
#' out = CG(
#'   b = b,
#'   Hq = \(b) Hq(b) + 2*b
#' )
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
          stop_if_nonPD = TRUE,
          silent = TRUE ){
  # NOTES: could use Sturm with bisection to detect minimum eigenvalue for Lanczos, rather than min(eigen(H)$values)

  r <- b - Hq(x)
  z = as.vector(Minv %*% r)
  p <- z
  rz_old <- sum(r * z)
  alpha = beta = rep(NA, max.it)
  status = 1    # 0 = non-PD;  1 = exceeds max.it;  2 = converged
  for(k in seq_len(max.it)){
    Ap <- Hq(p)
    denom <- sum(p * Ap)
    alpha[k] <- rz_old / denom
    # Check for non-PD after alpha, so can still use low-k solution if FALSE
    if( isTRUE(stop_if_nonPD) && (denom < -1e-8) ){
      status = 0
      break
    }
    x <- x + alpha[k] * p
    r <- r - alpha[k] * Ap
    z = as.vector(Minv %*% r)
    rz_new <- sum(r * z)
    discr <- sum(z * z)
    if (discr < e){
      status = 2
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
    beta = beta[seq_len(k)],
    status = status
  )
}

#' @title
#' Nonlinear minimizer using line search with approximated Newton solution
#'
#' @description
#' Nonlinear minimizer designed to use cheap Hessian-vector products using
#' [make_Hq].  The minimizer approximates a Newton update using truncated conjugate
#' gradient, while adapting a regularization parameter designed and using a
#' line-search for each step.
#'
#' @inheritParams CG
#' @inheritParams TMB::newton
#' @param par initial parameter vector
#' @param fn function to evaluate negative log-likelihood
#' @param gr function to evaluate gradient of negative log-likelihood
#' @param gr_tol early stopping condition for gradient of Newton solver
#' @param e_ratio early stopping condition for error of CG, relative to initial gradient
#' @param maxit_newton maximum iterations for Newton solver
#' @param maxit_CG maximum iterations for CG solution for each Newton iteration
#' @param line_steps number of steps to explore for linear-search given Newton update
#' @param c1 stopping condition for line search given CG solution in each Newton iteration, where
#'        \eqn{1>c1>0} seeks to prevent overshoot and resulting zigzagging.
#' @param beta updates in line search stepsize alpha when Armijo sufficient decrease condition fails
#' @param diagnostics whether to provide extra diagnostics for each Newton iteration
#' @param silent Be silent or print progress?
#'
#' @details
#' This minimizer approximates Newton steps \eqn{x_{i+1} = x_{i} - H(x_i)^{-1} g(x_i)},
#' where \eqn{H(x_i) = \nabla^2 f(x_i)} is the Hessian matrix and
#' \eqn{g(x_i) = \nabla f(x_i)} is the gradient of the negative log-likelihood
#' \eqn{f(x_i)}.
#' It then uses several strategies for numerical and computational efficiency.
#' \enumerate{
#' \item It approximates each Newton step \eqn{H(x_i)^{-1} g(x_i)} using a truncated
#'       conjugate gradient algorithm that involves Hessian-vector products \eqn{H(x_i) v}.
#'       It then recursively improves the solution until a desired accuracy is reached,
#'       controlled by \code{gr_tol}, often using many fewer steps than the full dimension
#'       \code{length(x_i)};
#' \item It approximates each Newton solution \eqn{H(x_i)^{-1} g(x_i)}
#'       without ever constructing or storing the Hessian \eqn{H(x_i)}, and instead computing
#'       Hessian-vector product \eqn{H(x_i) v = \nabla( v^T \nabla f(x_i) )} using
#'       autodifferentiation in RTMB;
#' \item After approximating each step direction \eqn{H(x_i)^{-1} g(x_i)},
#'       it uses a line-search algorithm while decreasing the step-size if needed
#'       to find a suitable decrease in the objective function, controlled
#'       \code{line_steps}, \code{beta}, and \code{c1};
#' \item The CG is monitored to identify whether
#'       the Hessian matrix is positive definite (PD), and if it is not PD then the CG may be
#'       terminated, controlled by \code{stop_if_nonPD}.
#' \item When using \code{smartsearch = TRUE}, if the Hessian is not PD, then
#'       a regularization \code{ustep} is decreased, corresponding to an increase in
#'       \eqn{t=\frac{1}{u}-1} where the Newton step is actually using \eqn{(H + tI)^{-1} g}
#'       where \eqn{g} is the gradient, and \eqn{1>u>0} corresponds to \eqn{t>0}.
#'       This increase or decrease in \code{ustep} (and associated decrease/increase
#'       in \eqn{t}) is copied from [TMB::newton()]. If \code{smartsearch = FALSE} then
#'       \eqn{ustep=1} and \eqn{t=0} such that CG uses the Hessian
#'       corresponding to unregularized Newton steps.
#'       This smartsearch behavior controlled by \code{u0}, \code{ustep},
#'       \code{power}, and \code{tol10}, and it corresponds to an adaptive "trust-region".
#' }
#'
#' @examples
#' library(Matrix)
#' library(RTMB)
#'
#' # Settings
#' n = 10^3
#' rho = 0.99
#'
#' # Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
#' P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
#' Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
#' x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
#' y = rpois( n = n, lambda = exp(1 + x) )
#' which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
#' y[-which_seen] = NA
#'
#' nll = function(p){
#'   -dgmrf(p$x, Q = Q, log = TRUE) - sum(dpois(y, lambda = exp(1 + p$x), log=TRUE), na.rm=TRUE)
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
#'   Hq = Hq,
#'   maxit_newton = 1e4,
#'   gr_tol = 1e-4
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
          line_steps = 20,
          smartsearch = TRUE,
          ustep = 1, ## Start out optimistic: Newton step
          power = 0.5, ## decrease=function(u)const*u^power
          u0 = 1e-4,  ## Increase u=0 to this value
          stop_if_nonPD = smartsearch,
          tol10 = 0.001,
          diagnostics = FALSE,
          silent = FALSE ){

  if(line_steps<1)stop("`line_steps` must be 1 or greater")
  start_time = Sys.time()
  norm <- function(x) sqrt(sum(x^2))
  x = par
  grad = as.vector(gr(x))
  nll = fn(x)
  t = 0  # Trust region regulization t >= 0
  fail = 0
  status_iter = nll_iter = grad_iter = ustep_iter = stepsize_iter = CG_iter = rep(NA, maxit_newton)

  ## Adaptive stepsize algorithm (smartsearch)
  phi <- function(u)1/u-1
  #invphi <- function(x)1/(x+1)           # IGNORED FOR NOW
  ## ========== Functions controling the algorithm
  ## Important requirements:
  ## 1. increase(u) and decrease(u) takes values in [0,1]
  ## 2. increase(u)>u and decrease(u)<u
  ## 3. increase(u)->1 when u->1
  ## 4. decrease(u)->0 when u->0
  ## Properties of algorithm:
  ## * ustep must converge towards 1 (because 1 <==> Positive definite hessian)

  ## power<1 - controls the boundary *repulsion*
  increase <- function(u)u0+(1-u0)*u^power
  ##decrease <- function(u)1-increase(1-u)
  ## Solve problem with accuracy when u apprach 0
  decrease <- function(u)ifelse(u>1e-10,1-increase(1-u),(1-u0)*power*u)
  ##plot(increase,0,1,ylim=c(0,1));plot(decrease,0,1,add=TRUE);abline(0,1)
  ustep = increase(ustep)

  for( newton_iter in seq_len(maxit_newton) ){
    # CG for H^-1 grad
    step = CG(
      b = grad,
      Hq = \(q) Hq(q,x) + phi(ustep)*q,
      max.it = maxit_CG,
      #e = e_ratio * sum(grad^2)
      e = max(e_ratio * sum(grad^2), 1e-8),
      stop_if_nonPD = stop_if_nonPD
    )
    CG_iter[newton_iter] = step$k
    status_iter[newton_iter] = step$status
    ustep_iter[newton_iter] = ustep

    # Update trust-region
    Tri = tridiag( step$alpha, step$beta[seq_len(length(step$alpha)-1)] )
    eigval = eigen(Tri, only.values = TRUE, symmetric = TRUE)$values
    min_eigval = min(0, min(eigval) - 1e-3)

    ## Line search using c1 and beta with Armijo sufficient decrease condition
    stepsize = 1
    step_status = 0
    for( i in seq_len(line_steps) ){
      x_test = x - stepsize * step$x
      nll_test = fn(x_test)
      if( !is.nan(nll_test) ){
        if( nll_test < (nll - c1 * stepsize * sum(grad * step$x)) ){
          step_status = 1
          break
        }
      }
      stepsize = stepsize * beta
    }
    stepsize_iter[newton_iter] = stepsize

    # First checks
    #if(newton_iter == 30 )browser()
    if( isTRUE(smartsearch) ){
      if( (step$status %in% c(1,2)) && (step_status==1) ){
        nll = nll_test
        x = x_test
        grad = as.vector(gr(x))
        # Contract
        #t = sqrt(t) * alpha^-0.5
        ustep = increase(ustep)
      }else{
        # Expand
        #t = t * alpha^-0.5
        ustep = decrease(ustep)
      }
      if( isFALSE(silent) ){
        cat("value:", nll,"mgc:",max(abs(grad)),"ustep:",ustep,"\n")
      }
    }else{
      nll = nll_test
      x = x_test
      grad = as.vector(gr(x))
      if( isFALSE(silent) ){
        cat("value:", nll,"mgc:",max(abs(grad)),"\n")
      }
    }

    # Check convervence
    nll_iter[newton_iter] = nll
    if( newton_iter > 10 ){
      tail10 <- tail(nll_iter[seq_len(newton_iter)],10)
      improve10 <- tail10[1] - tail10[length(tail10)]
      if(improve10 < tol10){
        if(isFALSE(silent))cat("Not improving much - will try early exit...")
        pd <- ifelse( status_iter[newton_iter] %in% c(1,2), TRUE, FALSE )
        if(isFALSE(silent))cat("PD hess?:",pd,"\n")
        if(pd) break
        fail <- fail+1
      }
    }
    #if(norm(par-parold)<step.tol){
    #  break
    #}
    #if(fail>5){
    #  stop("Newton drop out: Too many failed attempts.")
    #}

    # Update stuff
    grad_iter[newton_iter] = sqrt(sum(grad^2))
    if( grad_iter[newton_iter] < gr_tol ){
      break
    }
  }
  out = list(
    value = nll,
    par = x,
    grad = sqrt(sum(grad^2)),
    runtime = Sys.time() - start_time,
    newton_steps = newton_iter,
    CG_steps = sum(CG_iter, na.rm=TRUE)
  )
  if( isTRUE(diagnostics) ){
    out$diag = list(
      nll_iter = nll_iter[seq_len(newton_iter)],
      grad_iter = grad_iter[seq_len(newton_iter)],
      newton_iter = newton_iter,
      stepsize_iter = stepsize_iter[seq_len(newton_iter)],
      CG_iter = CG_iter[seq_len(newton_iter)],
      ustep_iter = ustep_iter[seq_len(newton_iter)],
      status_iter = status_iter[seq_len(newton_iter)]
    )
  }
  return(out)
}
