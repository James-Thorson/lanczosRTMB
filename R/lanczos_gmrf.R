

#' @title
#' Approximate the log-density of a GMRF using Lanczos
#'
#' @description
#' Calculate log-density of a GMRF using Lanczos recursion
#'
#' @inheritParams lanczos
#' @inheritParams lanczos_logdet
#' @inheritParams RTMB::dgmrf
#'
#' @examples
#' library(RTMB)
#' library(Matrix)
#' library(lanczosRTMB)
#'
#' # Create a precision for a 1D AR1 process
#' P = bandSparse( n = 100, k = c(-1,1), diagonals = list(rep(0.5,99),rep(0.5,99)) )
#' Q = (Diagonal(100) - 0.5*t(P) ) %*% (Diagonal(100) - 0.5*P)
#' x = RTMB:::rgmrf0(n = 1, Q )[,1]
#'
#' # Exact density using RTMB
#' dgmrf(x, Q = Q, log = TRUE)
#'
#' # Lanczos-approximated density
#' lanczos_dgmrf(x, Q = Q, k = 20, log = TRUE)
#'
#' @export
lanczos_dgmrf <-
function( x,
          Q,
          k,
          mu = rep(0,length(x)),
          m = 3,
          log = FALSE,
          seed = NULL,
          orthogonalize = TRUE ){

  delta = x - mu
  Hq = function(x,...) (Q %*% x)[,1]
  env = list(
    which_random = seq_along(x),
    x = rep(0, length(x))
  )
  attr(Hq,"env") = env
  logdet_m = lanczos_logdet(
    Hq = Hq,
    k = k,
    m = m,
    orthogonalize = orthogonalize,
    seed = seed
  )

  nll = -0.5*sum(delta * (t(delta) %*% Q)) + 0.5*logdet_m - length(x)/2*log(2*pi)
  if(isTRUE(log)){return(nll)}else{return(exp(nll))}
}





