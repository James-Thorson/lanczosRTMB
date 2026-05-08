
# Adapted from pracma::gmres
GMRES <-
function( b,
          Hq,
          x0 = rep(0, length(b)),
          errtol = 1e-06,
          kmax = length(b) + 1,
          reorth = 1 ){

  #stopifnot(is.numeric(A), is.numeric(b), is.matrix(A))

  zeros = pracma::zeros
  .givapp = pracma:::.givapp
  eye = pracma:::eye

  b <- as.matrix(b)
  n <- length(b)
  h <- zeros(kmax)
  v <- zeros(n, kmax)
  c <- zeros(kmax + 1, 1)
  s <- zeros(kmax + 1, 1)
  normF <- function(x) norm(as.matrix(x), type = "F")
  x <- as.matrix(x0)
  if (norm(x, "F") != 0) {
    r <- b - Hq(x)
  }
  else {
    r <- b
  }
  rho <- norm(r, "F")
  g <- rho * eye(kmax + 1, 1)
  errtol <- errtol * norm(b, "F")
  error <- c()
  error <- c(error, rho)
  niter <- 0
  if (rho < errtol) return(list(x = x, error = error, niter = niter))
  v[, 1] <- r/rho
  beta <- rho
  k <- 0
  while (rho > errtol && k < kmax) {
    k <- k + 1
    v[, k + 1] <- Hq(v[, k])
    normav <- normF(v[, k + 1])
    for (j in 1:k) {
      h[j, k] <- t(v[, j]) %*% v[, k + 1]
      v[, k + 1] <- v[, k + 1] - h[j, k] * v[, j]
    }
    h[k + 1, k] <- normF(v[, k + 1])
    normav2 <- h[k + 1, k]
    if ((reorth == 1 && normav + 0.001 * normav2 == normav) || reorth == 3) {
      for (j in 1:k) {
        hr <- t(v[, j]) %*% v[, k + 1]
        h[j, k] <- h[j, k] + hr
        v[, k + 1] = v[, k + 1] - hr * v[, j]
      }
      h[k + 1, k] <- normF(v[, k + 1])
    }
    if (h[k + 1, k] != 0){
      v[, k + 1] <- v[, k + 1]/h[k + 1, k]
    }
    if (k > 1){
      h[1:k, k] <- .givapp(c[1:(k - 1)], s[1:(k - 1)], h[1:k, k], k - 1)
    }
    nu <- normF(h[k:(k + 1), k])
    if (nu != 0) {
      c[k] <- Conj(h[k, k]/nu)
      s[k] <- -h[k + 1, k]/nu
      h[k, k] <- c[k] * h[k, k] - s[k] * h[k + 1, k]
      h[k + 1, k] <- 0
      g[k:(k + 1)] <- .givapp(c[k], s[k], g[k:(k + 1)], 1)
    }
    rho <- abs(g[k + 1])
    error <- c(error, rho)
  }
  y <- qr.solve(h[1:k, 1:k], g[1:k])
  x <- x0 + v[1:n, 1:k] %*% y
  return(list(x = x, error = error, niter = k))
}

