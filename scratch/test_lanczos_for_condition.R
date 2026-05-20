
library(lanczosRTMB)
library(Matrix)

# Settings
n = 10^4
rho = 0.99

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
eigval = eigen(Q)$values
mean(eigval^4)

#
Hq = function(q,x) as.vector(Q%*%q)
L = lanczos(
  Hq,
  q = sample( c(-1,1), size = n, replace = TRUE ),
  k = 30,
  x = rep(0, n)
)
T = tridiag( L$alpha, L$beta )
eigval2 = eigen(T)$values
mean(eigval2^4)

