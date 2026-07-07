
library(Matrix)
library(lanczosRTMB)

###############
# Lanczos for covariate diffusion
#  SPDE:  Dinv = I + kappa^-2 C^-1 G
#  AR:  Dinv = I - P
###############

n = 1e2 + 1
kappa = exp(- 2 / n )
P = bandSparse( n = n, k = c(-1,1), diag = list(rep(0.5,n),rep(0.5,n)) )
I = V = Diagonal(n = n)
invD = I - kappa * P
#Q = (I - kappa * t(P)) %*% solve(V) %*% (I - kappa * P)
x = rep(0, n)
x[ceiling(n/2)] = 1

# Exact
if( n <= 1001 ){
  D0 = expm( kappa*P - Diagonal(n=n) )
}else{
  D0 = sparseMatrix( i = 1, j = 1, x = 0, dims = c(n,n) )
}

# Implicit
if( n <= 10001 ){
  D = solve(invD)
}else{
  D = sparseMatrix( i = 1, j = 1, x = 0, dims = c(n,n) )
}

Hq = \(q,x,update_H) (invD %*% q)[,1]

# Lanczos
if( n <= 100001 ){
  L = lanczos(Hq, q = x, k = 1000, x = rep(0,n) )
  T = tridiag( L$alpha, L$beta )
  D1 = L$Q %*% solve(T) %*% t(L$Q)
}else{
  D1 = sparseMatrix( i = 1, j = 1, x = 0, dims = c(n,n) )
}

# CG
cg_solve = CG( x, Hq, x = x, e = 1e-4 )
cg_solve = CG( x, Hq, x = x, e = 0, max.it = 2 )

Y = cbind(D0 %*% x, D %*% x, D1 %*% x, cg_solve$x )
matplot(
  scale(Y, scale = colSums(Y), center = FALSE ),
  type = "l",
  lwd = 2,
  lty = "solid",
  col = rainbow(4)
)
legend( "topright", bty = "n", fill = rainbow(4),
        legend = c("expm","implicit","lanczos","CG") )



###############
# Initial exploration
###############

#n = 100
#kappa = exp(-10 / n )
#P = bandSparse( n = n, k = c(-1,1), diag = list(rep(0.5,n),rep(0.5,n)) )
#I = V = Diagonal(n = n)
#Q = (I - kappa * t(P)) %*% solve(V) %*% (I - kappa * P)
#
#invD = I - kappa * P
#D = solve(invD)
#
#Hq = \(q,x) (invD %*% q)[,1]
#q = sample( c(-1,1), size = n, replace = TRUE )
#L = lanczos(Hq, q = q, k = 30, x = rep(0,n) )
#T = tridiag( L$alpha, L$beta )
#
#D1 = L$Q %*% solve(T) %*% t(L$Q)
#
#cols = ceiling(seq(1,n,length = 6))
#matplot(D1[,cols], type = "l", col = viridisLite::viridis(n), lty = "solid" )

