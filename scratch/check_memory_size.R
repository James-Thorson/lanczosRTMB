
library(lanczosRTMB)
library(Matrix)
source( R'(C:\Users\James.Thorson\Desktop\Git\lanczosRTMB\R\lanczos.R)' )

# Settings
n = 10^5
rho = 0.99

# Simulate AR1 process approaching random walk (i.e., ill-conditioned inner problem)
P = bandSparse( n = n, k = c(-1), diagonals = list(rep(1,n)) )
Q = (Diagonal(n) - rho * t(P)) %*% (Diagonal(n) - rho * P)
x = RTMB:::rgmrf0( n= 1, Q = Q )[,1]
y = rpois( n = n, lambda = exp(1 + x) )
which_seen = sample( seq_len(n), size = n/10, replace = FALSE)
y[-which_seen] = NA

nll = function(p){
  -dgmrf(p$x, Q = Q, log = TRUE) - sum(dpois(y, lambda = exp(1 + p$x), log=TRUE), na.rm=TRUE)
}
parlist = list( x=rnorm(n) )

tape = MakeTape(nll, parlist)
gr = tape$jacfun()

# Before/after comparison at the process level
if(exists("Hq")) rm("Hq"); gc(); Sys.sleep(0.5)
mem_before <- ps::ps_memory_info()["rss"]
Hq = make_Hq( tape, x = unlist(parlist) )
gc(); Sys.sleep(0.5)  # sleep lets OS settle RSS
mem_after <- ps::ps_memory_info()["rss"]
cat("Approx C++ allocation:", (mem_after - mem_before) / 1e9, "GB\n")
#Hq(unlist(parlist))
