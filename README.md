# lanczosRTMB
Lanczos methods for RTMB

### Uses
_lanczosRTMB_ is designed to approximate properties of a hierarchical model without ever constructing or inverting the matrix of second derivatives ("Hessian matrix").  It can be used to approximate:

* Standard errors using a matrix-free delta method
* The log-determinant of the Hessian using stochastic trace estimation
* The log marginal likelihood using the Laplace approximation 
