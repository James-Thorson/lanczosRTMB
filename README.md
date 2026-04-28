## lanczosRTMB
[![docs (main)](https://img.shields.io/badge/docs-lanczosRTMB-orange.svg?colorB=E91E63)](https://james-thorson.github.io/lanczosRTMB/)

Lanczos methods for RTMB

### Uses
_lanczosRTMB_ is designed to approximate properties of a hierarchical model implemented in _RTMB_ without ever constructing or inverting the matrix of second derivatives ("Hessian matrix").  It can be used to approximate:

* Standard errors using a matrix-free delta method
* The log-determinant of the Hessian using stochastic trace estimation
* The log marginal likelihood using the Laplace approximation 
