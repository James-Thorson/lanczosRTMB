## lanczosRTMB

[![docs
(main)](https://img.shields.io/badge/docs-lanczosRTMB-orange.svg?colorB=E91E63)](https://james-thorson.github.io/lanczosRTMB/)

Lanczos methods for RTMB

### Uses

*lanczosRTMB* is designed to approximate properties of a hierarchical
model implemented in *RTMB* without ever constructing or inverting the
matrix of second derivatives (“Hessian matrix”). It can be used to
approximate:

- Standard errors using a matrix-free delta method
- The log-determinant of the Hessian using stochastic trace estimation
- The log marginal likelihood using the Laplace approximation

### Potential developments

There’s many avenues for further improvement. These include:

Adding a quasi-Newton inner optimizer using `Hq`, i.e., a variant of
conjugate gradient;

Moving the FD logic for the `obj$gr` to tape;

Exploring a conjugate gradient algorithm for an implicit solution to the
gradient of the log-determinant;
