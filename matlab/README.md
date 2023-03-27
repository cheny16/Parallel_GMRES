# GMRES - MATLAB version
It is just a sequatial MATLAB implementation of GMRES algorithm for solving $Ax=b$, as a comparison with the C++ version.

## Usage

Just put the `GMRES.m` file into the current working directory, pass matrix `A` and right-hand `b` to the routine:

```matlab
[X, residuals] = GMRES(A, b, tol, maxiter)
```

where `tol` is the convergence tolerance, `maxiter` is the max number if iterations.

After calling, `X` would be the solution, `residuals` would be the residual values.
