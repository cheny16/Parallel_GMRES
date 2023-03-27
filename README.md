# Parallel GMRES

<img src="https://img.shields.io/badge/build-passing-brightgreen" /> <img src="https://img.shields.io/badge/tests-passing-blue" />

Parallel GMRES for solving sparse linear systems $Ax=b$

## Algorithm (GMRES with Givens)
$Ax=b,\ where\ A\ is\ square\ matrix,\ given\ an\ initial\ guess\ x^{(0)}, tol\ \epsilon>0,\ and\ maxiter$

1. $r_0=b-Ax^{(0)}\ is\ the\ initial\ residual,\ \beta:=||r_0||_2,\ v_1:=r_0/\beta$

2. $For\ j=1,2,...,maxiter\ Do:$
3. $\ \ \ \ \ \ \ \ w=Av_j$ % Matrix-vector multiplication
4. $\ \ \ \ \ \ \ \ For\ i=1,2,...,j\ Do:$ % Arnoldi Part (MGS)
5. $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ h_{i,j}=dot(v_i,w)$
6. $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ w=w-h_{i,j}v_i$
7. $\ \ \ \ \ \ \ \ End\ for$
8. $\ \ \ \ \ \ \ \ h_{j+1,j}=||w||_2$
9. $\ \ \ \ \ \ \ \ If\ h_{j+1,j}\ \neq 0\ then:$
10. $\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ v_{j+1}=w/h_{j+1,j}$
11. $\ \ \ \ \ \ \ \ End\ if$
12. $\ \ \ \ \ \ \ \ [H, c, s]=GivensRotate(H, c, s, j)$ % Apply Givens rotations to upper Hessenberg matrix H
13. $\ \ \ \ \ \ \ \ Eliminate H(j+1, i)$
14. $\ \ \ \ \ \ \ \ Update\ the\ residual$
15. $\ \ \ \ \ \ \ \ If\ residual < \epsilon\ then:\ break$
16. $End\ for$
17. $Solve\ least-squares\ problem (back\ substitution)$

## Parts
- matlab/: MATLAB version 
- omp/: OpenMP version

## TODO
- [ ] MPI version
- [ ] Preconditioner
