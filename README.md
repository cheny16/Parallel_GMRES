# Parallel_GMRES
Parallel GMRES for solving sparse linear systems $Ax=b$

## Algorithm (GMRES with Householder QR)
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
12. $\ \ \ \ \ \ \ \ [Q, R]=HouseholderQR(H)$ % Apply Householder QR to upper Hessenberg matrix H
13. $\ \ \ \ \ \ \ \ Solve\ least-squares\ problem:\ Rx^*=Q^Tb\ (back\ substitution)$
14. $\ \ \ \ \ \ \ \ Update\ residual\ r_k(k)=abs(b-Ax^*)$
15. $\ \ \ \ \ \ \ \ If\ residual < \epsilon\ then:\ break$
16. $End\ for$
