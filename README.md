# Parallel GMRES
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

## Parallel Parts
- [ ] Parallel Arnoldi
- [x] Parallel QR

## Modules / Routines
- Matrix
  | Routines      | Brief           |
  | --------      | -----           |
  |allocate_matrix| Allocate matrix |
  |dealloc_matrix | Free matrix     |
  |read_matrix    | Read matrix from SuiteSparse Matrix file (.mtx) |
  |print_matrix   | Print matrix to stdout |
  |MMmultiply     | Matrix matrix multiplication |
  |MVmultiply     | Matrix vector multiplication |
- ParaQR
  | Routines      | Brief           |
  | --------      | -----           |
  |GenHHMatrix    | Generate Householder reflection matrix |
  |QRDecomp       | QR Decomposition with Householder reflections |
- GMRES
  | Routines      | Brief           |
  | --------      | -----           |
  |Arnoldi        | Arnoldi iterations (MGS) |
  |BackSub        | Back substitution for solving $Rx=Q^Tb$ |
  |CalResidual    | Calculate residual vector of the current sulotion |
  |GMRES          | GMRES iteration with Householder QR |

## How to Build 
### Requirements
- GCC >= 9.3.0
- CMake >= 3.15
- MPI (like OpenMPI)
- OpenBLAS >= 0.3.0
- Boost >= 1.32.0 (with Program-Options module)

### Buiding
```bash
$ git clone https://github.com/cheny16/Parallel_GMRES.git
$ cd Parallel_GMRES
$ git checkout main
$ mkdir build && cd build
$ cmake ..
$ make
```

## Usage
```bash
# in build/
$ mpirun -np <num_of_procs> ./ParallelGMRES --help    <display helping messages> \
                                            --input   <input matrix file (.mtx)> \
                                            --tol     <convergence tolerance>    \
                                            --maxiter <max iterations>           \
                                            --timing  <display timing>
```
For example:
```bash
mpirun -np 2 ./ParallelGMRES --input ./datas/saylr1.mtx --tol 0.0001 --maxiter 20 --timing
```

## TODO
- [ ] Preconditioner
- [ ] Restarting
