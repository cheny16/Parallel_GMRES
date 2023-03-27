# Parallel GMRES - OpenMP version

This is an implementation of OpenMP shared memory parallel GMRES for solving sparse linear system.

## How to Build 
### Requirements
- GCC >= 9.3.0
- CMake >= 3.15
- OpenBLAS >= 0.3.0 (using OpenBLAS)
- Intel MKL (using MKL)
- Boost >= 1.32.0 (with Program-Options module)

### Buiding
```bash
$ git clone https://github.com/cheny16/Parallel_GMRES.git
$ cd Parallel_GMRES
$ git checkout main
$ mkdir build && cd build
# if use MKL (default: OFF)
$ cmake .. -D USE_MKL=ON
$ make
```

## Usage
```bash
# in build/
$ ./gmres --help    <display helping messages> \
          --input   <input matrix file (.mtx)> \
          --tol     <convergence tolerance>    \
          --maxiter <max iterations>           \
          --verbose <display details>          \
          --timing  <display timing>
```
For example:

Apply `gmres` on `saylr1.mtx` matrix, with `tolerance` is `0.00001`, and `maxiter` is `200`, to show details and timing, and output the results to `saylr1.out`:

`./gmres --i ../datas/saylr1.mtx --tol 0.00001 --m 200 --verbose --timing > saylr1.out`

## Modules / Routines
- Matrix
  | Routines      | Brief           |
  | --------      | -----           |
  |AllocMatrix    | Allocate matrix |
  |FreeMatrix     | Free matrix     |
  |ReadMatrix     | Read matrix from SuiteSparse Matrix file (.mtx) |
  |PrintMatrix    | Print matrix to stdout |
  |AllocVector    | Allocate vector |
  |FreeVector     | Free vector |
- Givens
  | Routines        | Brief           |
  | --------        | -----           |
  |GenGivensRotator | Generate Givens rotator |
  |GivensRotate     | Apply Givens rotation   |
- GMRES
  | Routines      | Brief           |
  | --------      | -----           |
  |Arnoldi        | Arnoldi iterations (MGS) |
  |SolveLSP       | Solve least-squares problem (LSP) |
  |CalResidualVec | Calculate residual vector of the current sulotion |
  |GMRES          | GMRES iteration with Givens rotations |