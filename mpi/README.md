# Parallel GMRES - MPI version

To be finished...

## Parallelizable Parts

| GMRES Kernel                                | Parallel/Sequential | MPI Routines |
| ------------------------------------------- | ------------------- | ------------ |
| Arnoldi Part - sparse matrix vector product | parallel            | Allgatherv   |
| Arnoldi Part - dot product                  | parallel            | Allreduce    |
| Givens rotations                            | sequential          | -            |
| Solving the LSP                             | sequential          | -            |

