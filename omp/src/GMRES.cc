/**
 * @file GMRES.cc
 * @author cheny16
 * @brief GMRES source file
 * @version 1.0
 * @date 2023-03-25
 * 
 */
#include "../include/GMRES.hpp"
#include "../include/Matrix.hpp"
#include "../include/Givens.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>

#ifdef USE_MKL
    #include <mkl.h>
    #include <mkl_cblas.h>
#else
    #include <cblas.h>
#endif

#ifdef _OPENMP
    #include <omp.h>
#endif

/**
 * @brief Arnoldi iterations (Modified Gram-Schmidt)
 * 
 * @param A [In] Matrix A
 * @param v [In/Out] Matrix v Orthogonal vectors (rows)*(iter+1)
 * @param H [In/Out] Upper Hessenberg matrix H
 * @param rows Rows of matrix
 * @param iter The current iter
 * @param maxiter Max iter
 */
void Arnoldi(const Matrix A, Matrix v, Matrix H, 
             const int rows, const int iter, const int maxiter) {
    // Calculate w = A * v[iter] - Krylov vector
    Vector w = AllocVector(rows);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rows, rows, 1.0, A, rows, &v[iter * rows], 1, 0.0, w, 1);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (auto i = 0; i <= iter; i++) {
        double dot = cblas_ddot(rows, w, 1, &v[i * rows], 1);
        H[i * maxiter + iter] = dot;
        cblas_daxpy(rows, -dot, &v[i * rows], 1, w, 1);
    }

    double norm_w = cblas_dnrm2(rows, w, 1);
    H[(iter + 1) * maxiter + iter] = norm_w;
    if (norm_w != 0) {
        cblas_dscal(rows, 1.0 / norm_w, w, 1);
    }

    cblas_dcopy(rows, w, 1, &v[(iter + 1) * rows], 1);
    
    FreeVector(w);
}

/**
 * @brief Solve Least-Square Problem
 * 
 * @param H [In] Hessenberg matrix H
 * @param v [In] Arnoldi vectors
 * @param e [In] Residuals
 * @param x [In/Out] Solution of x
 * @param rows Rows
 * @param iter The curent iter
 * @param maxiter Max iter
 */
void SolveLSP(const Matrix H, const Matrix v, const Vector e, Vector x, 
              const int rows, const int iter, const int maxiter) {
    Vector y = AllocVector(iter + 1);
    // solve y using back substitution
    // The last line
    y[iter] = e[iter] / H[iter * maxiter + iter];
    // The rest lines
    for (auto i = iter - 1; i >= 0; i--) {
        double sum = 0.0;

#ifdef _OPENMP
        // #pragma omp parallel for reduction (+:sum)
        #pragma omp simd reduction (+:sum)
#endif 
        for (auto j = i + 1; j <= iter; j++) {
            sum += H[i * maxiter + j] * y[j];
        }
        y[i] = (e[i] - sum) / H[i * maxiter + i];
    }

    // solve x = x0 + v*y
    // Size of v: rows * (iter+1)
    cblas_dgemv(CblasColMajor, CblasNoTrans, rows, iter + 1, 1.0, v, rows, y, 1, 1.0, x, 1);
}

/**
 * @brief Calculate residual vector of the current sulotion
 * 
 * @param A [In] Matrix A
 * @param b [In] Vector b
 * @param x [In] Vector x The current solution
 * @param r [In/Out] Vector Res Residual Vector
 * @param rows Rows of matrix
 */
void CalResidualVec(const Matrix A, const Vector b, const Vector x, Vector r, const int rows) {
    // Calculate r = b - A*x
    // 1st: r = -A*x
    cblas_dgemv(CblasRowMajor, CblasNoTrans, rows, rows, -1.0, A, rows, x, 1, 0.0, r, 1);
    // 2nd: r = r + b
    // cblas_daxpy(rows, 1.0, b, 1, r, 1);
#ifdef _OPENMP
    #pragma omp parallel for
#endif    
    for (auto i = 0; i < rows; i++) {
        // r[i] = fabs(b[i] - r[i]);
        r[i] = b[i] - r[i];
    }
}

/**
 * @brief GMRES Method for Ax = b
 * 
 * @param A [In] Matrix A (row-major)
 * @param b [In] Vector b: Right-hand side
 * @param x [In/Out] Vector x: Initial guess
 * @param rows Rows of matrix
 * @param tol Convergence tolerance
 * @param maxiter The number of maxiter
 * @param verbose if_verbose
 */
void GMRES(const Matrix A, const Vector b, Vector x, 
           const int rows, const double tol, const int maxiter, bool verbose) {
    Vector e = AllocVector(maxiter+1);            /* Residual values */
    double residuals[maxiter+1];

    // Calculate initial residual vector of the initial guess of x0
    Vector r = AllocVector(rows);                 /* Initial residual r0 */
    CalResidualVec(A, b, x, r, rows);

    // beta = norm2(r0)
    double beta   = cblas_dnrm2(rows, r, 1);
    double norm_b = cblas_dnrm2(rows, b, 1);
    residuals[0] = beta / norm_b;

    if (residuals[0] <= tol) {
        std::cerr << "Convergence before iterations\n";
    }
    e[0] = beta;

    if (verbose) {
        std::cerr << "\nResidual before iterations is: " << residuals[0] << "\n";
    }

    Matrix H = AllocMatrix(maxiter+1, maxiter);   /* Hessenberg matrix */
    Matrix v = AllocMatrix(rows, maxiter+1);      /* Arnoldi vectors */

    Vector c = AllocVector(maxiter);              /* Givens rotator */
    Vector s = AllocVector(maxiter);              /* Givens rotator */
    Vector y = AllocVector(maxiter+1);            /* Solutions of LSP */

    // Scale r
    cblas_dscal(rows, 1.0 / beta, r, 1);
    // Copy r to v
    cblas_dcopy(rows, r, 1, v, 1);

    // Main GMRES iteration
    auto iter = 0;
    for ( ; iter < maxiter; iter++) {
        // 1. Arnoldi part
        Arnoldi(A, v, H, rows, iter, maxiter);/* New Arnoldi vector for next iter */
        
        // 2. Givens rotations
        // Apply old rotationer to the j-th column
        for (auto j = 0; j < iter; j++) {
            GivensRotate(c[j], s[j], &H[j * maxiter + iter], &H[(j+1) * maxiter + iter]);
        }

        // Generate new Givens rotator
        GenGivensRotator(H[iter * maxiter + iter], H[(iter + 1) * maxiter + iter], &c[iter], &s[iter]);

        // Apply Givens rotator to res_vals (update residual)
        H[iter * maxiter + iter] = c[iter] * H[iter * maxiter + iter] + \
                                   s[iter] * H[(iter + 1) * maxiter + iter];
        H[(iter + 1) * maxiter + iter] = 0.0;

        // Update the residual
        e[iter + 1] = -s[iter] * e[iter];
        e[iter    ] =  c[iter] * e[iter];

        double res_value = fabs(e[iter + 1]) / norm_b;
        residuals[iter] = res_value;

        if (verbose) {
            std::cerr << "Residual in " << iter+1 << " iter is: " << res_value << "\n";
        }
        std::cout << iter + 1 << " " << std::setprecision(10) << res_value << "\n";

        // 3. Check for convergence
        if (res_value <= tol) {
            if (verbose) std::cerr << "\nConvergence in " << iter+1 << " steps\n";

            // 4. Solve the upper triangular linear system
            SolveLSP(H, v, e, x, rows, iter, maxiter);

            break;
        }
        // Next iter
    } // Main GMRES iteration

    if (iter >= maxiter) {
        std::cerr << "\nNo convergence after max iter\n";
    }

    FreeMatrix(H); FreeMatrix(v);
    FreeVector(e); FreeVector(r); FreeVector(c); FreeVector(s); FreeVector(y);
}