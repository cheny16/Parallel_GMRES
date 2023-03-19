/**
 * @file GMRES.cc
 * @author cheny16
 * @brief GMRES source file
 * @version 0.1
 * @date 2023-03-19
 * 
 */
#include "../include/GMRES.hpp"
#include "../include/Matrix.hpp"
#include "../include/ParaQR.hpp"

#include <cblas.h>
#include <cmath>
/**
 * @brief Arnoldi iterations (MGS)
 * 
 * @param A Matrix A
 * @param V Vector V
 * @param H [In/Out] Upper Hessenberg matrix H
 * @param rows Rows of matrix
 * @param iter The current iter
 * @return Vector 
 */
Vector Arnoldi(const Matrix A, const Vector V, Matrix H, const int rows, const int iter) {
    Vector W = MVmultiply(A, V, rows, rows, rows);
    for (auto i = 0; i < iter; i++) {
        float dot = cblas_sdsdot(rows, 1, V, 1, W, 1);
        H[i][iter] = dot;
        for (auto j = 0; j < rows; j++) {
            W[j] = W[j] - dot*V[j]; 
        }
    }
    return W;
}

/**
 * @brief Back substitution for solving Rx=QTb
 * 
 * @param R [In] Matrix R
 * @param Q [In] Matrix Q
 * @param B [In] Vector B
 * @param X [In/Out] Vector X Solution
 * @param rows Rows of matrix
 */
void BackSub(const Matrix R, const Matrix Q, const Vector B, Vector X, const int rows) {
    // The last line
    X[rows-1] = B[rows-1] / R[rows-1][rows-1];
    // The rest lines
    for (auto i = 1; i < rows; i++) {
        float sum = 0.0;
        for (auto j = 0; j < i; j++) {
            sum = sum + R[rows-1-i][rows-1-j] * X[rows-1-j];
        }
        sum = B[rows-1-i] - sum;
        X[rows-1-i] = sum / R[rows-1-i][rows-1-i];
    }
}

/**
 * @brief Calculate residual vector of the current sulotion
 * 
 * @param A [In] Matrix A
 * @param B [In] Vector B
 * @param X [In] Vector X The current solution
 * @param R [In/Out] Vector Res Residual Vector
 * @param rows Rows of matrix
 */
void CalResidual(const Matrix A, const Vector B, const Vector X, Vector Res, const int rows) {
    // Res = B - A*X
    Res = MVmultiply(A, X, rows, rows, rows);
    for (auto i = 0; i < rows; i++) {
        Res[i] = std::abs(B[i] - Res[i]);
    }
}

/**
 * @brief GMRES iteration with Householder QR
 * 
 * @param A [In] Matrix A
 * @param B [In] Vector B
 * @param X [In/Out] Vector X
 * @param rows Rows of matrix
 * @param tol Convergence tolerance
 * @param maxiter The number of maxiter
 * @param comm MPI_Comm
 */
void GMRES(const Matrix A, const Vector B, Vector X, 
           const int rows, const float tol, const int maxiter, MPI_Comm comm) {
    float residual[rows]{0};
    // Initial residual of the initial guess of x0
    Vector V = MVmultiply(A, X, rows, rows, rows); /* Initial residual */
    CalResidual(A, B, X, V, rows);

    float normR0 = cblas_snrm2(rows, V, 1);
    residual[0] = normR0;

    // v1 = r0/normR0
    for (auto i = 0; i < rows; i++) {
        V[i] = V[i] / normR0;
    }
    Matrix H = allocate_matrix(rows+1, rows);
    // Main GMRES iteration
    for (auto iter = 0; iter < maxiter; iter++) {
        // 1. Arnoldi part
        Vector W = Arnoldi(A, V, H, rows, iter);
        // h(j+1,j)=norm2(w)
        float normW = cblas_snrm2(rows, W, 1);
        H[iter+1][iter] = normW;
        if (normW != 0) {
            for (auto i = 0; i < rows; i++) {
                V[i] = W[i] / normW;
            }
        }
        else 
            break;
        
        // 2. Householder QR
        // Allocate and initialize
        Matrix Q = allocate_matrix(rows ,rows);
        Matrix R = allocate_matrix(rows, rows);
        for (auto i = 0; i < rows; i++) {
            for (auto j = 0; j < rows; j++) {
                R[i][j] = A[i][j];
            }
        }
        QRDecomp(Q, R, rows, comm);

        // 3. Back substitution for Rx=QTb
        BackSub(R, Q, B, X, rows);

        // 4. Update residual
        CalResidual(A, B, X, V, rows);
        float normV = cblas_snrm2(rows, V, 1);
        residual[iter] = normV;

        // 5. Check convergence
        if (normV < tol) {
            std::cerr << "Convergence in " << iter << " steps\n";
            break;
        }
    }
    std::cerr << "No convergence after max iter\n";
}