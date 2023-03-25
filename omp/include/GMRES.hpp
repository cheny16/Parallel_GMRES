/**
 * @file GMRES.hpp
 * @author cheny16
 * @brief GMRES header file
 * @version 0.1
 * @date 2023-03-19
 * 
 */
#ifndef GMRES_HPP
#define GMRES_HPP

#include "Matrix.hpp"

void Arnoldi(const Matrix A, Matrix v, Matrix H, 
             const int rows, const int iter, const int maxiter);

void SolveLSP(const Matrix H, const Matrix v, const Vector e, Vector x, 
              const int rows, const int iter, const int maxiter);

void CalResidualVec(const Matrix A, const Vector b, const Vector x, Vector r, const int rows);

void GMRES(const Matrix A, const Vector b, Vector x, 
           const int rows, const double tol, const int maxiter, bool verbose=false);

#endif /* GMRES_HPP */