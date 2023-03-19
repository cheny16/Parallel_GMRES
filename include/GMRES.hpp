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
#include "mpi.h"

Vector Arnoldi(const Matrix A, const Vector V, Matrix H, const int rows, const int iter);

void BackSub(const Matrix R, const Matrix Q, const Vector B, Vector X, const int rows);

void CalResidual(const Matrix A, const Vector B, const Vector X, Vector R, const int rows);

void GMRES(const Matrix A, const Vector B, Vector X, 
           const int rows, const float tol, const int maxiter, MPI_Comm comm);

#endif /* GMRES_HPP */