/**
 * @file QR.hpp
 * @author cheny16
 * @brief QR header file
 * @version 0.1
 * @date 2023-03-17
 * 
 */
#ifndef QR_HPP
#define QR_HPP

#include "Matrix.hpp"
#include "mpi.h"

const int MAX_PROCS = 50;

Matrix GenHHMatrix(float *V, int rows);

void QRDecomp(Matrix Q, Matrix R, int rows, MPI_Comm comm);

#endif /* QR_HPP */