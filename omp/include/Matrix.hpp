/**
 * @file Matrix.hpp
 * @author cheny16
 * @brief Matrix header file
 * @version 1.0
 * @date 2023-03-24
 * 
 */
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <string>

typedef double* Matrix;
typedef double* Vector;

Matrix AllocMatrix(const int rows, const int cols);
Matrix ReadMatrix(const std::string filename, int & rows, int & cols);

Matrix AllocVector(const int rows);

void FreeMatrix(Matrix mat);
void FreeVector(Vector vec);

void PrintMatrix(const Matrix mat, const int rows, const int cols);


#endif /* MATRIX_HPP */