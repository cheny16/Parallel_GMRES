/**
 * @file Matrix.hpp
 * @author cheny16
 * @brief Matrix header
 * @version 0.1
 * @date 2023-03-17
 * 
 */
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <string>

typedef float** Matrix;
typedef float*  Vector;

Matrix allocate_matrix(int rows, int cols);
Matrix read_matrix(const std::string filename, int & rows, int & cols);

void dealloc_matrix(Matrix mat);

void print_matrix(Matrix mat, int rows, int cols);

Matrix MMmultiply(const Matrix A, const Matrix B, 
                  const int A_rows_start, const int A_rows_end, const int A_cols_start, const int A_cols_end, 
                  const int B_rows_start, const int B_rows_end, const int B_cols_start, const int B_cols_end, 
                  bool if_a_trans=false,  bool if_b_trans=false);

Vector MVmultiply(const Matrix A, const Vector B, 
                  const int A_rows, const int A_cols, const int B_rows);

#endif /* MATRIX_HPP */