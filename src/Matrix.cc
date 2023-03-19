/**
 * @file Matrix.cc
 * @author cheny16
 * @brief Matrix interfaces
 * @version 0.1
 * @date 2023-03-19
 * 
 */
#include "../include/Matrix.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cblas.h>

/**
 * @brief Allocate matrix
 * 
 * @param rows Rows
 * @param cols Cols
 * @return Matrix
 */
Matrix allocate_matrix(int rows, int cols) {
    float * data = new float[rows * cols];
    float ** mat = new float * [rows];
    for (int i = 0; i < rows; i++) {
        mat[i] = data + i*cols;
    }
    return mat;
}

/**
 * @brief Read matrix from input file
 * 
 * @param filename Filename of SuiteSparse Matrix file (.mtx)
 * @param rows Rows
 * @param cols Cols
 * @return Matrix 
 */
Matrix read_matrix(const std::string filename, int & rows, int & cols) {
    std::ifstream file;
    file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try {
        file.open(filename);
    }
    catch (std::ifstream::failure e) {
        std::cerr << "Error opening file: " << e.what() << "\n";
    }

    int num_rows, num_cols;
    int lines;

    file >> num_rows >> num_cols >> lines;
    if (num_rows == num_cols) {
        rows = num_rows;
        cols = num_cols;
    }
    
    int i, j;
    float value;
    Matrix mat = allocate_matrix(rows, cols);

    for(auto line = 0; line < lines; line++){
        file >> i >> j >> value;
        if (i < rows-1 && j < cols-1) {
            mat[i][j] = value;
        }
    }
    file.close();
    return mat;
}

/**
 * @brief De-alloc matrix
 * 
 * @param mat Matrix
 * @param rows Rows
 * @param cols Cols
 */
void dealloc_matrix(Matrix mat) {
    delete [] mat[0];
    delete [] mat;
}

void print_matrix(Matrix mat, int rows, int cols) {
    for (int i = 0; i < static_cast<int>(rows); i++) {
        for (int j = 0; j < static_cast<int>(cols); j++) {
            std::cout << mat[i][j] << std::setw(2);
        }
        std::cout << "\n";
    }
}

/**
 * @brief Multiply two matrix
 * 
 * @param A Matrix A
 * @param B Matrix B
 * @param A_rows_start Start of A rows
 * @param A_rows_end End of A rows
 * @param A_cols_start Start of A cols
 * @param A_cols_end End of A cols
 * @param B_rows_start Start of B rows
 * @param B_rows_end End of B rows
 * @param B_cols_start Start of B cols
 * @param B_cols_end End of B cols
 * @param if_a_trans If TRANSPOSE A
 * @param if_b_trans If TRANSPOSE B
 * @return Matrix C = A*B
 */
Matrix MMmultiply(const Matrix A, const Matrix B, 
                  const int A_rows_start, const int A_rows_end, const int A_cols_start, const int A_cols_end, 
                  const int B_rows_start, const int B_rows_end, const int B_cols_start, const int B_cols_end, 
                  bool if_a_trans,  bool if_b_trans) {
    int A_rows = A_rows_end - A_rows_start;
    int A_cols = A_cols_end - A_cols_start;
    int B_rows = B_rows_end - B_rows_start;
    int B_cols = B_cols_end - B_cols_start;

    if (A_cols != B_rows) {
        std::cerr << "[ERROR] Dimensions of two matrix don't match!\n";
        exit(EXIT_FAILURE);
    }

    float *aa = new float[A_rows*A_cols];
    float *bb = new float[B_rows*B_cols];
    float *cc = new float[A_rows*B_cols];

    for (int i = 0; i < A_rows; i++) {
        for (int j = 0; j < A_cols; j++) {
            aa[i*A_cols + j] = A[A_rows_start + i][A_cols_start + j];
        }
    }
    for (int i = 0; i < B_rows; i++) {
        for (int j = 0; j < B_cols; j++) {
            bb[i*B_cols + j] = B[B_rows_start + i][B_cols_start + j];
        }
    }

    CBLAS_TRANSPOSE aTrans = {CblasNoTrans}, bTrans = {CblasNoTrans};
    if (if_a_trans) aTrans = CblasTrans;
    if (if_b_trans) bTrans = CblasTrans;

    cblas_sgemm(CblasRowMajor, aTrans, bTrans, A_rows, B_cols, A_cols, 1, aa, A_cols, bb, B_cols, 0, cc, B_cols);

    Matrix C = allocate_matrix(A_rows, B_cols);
    for (int i = 0; i < A_rows; i++) {
        for (int j = 0; j < B_cols; j++) {
            C[i][j] = cc[i * B_cols + j];
        }
    }

    delete [] aa; delete [] bb; delete [] cc;

    return C;
}

/**
 * @brief Matrix Vector Multiply
 * 
 * @param A Matrix A
 * @param B Vector B
 * @param A_rows Rows of A
 * @param A_cols Cols of A
 * @param B_rows Rows of B
 * @return Vector Vector C = A*B
 */
Vector MVmultiply(const Matrix A, const Vector B, 
                  const int A_rows, const int A_cols, const int B_rows) {
    if (A_cols != B_rows) {
        std::cerr << "[ERROR] Dimensions of two matrix don't match!\n";
        exit(EXIT_FAILURE);
    }

    float *aa = new float[A_rows*A_cols];

    for (int i = 0; i < A_rows; i++) {
        for (int j = 0; j < A_cols; j++) {
            aa[i*A_cols + j] = A[i][j];
        }
    }

    Vector C = new float[A_rows];
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, 1, A_cols, 1, aa, A_cols, B, 1, 0, C, 1);

    delete [] aa;

    return C;
}