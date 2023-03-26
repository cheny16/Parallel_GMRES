/**
 * @file Matrix.cc
 * @author cheny16
 * @brief Matrix interfaces ("Class-like")
 * @version 1.0
 * @date 2023-03-24
 * 
 */
#include "../include/Matrix.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cblas.h>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * @brief Allocate matrix
 * 
 * @param rows Rows
 * @param cols Cols
 * @return Matrix
 */
Matrix AllocMatrix(const int rows, const int cols) {
    double *mat = new double[rows * cols];
    if (mat == nullptr) {
        std::cerr << "Error allocating matrix!\n";
        exit(EXIT_FAILURE);
    }
#ifdef _OPENMP
    #pragma omp parallel for simd
#endif
    for (auto i = 0; i < rows*cols; i++) {
        mat[i] = 0.0f;
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
Matrix ReadMatrix(const std::string filename, int & rows, int & cols) {
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
    double value;
    Matrix mat = AllocMatrix(rows, cols);

    for(auto line = 0; line < lines; line++){
        file >> i >> j >> value;
        try {
            mat[(i-1) * cols + (j-1)] = value;
        }
        catch (const std::out_of_range oor) {
            std::cerr << "Out of Range error: " << oor.what() << "\n";
        }

    }
    file.close();
    return mat;
}

/**
 * @brief Allocate vector
 * 
 * @param rows Rows
 * @param cols Cols
 * @return Vector
 */
Matrix AllocVector(const int rows) {
    double *vec = new double[rows];
    if (vec == nullptr) {
        std::cerr << "Error allocating vector!\n";
        exit(EXIT_FAILURE);
    }
#ifdef _OPENMP
    #pragma omp parallel for simd
#endif
    for (auto i = 0; i < rows; i++) {
        vec[i] = 0.0f;
    }
    return vec;
}

/**
 * @brief Free matrix
 * 
 * @param mat Matrix
 */
void FreeMatrix(Matrix mat) {
    delete [] mat;
}

/**
 * @brief Free vector
 * 
 * @param vec Vector
 */
void FreeVector(Vector vec) {
    delete [] vec;
}

/**
 * @brief Print matrix
 * 
 * @param mat Matrix
 * @param rows Rows
 * @param cols Cols
 */
void PrintMatrix(const Matrix mat, const int rows, const int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            try {
                fprintf(stderr, "%6.2f", mat[i * cols + j]);
            }
            catch (const std::out_of_range oor) {
                std::cerr << "Out of Range error: " << oor.what() << "\n";
            }
        }
        fprintf(stderr, "\n");
    }
}

