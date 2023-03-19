/**
 * @file ParaQR.cc
 * @author cheny16
 * @brief Parallel QR
 * @version 0.1
 * @date 2023-03-19
 * 
 */
#include "../include/ParaQR.hpp"
#include "../include/Matrix.hpp"
#include <cblas.h>

/**
* @brief Generate Householder reflection matrix P = I - beta*VV'
 * 
 * @param V Vector
 * @param rows Rows of Vector
 * @return Matrix 
 */
Matrix GenHHMatrix(Vector V, int rows) {
    float * I = new float[rows*rows];
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            if (i == j) I[i * rows + j] = 1;
            else I[i * rows + j] = 0;
        }
    }
    cblas_ssyr(CblasRowMajor, CblasLower, rows, -2, V, 1, I, rows);

    Matrix P = allocate_matrix(rows, rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            if (i >= j) P[i][j] = I[i * rows + j];
            else P[i][j] = I[j * rows + i];
        }
    }

    delete [] I;
    return P;
}

/**
 * @brief QR Decomposition with Householder reflections
 * 
 * @param Q [In/Out] Matrix Q
 * @param R [In/Out] Matrix R
 * @param rows Rows of A
 * @param comm MPI_Comm
 */
void QRDecomp(Matrix Q, Matrix R, int rows, MPI_Comm comm) {
        int myrank;
    int mysize;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &mysize);

    if (myrank == 0) {
        /* Initialize Q(=I) */
        for (int i = 0 ; i < rows; i++) {
            for (int j = 0; j < rows; j++) {
                if (i == j)
                    Q[i][j] = 1.0f;
                else
                    Q[i][j] = 0.0f;
            }
        }
    }

    Matrix P = nullptr;
    Matrix P_part = nullptr;

    int start_[MAX_PROCS] = {0}, end_[MAX_PROCS] = {0};
    int disp_1[MAX_PROCS] = {0}, sendcounts_1[MAX_PROCS] = {0};
    int disp_2[MAX_PROCS] = {0}, sendcounts_2[MAX_PROCS] = {0};

    for (int i = 0; i < rows; i++) {
        float *vec = new float[rows - i];
        float vec_norm;
        P = allocate_matrix(rows-i, rows-i);

        if (myrank == 0) {
            int q = (rows - i) / mysize;

            for (int k = 0; k < mysize; k++) {
                start_[k] = k * q;
                end_[k] = (k + 1) * q;
            }
            end_[mysize - 1] = rows - i;

            for (int k = 0; k < mysize; k++) {
                sendcounts_1[k] = (end_[k] - start_[k]) * (rows - i);
                disp_1[k] = start_[k] * (rows - i);
                sendcounts_2[k] = (end_[k] - start_[k]) * rows;
                disp_2[k] = start_[k] * rows;
            }

            /* Cauculate vector */
            for (int j = i; j < rows; j++) {
                vec[j - i] = -R[j][i];
            }

            float x_norm = cblas_snrm2(rows-i, vec, 1);

            if (vec[0] < 0)
                vec[0] = vec[0] + x_norm;
            else
                vec[0] = vec[0] - x_norm;
            vec_norm = cblas_snrm2(rows-i, vec, 1);

            if (vec_norm > 0) {
                /* Normalize vector */
                for (int j = 0; j < rows-i; j++) {
                    vec[j] /= vec_norm;
                }
                /* Calculate Householder Matrix P */
                P = GenHHMatrix(vec, rows-i);
            }
        }

        MPI_Bcast(start_, mysize, MPI_INT, 0, comm);
        MPI_Bcast(end_, mysize, MPI_INT, 0, comm);
        MPI_Bcast(sendcounts_1, mysize, MPI_INT, 0, comm);
        MPI_Bcast(disp_1, mysize, MPI_INT, 0, comm);
        MPI_Bcast(sendcounts_2, mysize, MPI_INT, 0, comm);
        MPI_Bcast(disp_2, mysize, MPI_INT, 0, comm);
        MPI_Bcast(&vec_norm, 1, MPI_FLOAT, 0, comm);

        int sublines = end_[myrank] - start_[myrank];

        P_part = allocate_matrix(sublines, rows-i);

        MPI_Scatterv(&P[0][0], sendcounts_1, disp_1, MPI_FLOAT, &P_part[0][0], sendcounts_1[myrank], MPI_FLOAT, 0, comm);
        MPI_Bcast(&Q[0][0], rows * rows, MPI_FLOAT, 0, comm);
        MPI_Bcast(&R[0][0], rows * rows, MPI_FLOAT, 0, comm);
        
        if (vec_norm > 0) {
            //R
            Matrix R_sub = allocate_matrix(rows-i, rows-i);
            Matrix R_subpart = MMmultiply(P_part, R, 0, sublines, 0, rows - i, i, rows, i, rows);

            MPI_Gatherv(&R_subpart[0][0], sendcounts_1[myrank], MPI_FLOAT, &R_sub[0][0], sendcounts_1, disp_1, MPI_FLOAT, 0, comm);

            /* Update original R */
            if (myrank == 0) {
                for (int k = 0; k < rows - i; k++) {
                    for (int l = 0; l < rows - i; l++) {
                        R[k + i][l + i] = R_sub[k][l];
                    }
                }
            }

            //Q
            Matrix Q_sub = allocate_matrix(rows-i, rows);

            Matrix Q_subpart = MMmultiply(P_part, Q, 0, sublines, 0, rows - i, i, rows, 0, rows);

            MPI_Gatherv(&Q_subpart[0][0], sendcounts_2[myrank], MPI_FLOAT, &Q_sub[0][0], sendcounts_2, disp_2, MPI_FLOAT, 0, comm);

            /* Update original Q */
            if (myrank == 0) {
                for (int k = 0; k < rows - i; k++) {
                    for (int l = 0; l < rows; l++) {
                        Q[k + i][l] = Q_sub[k][l];
                    }
                }
            }

            dealloc_matrix(R_sub);
            dealloc_matrix(Q_sub);
        }
    }
    if (P != nullptr) dealloc_matrix(P);
}