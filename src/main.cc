/**
 * @file main.cc
 * @author cheny16
 * @brief Main function
 * @version 0.1
 * @date 2023-03-19
 * 
 */
#include "../include/Matrix.hpp"
#include "../include/GMRES.hpp"
#include "../include/ParaQR.hpp"

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <random>
#include <chrono>
#include <numeric>
#include <cblas.h>
#include <mpi.h>

namespace po = boost::program_options;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int worldRank, worldSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldSize);

    std::string inputFile;  /* Input file */
    float tol;              /* Convergence tolerance */
    int maxiter;            /* Max iter */
    bool verbose;           /* If verbose */
    bool timing;            /* If show timing */

    // Parse arguments
    if (worldRank == 0) {
        try {
            po::options_description desc("Allowed options");
            desc.add_options()
                ("help", "show help message")
                ("input", po::value<std::string>(), "input matrix file (.mtx)")
                ("tol", po::value<float>(), "convergence tolerance")
                ("maxiter", po::value<int>(), "max iterations")
                ("verbose", po::bool_switch()->default_value(false), "display details")
                ("timing", po::bool_switch()->default_value(false), "display timings");

            po::variables_map vm;
            po::store(po::parse_command_line(argc, argv, desc), vm);
            if (vm.count("help")) {
                std::cout << desc << "\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            if (!vm.count("input")) {
                std::cerr << "No input file. Please specify input file.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            if (!vm.count("tol")) {
                std::cerr << "No tolerance value. Please specify tolerance.\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            inputFile = vm["input"].as<std::string>();
            tol       = vm["tol"].as<float>();
            maxiter   = vm["maxiter"].as<int>();
            verbose   = vm["verbose"].as<bool>();
            timing    = vm["timing"].as<bool>();
        } 
        catch(std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    
    int rows;               /* Rows */
    int cols;               /* Cols */
    Matrix matA = nullptr;  /* Matrix A */

    // Read matrix from file
    if (worldRank == 0) {
        matA = read_matrix(inputFile, rows, cols);
    }

    // Bcast rows & cols
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

    Vector vecB = new float[rows];    /* Vector b */
    Vector X    = new float[rows];    /* Initial guess of X (x0) */
    for (auto i = 0; i < rows; ++i) {
        vecB[i] = 1.0f;
        X[i]    = 1.0f;
    }

    // Bcast matrix
    if (worldRank != 0) {
        matA = allocate_matrix(rows, cols);
    }
    MPI_Bcast(&matA[0][0], rows*rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // GMRES
    GMRES(matA, vecB, X, rows, tol, maxiter, MPI_COMM_WORLD, verbose);

    return 0;
}