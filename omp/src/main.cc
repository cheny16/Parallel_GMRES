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

#include <boost/program_options.hpp>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <random>
#include <chrono>
#include <numeric>
#include <cblas.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace po = boost::program_options;

int main(int argc, char **argv) {
    std::string inputFile;  /* Input file */
    float tol;              /* Convergence tolerance */
    int maxiter;            /* Max iter */
    bool verbose;           /* If verbose */
    bool timing;            /* If show timing */

    // Parse arguments
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
            std::cerr << desc << "\n";
            return -1;
        }

        if (!vm.count("input")) {
            std::cerr << "No input file. Please specify input file.\n";
            return -1;
        }

        if (!vm.count("tol")) {
            std::cerr << "No tolerance value. Please specify tolerance.\n";
            return -1;
        }

        inputFile = vm["input"].as<std::string>();
        tol       = vm["tol"].as<float>();
        maxiter   = vm["maxiter"].as<int>();
        verbose   = vm["verbose"].as<bool>();
        timing    = vm["timing"].as<bool>();
    } 
    catch(std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return -1;
    }
    
    int rows;               /* Rows */
    int cols;               /* Cols */
    Matrix A = nullptr;     /* Matrix A */

    // Read matrix from file
    A = ReadMatrix(inputFile, rows, cols);

    if (verbose) {
        std::cerr << "Matrix A:\n";
        PrintMatrix(A, rows, rows);
    }

    Vector b = AllocVector(rows);     /* Vector b */
    Vector x = AllocVector(rows);     /* Initial guess of X (x0) */
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (auto i = 0; i < rows; i++) {
        b[i] = 1.0f;
        x[i] = 0.0f;
    }

    // Call GMRES
    GMRES(A, b, x, rows, tol, maxiter, verbose);

    std::cerr << "Solution:\n";
    PrintMatrix(x, rows, 1);

    return 0;
}
