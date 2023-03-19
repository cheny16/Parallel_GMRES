set(HEADER_DIR "../include")

add_library(Matrix Matrix.cc ${HEADER_DIR}/Matrix.hpp)
add_library(ParaQR ParaQR.cc ${HEADER_DIR}/ParaQR.hpp)
add_library(GMRES  GMRES.cc  ${HEADER_DIR}/GMRES.hpp)

# Find MPI
find_package(MPI REQUIRED)
if (NOT MPI_FOUND)
    message(FATAL_ERROR "[Fatal error]: MPI required.\n")
endif()

# Find BLAS-LAPACK
find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
    find_path(BLAS_INCLUDE_DIRS cblas.h
                                /usr/include
                                /usr/local/include
                                $ENV{CMAKE_PREFIX_PATH})
    message(STATUS BLAS found: ${BLAS_INCLUDE_DIRS})
endif()

target_include_directories(Matrix PRIVATE ${HEADER_DIR} ${MPI_CXX_INCLUDE_PATH} ${BLAS_INCLUDE_DIRS})
target_include_directories(ParaQR PRIVATE ${MPI_CXX_INCLUDE_PATH}
                                  PUBLIC ${HEADER_DIR} ${BLAS_INCLUDE_DIRS})
target_include_directories(GMRES  PRIVATE ${MPI_CXX_INCLUDE_PATH}
                                  PUBLIC ${HEADER_DIR} ${BLAS_INCLUDE_DIRS})
# target_link_libraries(ParaQR PUBLIC Matrix)
# target_link_libraries(GMRES PUBLIC Matrix ParaQR)
