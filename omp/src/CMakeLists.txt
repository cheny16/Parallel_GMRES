set(HEADER_DIR "../include")

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall")

add_library(Matrix Matrix.cc ${HEADER_DIR}/Matrix.hpp)
add_library(Givens Givens.cc ${HEADER_DIR}/Givens.hpp)
add_library(GMRES  GMRES.cc  ${HEADER_DIR}/GMRES.hpp)

# Find OpenMP
find_package(OpenMP REQUIRED)
if (NOT OpenMP_FOUND)
    message(FATAL_ERROR "[Error]: OpenMP required.\n")
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

# include
target_include_directories(Matrix PRIVATE ${HEADER_DIR})
target_include_directories(Givens PRIVATE ${HEADER_DIR} ${BLAS_INCLUDE_DIRS})
target_include_directories(GMRES  PRIVATE ${HEADER_DIR} ${BLAS_INCLUDE_DIRS})

# link
target_link_libraries(GMRES PUBLIC Matrix Givens PRIVATE m ${OpenMP_CXX_FLAGS})

# dependencies
add_dependencies(GMRES Matrix Givens)