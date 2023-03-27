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

# Find CBLAS
option(USE_MKL "Use Intel MKL as the BLAS implementation" ON)

if(USE_MKL)
    set(BLA_VENDOR Intel10_64lp)
    add_definitions(-DUSE_MKL)
    message(STATUS "Using Intel MKL")
else()
    set(BLA_VENDOR OpenBLAS)
    message(STATUS "Using OpenBLAS")
endif()

find_package(BLAS)

if(USE_MKL)
    include_directories(${MKLROOT}/include)
endif()

# include
target_include_directories(Matrix PRIVATE ${HEADER_DIR})
target_include_directories(Givens PRIVATE ${HEADER_DIR} ${BLAS_INCLUDE_DIRS})
target_include_directories(GMRES  PRIVATE ${HEADER_DIR} ${BLAS_INCLUDE_DIRS})

# link
target_link_libraries(GMRES PUBLIC Matrix Givens PRIVATE m ${OpenMP_CXX_FLAGS} ${BLAS_LIBRARIES})

# dependencies
add_dependencies(GMRES Matrix Givens)