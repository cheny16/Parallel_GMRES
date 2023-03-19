cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME ParallelGMRES)
set(SOURCE_DIR src)

# Project name and other infomations
project(${PROJECT_NAME} VERSION 1.0)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O0 -g -lblas")

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_BUILD_TYPE RelWithDebInfo)

# Find BLAS-LAPACK
find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    find_path(BLAS_INCLUDE_DIRS cblas.h
                                /usr/include
                                /usr/local/include
                                $ENV{CMAKE_PREFIX_PATH})
    message(STATUS BLAS found: ${BLAS_INCLUDE_DIRS})
endif()

SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREADED OFF)

# Find Boost package for program_options
find_package(Boost 1.32.0 REQUIRED COMPONENTS program_options)
if (NOT Boost_FOUND)
    message(FATAL_ERROR "[Fatal error]: Boost (version >= 1.32.0) required.\n")
endif()
# Find MPI for parallel program
find_package(MPI REQUIRED)
if (NOT MPI_FOUND)
    message(FATAL_ERROR "[Fatal error]: MPI required.\n")
endif()

# Add the executable
add_executable(${PROJECT_NAME} ${SOURCE_DIR}/main.cc)
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_CXX_INCLUDE_PATH}
                           ${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}" "${SOURCE_DIR}")

target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_link_libraries(${PROJECT_NAME} PUBLIC 
                        Matrix 
                        ParaQR 
                        GMRES 
                        Boost::program_options 
                        ${MPI_CXX_LIBRARIES}
                        ${MPI_CXX_LINK_FLAGS} 
                        ${lapackblas_libraries})

# Add modules
link_directories(${SOURCE_DIR})
add_subdirectory(${SOURCE_DIR})

# Add clean target
add_custom_target(distclean COMMAND rm -rf CMakeCache.txt CMakeFiles
                  Makefile cmake_install.cmake CTestTestfile.cmake Testing)

# option(BUILD_TESTING "Build testing tree." OFF)
# include(CTest)

# if(BUILD_TESTING)
#     add_subdirectory(tests)
#     file(COPY ./datas DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/tests/)
# endif()
