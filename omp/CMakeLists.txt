cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME gmres)
set(SOURCE_DIR src)

# Project name and other infomations
project(${PROJECT_NAME} VERSION 1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall")

# Find CBLAS
option(USE_MKL "Use Intel MKL as the BLAS implementation" OFF)

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

SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREADED OFF)

# Find Boost package for program_options
find_package(Boost 1.32.0 REQUIRED COMPONENTS program_options)
if (NOT Boost_FOUND)
    message(FATAL_ERROR "[Fatal error]: Boost (version >= 1.32.0) required.\n")
endif()

# Find MPI for parallel program
find_package(OpenMP REQUIRED)
if (NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "[Error]: OpenMP required.\n")
endif()

# Add the executable
add_executable(${PROJECT_NAME} ${SOURCE_DIR}/main.cc)
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}" "${SOURCE_DIR}")

target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(${PROJECT_NAME} PUBLIC 
                        Matrix 
                        Givens 
                        GMRES 
                        Boost::program_options 
                        OpenMP::OpenMP_CXX 
                        ${BLAS_LIBRARIES})

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
