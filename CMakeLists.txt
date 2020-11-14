cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

#
# Differentiable Scattering Matrix
# Developed by 
# Changxi Zheng @ Columbia University
# Ziwei Zhu @ Columbia University
#

project(DifferentiableScatteringMatrix)

# General Configuration: find all the required libraries.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(sploosh_utils)
# compile option
option(BUILD_DEBUG "Turn on the debug mode" OFF)
option(BUILD_TESTS "Build unit test cases" ON)
option(USE_FAST_DCOMPLEX_PRODUCT "Build unit test cases" OFF)
option(USE_DCOMPLEX_DERIVATIVE "Use complex perturbation to compute matrix exponential derivative" OFF)

#===================================================================
## Compiler
# set compiler flags for debug/release
if ( BUILD_DEBUG )
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
else ()
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()

if ( USE_FAST_DCOMPLEX_PRODUCT )
    message(WARNING "Use fast (approximated) DComplex multiplication.")
    add_definitions(-DFAST_DCOMPLEX_PRODUCT)
endif ()

if ( USE_DCOMPLEX_DERIVATIVE )
    message(WARNING "Use DComplex perturbation to compute matrix exponential derivative.")
    add_definitions(-DUSE_DCOMPLEX_DERIVATIVE)
endif ()

add_definitions(-Wall)

config_compiler_and_linker()

# Try to use C++14 Standards
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++17 COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG(-std=c++14 COMPILER_SUPPORTS_CXX14)
if (COMPILER_SUPPORTS_CXX17)
    message(STATUS "Use C++17 standard")
    enable_cpp17()
    add_definitions(-DUSE_CXX17)
elseif (COMPILER_SUPPORTS_CXX14)
    message(STATUS "Use C++14 standard")
    enable_cpp14()
    add_definitions(-DUSE_CXX14)
else()
    message(STATUS "Use C++11 standard")
    enable_cpp11()
endif()

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang$")
    add_definitions(-Wno-unused-local-typedef)
endif()

set(ROOT_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE STRING "Root build directory")
set(ROOT_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" CACHE STRING "Root source directory")
set(ROOT_BIN_DIR "${ROOT_BUILD_DIR}/bin" CACHE STRING "Root binary directory")
set(ROOT_TESTS_DIR "${ROOT_BUILD_DIR}/tests" CACHE STRING "Root tests directory")
mark_as_advanced(ROOT_BUILD_DIR)
mark_as_advanced(ROOT_SRC_DIR)
mark_as_advanced(ROOT_BIN_DIR)
mark_as_advanced(ROOT_TESTS_DIR)

#===================================================================
find_package(TBB 2.2 REQUIRED)
if (TBB_FOUND)
    message(STATUS "Found Intel TBB version ${TBB_VERSION}: ${TBB_INCLUDE_DIRS}")
endif ()

find_package(MKL REQUIRED)
if(MKL_FOUND)
    message(STATUS "MKL libraries: ${MKL_LIBRARIES}")
    message(STATUS "MKL include directory: ${MKL_INCLUDE_DIR}")
    message(STATUS "MKL OpenMP type: ${MKL_OPENMP_TYPE}")
    message(STATUS "MKL OpenMP library: ${MKL_OPENMP_LIBRARY}")
    add_definitions(-DEIGEN_USE_MKL_ALL)
endif()

find_package(Eigen3 3.3 CONFIG REQUIRED)
if (EIGEN3_FOUND)
   message(STATUS "Found Eigen version ${EIGEN3_VERSION_STRING}: ${EIGEN3_INCLUDE_DIRS}")
endif (EIGEN3_FOUND)

# find_package(CBLAS REQUIRED)
# if ( CBLAS_Accelerate_LIBRARY )
#     # this is possibly available on Mac
#     set(NUMERIC_LIBS ${CBLAS_Accelerate_LIBRARY})
#     add_definitions(-DUSE_MAC_VECLIB)
# elseif ( CBLAS_cblas_LIBRARY )
#     include_directories(${CBLAS_cblas_f77blas_atlas_INCLUDE})
#     set(NUMERIC_LIBS ${CBLAS_cblas_LIBRARY})
# else()
#     message(SEND_ERROR "Fail to find usable CBLAS library")
# endif()

#===================================================================

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}/src)

#add_subdirectory(src/tools)

if ( BUILD_TESTS )
    add_subdirectory(src/tests)
endif ()
