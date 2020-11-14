# =================================================================
#
# Utility macros for CMake
#  Changxi Zheng (Columbia University)
# 
# =================================================================
macro(config_compiler_and_linker)
    if (CMAKE_CXX_COMPILER MATCHES ".*icpc$")
        if ( NOT USE_DEBUG )
            #set(CMAKE_CXX_FLAGS_RELEASE "-O3 -no-prec-div -xHost -opt-prefetch -unroll-aggressive -DNDEBUG")
            set(CMAKE_CXX_FLAGS_RELEASE "-O3 -no-prec-div -unroll-aggressive -DNDEBUG")
        endif ()

        set(OPENMP_LIB iomp5)
    endif ()
endmacro()

#macro(disable_eigen_warnings)
#    set(EIGEN_DISABLE_WARNINGS 1)
#    add_definitions(-wd161)
#endmacro()

macro(enable_cpp11)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG(-std=c++0x COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    elseif(COMPILER_SUPPORTS_CXX0X)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
endmacro()

macro(enable_cpp14)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endmacro()

macro(enable_cpp17)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endmacro()

# Macros to build the targets
include( CMakeParseArguments )

# Usage:
# 
# add_exe( name source1 [source2 ...]
#          [LINK_LIBRARIES lib1 ...]
#          [INCLUDE_DIRS dir1 ...]
#          [OUT_DIR "output dir"] )
function(add_exe name)
    cmake_parse_arguments(_exe "" "OUT_DIR" "LINK_LIBRARIES;INCLUDE_DIRS" ${ARGN})
    set(_exe_srcs ${_exe_UNPARSED_ARGUMENTS})
    message(STATUS "Add a executable target ${name}")

    add_executable(${name} ${_exe_srcs})
    if (_exe_OUT_DIR)
        set_target_properties(${name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY ${_exe_OUT_DIR})
    endif()

    if (_exe_LINK_LIBRARIES)
        target_link_libraries(${name} ${_exe_LINK_LIBRARIES})
    endif()

    if (_exe_INCLUDE_DIRS)
        include_directories(${_exe_INCLUDE_DIRS})
    endif()
endfunction()
