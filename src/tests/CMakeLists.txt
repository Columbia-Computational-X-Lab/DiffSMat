 include_directories(
#     ${TBB_INCLUDE_DIR}
     ${MKL_INCLUDE_DIR}
     ${EIGEN3_INCLUDE_DIRS}
     )

add_exe(test_dcomplex_eigen
        test_dcomplex_eigen.cpp
        ${ROOT_SRC_DIR}/dcomplex/DComplex.hpp
    LINK_LIBRARIES
        ${MKL_LIBRARIES}
        ${TBB_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")

add_exe(test_dcomplex_deriv
        test_dcomplex_deriv.cpp
        ${ROOT_SRC_DIR}/dcomplex/DComplex.hpp
    LINK_LIBRARIES
        ${TBB_LIBRARIES}
        ${MKL_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")

add_exe(test_mat_exp
        test_mat_exp.cpp
        ${ROOT_SRC_DIR}/dcomplex/DComplex.hpp
        ${ROOT_SRC_DIR}/dcomplex/DComplexMath.hpp
    LINK_LIBRARIES
        ${MKL_LIBRARIES}
        ${TBB_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")

file(GLOB RCWA_SRC_FILES 
    "${ROOT_SRC_DIR}/rcwa/*.h"
    "${ROOT_SRC_DIR}/rcwa/*.cpp")

add_exe(test_rcwa
        test_rcwa.cpp
        ${RCWA_SRC_FILES}
        ${ROOT_SRC_DIR}/core/RCWAScatterMatrix.h
        ${ROOT_SRC_DIR}/core/RCWAScatterMatrix.cpp
    LINK_LIBRARIES
        ${MKL_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")

add_exe(test_mat_deriv
        test_matrix_deriv.cpp
        ${ROOT_SRC_DIR}/core/MatrixDerivative.hpp
    LINK_LIBRARIES
        ${TBB_LIBRARIES}
        ${MKL_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")

add_exe(test_smat_deriv
        test_smat_deriv.cpp
        ${ROOT_SRC_DIR}/core/MatrixDerivative.hpp
        ${ROOT_SRC_DIR}/core/ScatterMatrixDerivative.hpp
        ${ROOT_SRC_DIR}/core/RCWAScatterMatrix.h
        ${ROOT_SRC_DIR}/core/RCWAScatterMatrix.cpp
        ${RCWA_SRC_FILES}
    LINK_LIBRARIES
        ${MKL_LIBRARIES}
        ${TBB_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")

add_exe(test_mat_product
        test_mat_product.cpp
        ${ROOT_SRC_DIR}/dcomplex/DComplex.hpp
        ${ROOT_SRC_DIR}/dcomplex/DComplexMatrixFunc.hpp
    LINK_LIBRARIES
        ${MKL_LIBRARIES}
        ${TBB_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")

add_exe(test_mat_exp_deriv
        test_mat_exp_deriv.cpp
        ${ROOT_SRC_DIR}/dcomplex/DComplex.hpp
        ${ROOT_SRC_DIR}/dcomplex/DComplexMatrixFunc.hpp
        ${ROOT_SRC_DIR}/core/MatrixDerivative.hpp
    LINK_LIBRARIES
        ${MKL_LIBRARIES}
        ${TBB_LIBRARIES}
    OUT_DIR "${ROOT_TESTS_DIR}")