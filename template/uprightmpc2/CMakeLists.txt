cmake_minimum_required(VERSION 3.5.0)

project(uprightmpc2)

set(PYBIND11_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../controls/thirdparty/pybind11" CACHE STRING "Pybind location")
set(EIGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../SDK/thirdparty/eigen" CACHE STRING "Eigen location")

set(BPYBIND ON CACHE BOOL "")

set(OSQP_DIR ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB osqp_src "${OSQP_DIR}/*.c")
# emosqp options
set(EMBEDDED 2)
set(PRINTING OFF CACHE BOOL "")
set(PROFILING OFF CACHE BOOL "")
set(DFLOAT ON CACHE BOOL "")
set(DLONG OFF CACHE BOOL "")
set(QDLDL_FLOAT_TYPE "float")
set(QDLDL_INT_TYPE "int")
set(QDLDL_BOOL_TYPE "int")
include_directories(${OSQP_DIR})
# No eigen or python or anything dependency
add_library(uprightmpc2 STATIC ${osqp_src})

# # STUPID C MATLAB VERSION -------------------------------------
# find_package(Matlab)
# set(ML_DIR ${CMAKE_CURRENT_SOURCE_DIR}/c)
# set(ML_SRC 
#   ${ML_DIR}/eigenc.c 
#   ${ML_DIR}/funapprox.c
#   ${ML_DIR}/wlqp.c
#   ${ML_DIR}/wlcontroller.c)

# add_definitions(-DBASIC_BLAS) # use this to build locally, but not on the target
# link_directories(${Matlab_ROOT_DIR}/extern/lib/win64/microsoft)

# # all at once
# matlab_add_mex(NAME wlControllerUpdate SRC wlmexc.c ${ML_SRC} ${osqp_src} LINK_TO libmwblas)
# target_include_directories(wlControllerUpdate PUBLIC ${ML_DIR})
# target_compile_options(wlControllerUpdate PRIVATE -DBASIC_BLAS)

# add_custom_command(TARGET wlControllerUpdate POST_BUILD
#   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:wlControllerUpdate> ${CMAKE_CURRENT_SOURCE_DIR}
#   COMMENT "Copied mex")
  
# C++ version ---------------------------------------------------------
if (${BPYBIND})
  include_directories(${EIGEN_DIR})
  add_subdirectory(${PYBIND11_DIR} "pybind11")
  pybind11_add_module(uprightmpc2py ${CMAKE_CURRENT_SOURCE_DIR}/py/uprightmpc2py.cpp ${CMAKE_CURRENT_SOURCE_DIR}/py/ematmult.cpp)
  target_link_libraries(uprightmpc2py PRIVATE uprightmpc2)

  add_custom_command(TARGET uprightmpc2py POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:uprightmpc2py> ${CMAKE_CURRENT_SOURCE_DIR}/..
    COMMENT "Copied python binding")
endif()
