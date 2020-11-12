cmake_minimum_required(VERSION 3.5.0)

project(wlcontroller)

# option(USEBLAS "Use BLAS and C instead of Eigen (C++)" OFF)
option(PYBIND "Build python bindings" ON)

set(PYBIND11_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../controls/thirdparty/pybind11" CACHE STRING "Pybind location")
set(EIGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../../SDK/thirdparty/eigen" CACHE STRING "Eigen location")

set(WLQP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/wlqp)
file(GLOB osqp_src "${WLQP_DIR}/*.c")
# emosqp options
set(EMBEDDED 2)
set(PRINTING OFF CACHE BOOL "")
set(PROFILING OFF CACHE BOOL "")
set(DFLOAT ON CACHE BOOL "")
set(DLONG OFF CACHE BOOL "")
set(QDLDL_FLOAT_TYPE "float")
set(QDLDL_INT_TYPE "int")
set(QDLDL_BOOL_TYPE "int")
# Generate header file with the global options
# configure_file(${WLQP_DIR}/configure/qdldl_types.h.in ${WLQP_DIR}/include/qdldl_types.h)
# configure_file(${WLQP_DIR}/configure/osqp_configure.h.in ${WLQP_DIR}/include/osqp_configure.h)
include_directories(${WLQP_DIR})

# STUPID C MATLAB VERSION -------------------------------------
find_package(Matlab)
set(ML_DIR ${CMAKE_CURRENT_SOURCE_DIR}/c)
set(ML_SRC 
  ${ML_DIR}/eigenc.c 
  ${ML_DIR}/funapprox.c
  ${ML_DIR}/wlqp.c
  ${ML_DIR}/wlcontroller.c)

add_definitions(-DBASIC_BLAS) # use this to build locally, but not on the target
link_directories(${Matlab_ROOT_DIR}/extern/lib/win64/microsoft)

# add_library(wlcontrollerMat STATIC ${ML_SRC} ${osqp_src})
# target_compile_options(wlcontrollerMat PRIVATE -DBASIC_BLAS)
# target_link_directories(wlcontrollerMat PUBLIC ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft)
# target_link_libraries(wlcontrollerMat libmwblas)
# target_include_directories(wlcontrollerMat PUBLIC ${Matlab_INCLUDE_DIRS})

# # add_executable(test ${ML_DIR}/main.c)
# # target_link_libraries(test PUBLIC wlcontrollerMat)

# matlab_add_mex(NAME wlControllerUpdate SRC wlmexc.c LINK_TO wlcontrollerMat)

# all at once
matlab_add_mex(NAME wlControllerUpdate SRC wlmexc.c ${ML_SRC} ${osqp_src} LINK_TO libmwblas)
target_include_directories(wlControllerUpdate PUBLIC ${ML_DIR})
target_compile_options(wlControllerUpdate PRIVATE -DBASIC_BLAS)

add_custom_command(TARGET wlControllerUpdate POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:wlControllerUpdate> ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Copied mex")
  
# C++ version ---------------------------------------------------------
include_directories(${EIGEN_DIR})
set(CPPSRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
set(SRC_CONTROLLER 
  ${CPPSRC_DIR}/wlqp.cpp
  ${CPPSRC_DIR}/wlcontroller.cpp
  )
set(SRC_CONTROLLER ${SRC_CONTROLLER} ${osqp_src})
add_library(wlcontroller STATIC ${SRC_CONTROLLER})

if(PYBIND)
  add_subdirectory(${PYBIND11_DIR} "pybind11")
  pybind11_add_module(wlqppy ${CMAKE_CURRENT_SOURCE_DIR}/pybind/wlqppy.cpp)
  target_link_libraries(wlqppy PRIVATE wlcontroller)

  add_custom_command(TARGET wlqppy POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:wlqppy> ${CMAKE_CURRENT_SOURCE_DIR}/..
    COMMENT "Copied python binding")
endif(PYBIND)
