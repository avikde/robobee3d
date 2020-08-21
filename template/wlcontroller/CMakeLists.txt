cmake_minimum_required(VERSION 3.5.0)

project(wlcontroller)

set(WLQP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/wlqp)
file(GLOB osqp_src "${WLQP_DIR}/src/osqp/*.c")
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
configure_file(${WLQP_DIR}/configure/qdldl_types.h.in ${WLQP_DIR}/include/qdldl_types.h)
configure_file(${WLQP_DIR}/configure/osqp_configure.h.in ${WLQP_DIR}/include/osqp_configure.h)
include_directories(${WLQP_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../SDK/thirdparty/eigen)

# # C++ version
# set(CPPSRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cpp)
# set(SRC_CONTROLLER 
#   ${CPPSRC_DIR}/wlqp.cpp
#   ${CPPSRC_DIR}/wlcontroller.cpp
#   )
# set(SRC_CONTROLLER ${SRC_CONTROLLER} ${osqp_src})
# add_library(wlcontroller STATIC ${SRC_CONTROLLER})

# add_executable(test ${CPPSRC_DIR}/main.cpp)
# target_link_libraries(test PUBLIC wlcontroller)

# STUPID C MATLAB VERSION
find_package(Matlab)
message("HI ${Matlab_ROOT_DIR}")
set(ML_DIR ${CMAKE_CURRENT_SOURCE_DIR}/c)
set(ML_SRC 
  ${ML_DIR}/eigenc.c 
  ${ML_DIR}/funapprox.c
  ${ML_DIR}/wlqp.c
  ${ML_DIR}/wlcontroller.c)
add_library(wlcontrollerMat STATIC ${ML_SRC} ${osqp_src})
target_link_directories(wlcontrollerMat PUBLIC ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft)
target_link_libraries(wlcontrollerMat libmwblas)
target_include_directories(wlcontrollerMat PUBLIC ${Matlab_INCLUDE_DIRS})

add_executable(test ${ML_DIR}/main.c)
target_link_libraries(test PUBLIC wlcontrollerMat)

# # set(SRC_MEX ${CMAKE_CURRENT_SOURCE_DIR}/wlmex.cpp)
# matlab_add_mex(NAME wlmex SRC ${SRC_MEX} LINK_TO wlcontroller)

