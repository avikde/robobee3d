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

set(SRC_CONTROLLER 
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/wlqp.cpp
  )
set(SRC_CONTROLLER ${SRC_CONTROLLER} ${osqp_src})

add_executable(test ${SRC_CONTROLLER})
