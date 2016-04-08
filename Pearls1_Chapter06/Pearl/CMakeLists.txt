###############################################################################
# Copyright (c) 2014, Jefferson Amstutz
# Copyright (c) 2014, SURVICE Engineering
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
#   * Neither the name of SURVICE Engineering nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
###############################################################################

project(ftevaluator)
cmake_minimum_required(VERSION 2.8)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

##############################################################
# ISA configuration
##############################################################

SET(XEON_ISA "AVX" CACHE STRING "Select highest ISA to support (SSE2, SSE3, SSSE3, SSE4.1, SSE4.2, AVX, AVX-I, AVX2")
SET(XEON_PHI_ISA OFF CACHE BOOL "Generate code path for Xeon Phi")

IF (XEON_ISA STREQUAL "SSE2")
  SET(ISA  1)
ENDIF ()

IF (XEON_ISA STREQUAL "SSE3")
  SET(ISA  2)
ENDIF ()

IF (XEON_ISA STREQUAL "SSSE3")
  SET(ISA  3)
ENDIF ()

IF (XEON_ISA STREQUAL "SSE4.1")
  SET(ISA  4)
ENDIF ()

IF (XEON_ISA STREQUAL "SSE4.2")
  SET(ISA  5)
ENDIF ()

IF (XEON_ISA STREQUAL "AVX")
  SET(ISA  6)
ENDIF ()

IF (XEON_ISA STREQUAL "AVX-I")
  SET(ISA  7)
ENDIF ()

IF (XEON_ISA STREQUAL "AVX2")
  SET(ISA  8)
ENDIF ()

SET(TARGET_SSE2  OFF)
SET(TARGET_SSE3  OFF)
SET(TARGET_SSSE3  OFF)
SET(TARGET_SSE41  OFF)
SET(TARGET_SSE42  OFF)
SET(TARGET_AVX  OFF)
SET(TARGET_AVXI  OFF)
SET(TARGET_AVX2  OFF)

IF (ISA GREATER 0)
  SET(TARGET_SSE2  ON)
  SET(ISPC_TARGETS "sse2")
ENDIF ()

IF (ISA GREATER 1)
  SET(TARGET_SSE3  ON)
ENDIF ()

IF (ISA GREATER 2)
  SET(TARGET_SSSE3  ON)
ENDIF ()

IF (ISA GREATER 3)
  SET(TARGET_SSE41  ON)
ENDIF ()

IF (ISA GREATER 4)
  SET(TARGET_SSE42  ON)
  SET(ISPC_TARGETS ${ISPC_TARGETS} "sse4")
ENDIF ()

IF (ISA GREATER 5)
  SET(TARGET_AVX  ON)
  SET(ISPC_TARGETS ${ISPC_TARGETS} "avx")
ENDIF ()

IF (ISA GREATER 6)
  SET(TARGET_AVXI  ON)
ENDIF ()

IF (ISA GREATER 7)
  SET(TARGET_AVX2  ON)
  SET(ISPC_TARGETS ${ISPC_TARGETS} "avx2")
ENDIF ()

SET(ISPC_TARGETS "${ISPC_TARGETS}")
STRING(REGEX REPLACE ";" "," ISPC_TARGETS "${ISPC_TARGETS}")

##############################################################
# Compiler
##############################################################

SET(COMPILER "GCC" CACHE STRING "Set to GCC, CLANG, or ICC")

IF (XEON_PHI_ISA)
  INCLUDE (icc_xeonphi)
ELSEIF (COMPILER STREQUAL "GCC")
  INCLUDE (gcc)
ELSEIF (COMPILER STREQUAL "CLANG")
  INCLUDE (clang)
ELSEIF (COMPILER STREQUAL "ICC")
  INCLUDE (icc)
ELSE ()
  MESSAGE(FATAL_ERROR "Unknown compiler specified: " ${COMPILER})
ENDIF ()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type: Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_VERBOSE_MAKEFILE false CACHE BOOL "Enables verbose mode.")

IF(XEON_PHI_ISA)
  SET(__XEON__ 0)
ELSE()
  SET(__XEON__ 1)
ENDIF()

INCLUDE (ispc)

##############################################################
# Application
##############################################################

SET(FTE_THREADING ON CACHE BOOL "Enable OpenMP threading for evaluations")

if(FTE_THREADING)
  add_definitions(-DEVAL_THREADING)
endif()

SET(FTE_ISPC_PACKET_SIZE 16 CACHE INT "ISPC evaluation packet size")

add_definitions(-DPACKET_SIZE=${FTE_ISPC_PACKET_SIZE})

SET(FTE_ISPC_EVAL_ENABLED ON CACHE BOOL "Use ISPC code path when evaluating")

if(FTE_ISPC_EVAL_ENABLED)
  add_definitions(-DUSE_ISPC)
endif()

SET(FTE_INIT_COMPONENTS OFF CACHE BOOL "Initialize component values to [0.0,1.0)")

if(FTE_INIT_COMPONENTS)
  add_definitions(-DINIT_COMPONENTS)
endif()

SET(FTE_EVAL_CHUNK_SIZE 4096 CACHE INT "Maximum size of components buffer in-flight (in MB)")

add_definitions(-DCHUNK_SIZE=${FTE_EVAL_CHUNK_SIZE})

include_directories(
  ${CMAKE_SOURCE_DIR}
  )

set(fte_HDRS
    ftParser.h
    ftEvaluator.h
    ftTypes.h
)

set(fte_SRCS
    ftEvaluator.cpp
    ftEvaluator.ispc
    ftParser.cpp
    ftNode.cpp
    ftNodeVisitor.cpp
    ftFileReader.cpp
    ftNamePool.cpp
)

add_ispc_executable(fte ${fte_HDRS} ${fte_SRCS} main.cpp)

###############################################################################
# Add 'make install' target

install(TARGETS ${EXEC_NAME} RUNTIME DESTINATION bin)

