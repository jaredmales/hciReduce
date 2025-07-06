cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0076 NEW)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project("hciReduce")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=c11")

find_package(OpenMP)
if(OPENMP_FOUND)
message("Found OpenMP")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_library(LIBRT rt)
if(LIBRT_FOUND)
message("Explicitly requiring -lrt")
endif()
find_library(LIBM m)
if(LIBM_FOUND)
message("Explicitly requiring -lm")
endif()
find_library(LIBPTHREAD pthread)
if(LIBPTHREAD_FOUND)
message("Explicitly requiring -lpthread")
endif()

#######################################################################
#
#                            Source
#
#######################################################################
add_subdirectory(src)

#######################################################################
#
#                            Tests
#
#######################################################################
