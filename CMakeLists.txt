cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0076 NEW)

set(CMAKE_VERBOSE_MAKEFILE ON)

include(FindPkgConfig)
include(CheckLanguage)
find_package(PkgConfig REQUIRED)

project("hciReduce")

pkg_check_modules(MXLIB REQUIRED IMPORTED_TARGET mxlib)
string(REPLACE ";" " " MXLIB_PC_CFLAGS "${MXLIB_CFLAGS}")
set(CMAKE_C_FLAGS "${MXLIB_PC_CFLAGS}")


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
