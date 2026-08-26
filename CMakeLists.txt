cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0076 NEW)

set(CMAKE_VERBOSE_MAKEFILE ON)

include(FindPkgConfig)
find_package(PkgConfig REQUIRED)

project("hciReduce" LANGUAGES CXX)

pkg_check_modules(MXLIB REQUIRED IMPORTED_TARGET mxlib)

option(HCIREDUCE_BUILD_TESTS "Whether or not to build hciReduce tests" OFF)
option(HCIREDUCE_BUILD_COVERAGE "Whether or not to enable coverage instrumentation" OFF)
set(HCIREDUCE_COVERAGE_TEST_TIMEOUT 300 CACHE STRING "CTest timeout (seconds) for the coverage target")
option(HCIREDUCE_BUILD_DOCS "Whether or not to build Doxygen docs as part of the default build" OFF)
option(HCIREDUCE_BUILD_CUDA_BENCHMARKS "Whether or not to build standalone CUDA performance benchmarks" OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    if(HCIREDUCE_BUILD_COVERAGE)
        set(_hcireduce_default_build_type Debug)
    else()
        set(_hcireduce_default_build_type Release)
    endif()
    set(CMAKE_BUILD_TYPE ${_hcireduce_default_build_type} CACHE STRING "Choose the build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
    message(STATUS "No build type selected; defaulting to ${_hcireduce_default_build_type}")
endif()

if(HCIREDUCE_BUILD_COVERAGE)
    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(FATAL_ERROR "HCIREDUCE_BUILD_COVERAGE currently requires GCC so lcov can read gcov data")
    endif()

    # Link this after PkgConfig::MXLIB so coverage options remain last.
    add_library(hcireduce-coverage-options INTERFACE)
    target_compile_definitions(hcireduce-coverage-options INTERFACE HCIREDUCE_BUILD_COVERAGE=1)
    target_compile_options(hcireduce-coverage-options INTERFACE --coverage -O0 -fno-fast-math)
    target_link_options(hcireduce-coverage-options INTERFACE --coverage)
endif()


#######################################################################
#
#                            Source
#
#######################################################################
add_subdirectory(src)

if(HCIREDUCE_BUILD_CUDA_BENCHMARKS)
    find_package(CUDAToolkit 13.0 REQUIRED)
    add_subdirectory(benchmarks)
endif()

#######################################################################
#
#                            Tests
#
#######################################################################

if(HCIREDUCE_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
else()
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()

add_custom_target(tests DEPENDS hcireduceTests)

#######################################################################
#
#                    Documentation and Coverage
#
#######################################################################

set(HCIREDUCE_COVERAGE_BUILD_DIR ${CMAKE_BINARY_DIR}/coverage-build)
set(HCIREDUCE_COVERAGE_INFO ${CMAKE_BINARY_DIR}/coverage.info)
set(HCIREDUCE_COVERAGE_FILTERED_INFO ${CMAKE_BINARY_DIR}/coverage_filtered.info)
set(HCIREDUCE_DOC_OUTPUT_DIR ${CMAKE_BINARY_DIR}/doc)
set(HCIREDUCE_DOC_HTML_DIR ${HCIREDUCE_DOC_OUTPUT_DIR}/html)
set(HCIREDUCE_COVERAGE_REPORT_DIR ${HCIREDUCE_DOC_HTML_DIR}/coverage)

find_package(Doxygen QUIET)

if(Doxygen_FOUND)
    set(HCIREDUCE_DOXYGEN_CONFIG_IN ${CMAKE_SOURCE_DIR}/doc/hcireduce_doxygen.in)
    set(HCIREDUCE_DOXYGEN_CONFIG ${CMAKE_BINARY_DIR}/hcireduce_doxygen)
    set(HCIREDUCE_DOXYGEN_OUTPUT_DIRECTORY ${HCIREDUCE_DOC_OUTPUT_DIR})
    set(HCIREDUCE_DOC_FIGURE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/doc/figures)
    set(HCIREDUCE_DOC_FIGURE_BUILD_DIR ${CMAKE_BINARY_DIR}/doc-figures)

    # Every standalone TikZ source produces same-named PDF and SVG assets. CONFIGURE_DEPENDS
    # makes newly added diagrams join the pipeline after CMake automatically regenerates.
    file(GLOB HCIREDUCE_DOC_TIKZ_SOURCES CONFIGURE_DEPENDS ${HCIREDUCE_DOC_FIGURE_SOURCE_DIR}/*.tex)
    find_program(HCIREDUCE_LATEXMK_EXECUTABLE latexmk)
    find_program(HCIREDUCE_PDFTOCAIRO_EXECUTABLE pdftocairo)
    set(HCIREDUCE_DOC_FIGURE_OUTPUTS)
    set(_hcireduce_doc_figure_extra_files)

    if(HCIREDUCE_DOC_TIKZ_SOURCES AND HCIREDUCE_LATEXMK_EXECUTABLE AND HCIREDUCE_PDFTOCAIRO_EXECUTABLE)
        set(HCIREDUCE_DOC_FIGURE_IMAGE_PATH ${HCIREDUCE_DOC_FIGURE_BUILD_DIR})
        foreach(_hcireduce_tikz_source IN LISTS HCIREDUCE_DOC_TIKZ_SOURCES)
            get_filename_component(_hcireduce_figure_stem ${_hcireduce_tikz_source} NAME_WE)
            set(_hcireduce_figure_pdf ${HCIREDUCE_DOC_FIGURE_BUILD_DIR}/${_hcireduce_figure_stem}.pdf)
            set(_hcireduce_figure_svg ${HCIREDUCE_DOC_FIGURE_BUILD_DIR}/${_hcireduce_figure_stem}.svg)
            add_custom_command(
                OUTPUT ${_hcireduce_figure_pdf} ${_hcireduce_figure_svg}
                COMMAND ${CMAKE_COMMAND} -E make_directory ${HCIREDUCE_DOC_FIGURE_BUILD_DIR}
                COMMAND ${HCIREDUCE_LATEXMK_EXECUTABLE} -pdf -interaction=nonstopmode -halt-on-error
                        -outdir=${HCIREDUCE_DOC_FIGURE_BUILD_DIR} ${_hcireduce_tikz_source}
                COMMAND ${HCIREDUCE_PDFTOCAIRO_EXECUTABLE} -svg ${_hcireduce_figure_pdf} ${_hcireduce_figure_svg}
                COMMAND ${HCIREDUCE_LATEXMK_EXECUTABLE} -c -outdir=${HCIREDUCE_DOC_FIGURE_BUILD_DIR}
                        ${_hcireduce_tikz_source}
                DEPENDS ${_hcireduce_tikz_source}
                WORKING_DIRECTORY ${HCIREDUCE_DOC_FIGURE_SOURCE_DIR}
                COMMENT "Rendering documentation diagram ${_hcireduce_figure_stem}"
                VERBATIM
            )
            list(APPEND HCIREDUCE_DOC_FIGURE_OUTPUTS ${_hcireduce_figure_pdf} ${_hcireduce_figure_svg})
            list(APPEND _hcireduce_doc_figure_extra_files ${_hcireduce_figure_svg} ${_hcireduce_figure_pdf})
        endforeach()
    else()
        set(HCIREDUCE_DOC_FIGURE_IMAGE_PATH ${HCIREDUCE_DOC_FIGURE_SOURCE_DIR})
        foreach(_hcireduce_tikz_source IN LISTS HCIREDUCE_DOC_TIKZ_SOURCES)
            get_filename_component(_hcireduce_figure_stem ${_hcireduce_tikz_source} NAME_WE)
            list(APPEND _hcireduce_doc_figure_extra_files
                 ${HCIREDUCE_DOC_FIGURE_SOURCE_DIR}/${_hcireduce_figure_stem}.svg
                 ${HCIREDUCE_DOC_FIGURE_SOURCE_DIR}/${_hcireduce_figure_stem}.pdf)
        endforeach()
        if(HCIREDUCE_DOC_TIKZ_SOURCES)
            message(WARNING "latexmk or pdftocairo was not found; documentation will use pre-rendered diagram assets")
        endif()
    endif()

    string(JOIN " \\\n                          " HCIREDUCE_DOC_FIGURE_HTML_EXTRA_FILES
           ${_hcireduce_doc_figure_extra_files})

    configure_file(${HCIREDUCE_DOXYGEN_CONFIG_IN} ${HCIREDUCE_DOXYGEN_CONFIG} @ONLY)

    add_custom_target(doc_figures DEPENDS ${HCIREDUCE_DOC_FIGURE_OUTPUTS})

    if(HCIREDUCE_BUILD_DOCS)
        set(_hcireduce_docs_all ALL)
    else()
        set(_hcireduce_docs_all)
    endif()

    add_custom_target(docs ${_hcireduce_docs_all}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${HCIREDUCE_DOC_OUTPUT_DIR}
        COMMAND ${DOXYGEN_EXECUTABLE} ${HCIREDUCE_DOXYGEN_CONFIG}
        COMMAND ${CMAKE_COMMAND} -DHCIREDUCE_DOC_HTML_DIR=${HCIREDUCE_DOC_HTML_DIR} -P ${CMAKE_SOURCE_DIR}/cmake/restructure_doxygen_navtree.cmake
        COMMAND ${CMAKE_COMMAND} -DCOVERAGE_REPORT_DIR=${HCIREDUCE_COVERAGE_REPORT_DIR} -DPROJECT_NAME=hciReduce -P ${CMAKE_SOURCE_DIR}/cmake/ensure_coverage_placeholder.cmake
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
        DEPENDS ${HCIREDUCE_DOXYGEN_CONFIG}
                ${CMAKE_SOURCE_DIR}/doc/hcireduce_header.html
                ${CMAKE_SOURCE_DIR}/doc/hcireduce-custom.css
                ${CMAKE_SOURCE_DIR}/doc/hcireduceDoxygenLayout.xml
                ${CMAKE_SOURCE_DIR}/doc/hcireduce.bib
                ${CMAKE_SOURCE_DIR}/doc/main.dox
                ${CMAKE_SOURCE_DIR}/doc/groupdefs.dox
                ${CMAKE_SOURCE_DIR}/doc/introduction.dox
                ${CMAKE_SOURCE_DIR}/doc/common.dox
                ${CMAKE_SOURCE_DIR}/doc/common_config.dox
                ${CMAKE_SOURCE_DIR}/doc/klip.dox
                ${CMAKE_SOURCE_DIR}/doc/klip_config.dox
                ${CMAKE_SOURCE_DIR}/doc/p4.dox
                ${CMAKE_SOURCE_DIR}/doc/p4_config.dox
                ${CMAKE_SOURCE_DIR}/doc/hciAnalyze.dox
                ${CMAKE_SOURCE_DIR}/doc/hciAnalyze_config.dox
                ${CMAKE_SOURCE_DIR}/doc/programming.dox
                ${CMAKE_SOURCE_DIR}/doc/klip_algorithm.dox
                ${CMAKE_SOURCE_DIR}/doc/p4_algorithm.dox
                ${CMAKE_SOURCE_DIR}/doc/testing.dox
                ${CMAKE_SOURCE_DIR}/doc/building_tests.dox
                ${CMAKE_SOURCE_DIR}/doc/coverage_report.dox
                ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css/doxygen-awesome.css
                ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css/doxygen-awesome-sidebar-only.css
                ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css/doxygen-awesome-sidebar-only-darkmode-toggle.css
                ${CMAKE_SOURCE_DIR}/doc/hcireduce-custom.css
                ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css/doxygen-awesome-darkmode-toggle.js
                ${CMAKE_SOURCE_DIR}/doc/doxygen-awesome-css/doxygen-awesome-fragment-copy-button.js
                ${HCIREDUCE_DOC_FIGURE_OUTPUTS}
                ${CMAKE_SOURCE_DIR}/cmake/ensure_coverage_placeholder.cmake
                ${CMAKE_SOURCE_DIR}/cmake/restructure_doxygen_navtree.cmake
        USES_TERMINAL
        VERBATIM
    )
    add_dependencies(docs doc_figures)
else()
    message(WARNING "doxygen not found. The docs target will not be available.")
endif()

find_program(HCIREDUCE_LCOV_EXECUTABLE lcov)
find_program(HCIREDUCE_GENHTML_EXECUTABLE genhtml)

if(HCIREDUCE_LCOV_EXECUTABLE AND HCIREDUCE_GENHTML_EXECUTABLE)
    add_custom_target(coverage
        COMMAND ${CMAKE_COMMAND} -S ${CMAKE_SOURCE_DIR} -B ${HCIREDUCE_COVERAGE_BUILD_DIR} -G "${CMAKE_GENERATOR}" -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_BUILD_TYPE=Debug -DHCIREDUCE_BUILD_TESTS=ON -DHCIREDUCE_BUILD_COVERAGE=ON
        COMMAND ${CMAKE_COMMAND} --build ${HCIREDUCE_COVERAGE_BUILD_DIR} --target hcireduceTests -j
        COMMAND ${HCIREDUCE_LCOV_EXECUTABLE} --directory ${HCIREDUCE_COVERAGE_BUILD_DIR} --zerocounters
        COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${HCIREDUCE_COVERAGE_BUILD_DIR} --output-on-failure --timeout ${HCIREDUCE_COVERAGE_TEST_TIMEOUT}
        COMMAND ${HCIREDUCE_LCOV_EXECUTABLE} --directory ${HCIREDUCE_COVERAGE_BUILD_DIR} --capture --output-file ${HCIREDUCE_COVERAGE_INFO}
        COMMAND ${HCIREDUCE_LCOV_EXECUTABLE} --extract ${HCIREDUCE_COVERAGE_INFO} ${CMAKE_SOURCE_DIR}/src/* --ignore-errors unused,unused --ignore-errors inconsistent,inconsistent --output-file ${HCIREDUCE_COVERAGE_FILTERED_INFO}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${HCIREDUCE_COVERAGE_REPORT_DIR}
        COMMAND ${HCIREDUCE_GENHTML_EXECUTABLE} ${HCIREDUCE_COVERAGE_FILTERED_INFO} --output-directory ${HCIREDUCE_COVERAGE_REPORT_DIR} --title hciReduce --hierarchical --merge-aliases --suppress-aliases --filter function --demangle-cpp --css-file ${CMAKE_SOURCE_DIR}/tests/coverage/gcov.css
        USES_TERMINAL
        VERBATIM
    )
    if(TARGET docs)
        add_dependencies(coverage docs)
    endif()
else()
    message(WARNING "lcov/genhtml not found. The coverage target will not be available.")
endif()

add_custom_target(coverage_clean
    COMMAND ${CMAKE_COMMAND} -E rm -f ${HCIREDUCE_COVERAGE_INFO}
    COMMAND ${CMAKE_COMMAND} -E rm -f ${HCIREDUCE_COVERAGE_FILTERED_INFO}
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${HCIREDUCE_COVERAGE_REPORT_DIR}
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${HCIREDUCE_COVERAGE_BUILD_DIR}
)

if(TARGET docs)
    add_custom_target(docs_clean
        COMMAND ${CMAKE_COMMAND} -E rm -rf ${HCIREDUCE_DOC_OUTPUT_DIR}
        COMMAND ${CMAKE_COMMAND} -E rm -rf ${HCIREDUCE_DOC_FIGURE_BUILD_DIR}
    )
endif()
