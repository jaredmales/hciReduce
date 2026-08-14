#################################################################
#                    hciReduce test configuration
#################################################################

set(HCIREDUCE_TEST_MAIN_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/testMain.cpp)

add_library(hcireduceTestSupport STATIC common/HCIobservation_test_fixture.cpp)
target_include_directories(hcireduceTestSupport PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hcireduceTestSupport PUBLIC hcireduce-shared)

# Keep one production test source per executable. This bounds compile time and
# makes a failing subsystem immediately visible in CTest output.
set(HCIREDUCE_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/apps/klipReduce_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/ADIobservation_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCI_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_config_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_combine_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_coadd_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_logistics_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_output_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_preprocess_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_preprocess_physics_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_read_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/HCIobservation_RDI_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/common/KLIPreduction_test.cpp
)

set(HCIREDUCE_TEST_TARGETS)

function(hcireduce_add_test_executable TEST_SOURCE)
    file(RELATIVE_PATH _hcireduce_test_rel "${CMAKE_CURRENT_SOURCE_DIR}" "${TEST_SOURCE}")
    string(REPLACE "/" "_" _hcireduce_test_name "${_hcireduce_test_rel}")
    string(REPLACE "." "_" _hcireduce_test_name "${_hcireduce_test_name}")
    set(_hcireduce_test_target "hcireduceTest_${_hcireduce_test_name}")

    add_executable(${_hcireduce_test_target} ${HCIREDUCE_TEST_MAIN_SOURCE} ${TEST_SOURCE})
    target_include_directories(${_hcireduce_test_target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(${_hcireduce_test_target} PRIVATE hcireduceTestSupport)

    if(NOT HCIREDUCE_BUILD_TESTS)
        set_target_properties(
            ${_hcireduce_test_target}
            PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE
        )
    endif()

    if(HCIREDUCE_BUILD_TESTS AND BUILD_TESTING)
        add_test(NAME ${_hcireduce_test_target} COMMAND ${_hcireduce_test_target})
    endif()

    set(HCIREDUCE_TEST_TARGETS ${HCIREDUCE_TEST_TARGETS} ${_hcireduce_test_target} PARENT_SCOPE)
endfunction()

foreach(_hcireduce_test_source IN LISTS HCIREDUCE_TEST_SOURCES)
    hcireduce_add_test_executable("${_hcireduce_test_source}")
endforeach()

add_custom_target(hcireduceTests DEPENDS hcireduce-shared hcireduceTestSupport ${HCIREDUCE_TEST_TARGETS})

if(HCIREDUCE_BUILD_TESTS AND BUILD_TESTING)
    add_custom_target(hcireduceTestRun
        COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR} --output-on-failure
        DEPENDS hcireduceTests
        USES_TERMINAL
        VERBATIM
    )
endif()

set(HCIREDUCE_ONE_TEST
    ""
    CACHE STRING
          "Single test source to build, relative to tests/ (e.g. common/HCIobservation_test.cpp)"
)

if(HCIREDUCE_ONE_TEST)
    if(IS_ABSOLUTE "${HCIREDUCE_ONE_TEST}")
        set(HCIREDUCE_ONE_TEST_SOURCE "${HCIREDUCE_ONE_TEST}")
    else()
        set(HCIREDUCE_ONE_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${HCIREDUCE_ONE_TEST}")
    endif()

    if(HCIREDUCE_ONE_TEST_SOURCE MATCHES "\\.o$")
        string(REGEX REPLACE "\\.o$" ".cpp" HCIREDUCE_ONE_TEST_SOURCE "${HCIREDUCE_ONE_TEST_SOURCE}")
    elseif(HCIREDUCE_ONE_TEST_SOURCE MATCHES "\\.hpp$")
        string(REGEX REPLACE "\\.hpp$" ".cpp" HCIREDUCE_ONE_TEST_SOURCE "${HCIREDUCE_ONE_TEST_SOURCE}")
    elseif(NOT HCIREDUCE_ONE_TEST_SOURCE MATCHES "\\.cpp$")
        string(APPEND HCIREDUCE_ONE_TEST_SOURCE ".cpp")
    endif()

    if(NOT EXISTS "${HCIREDUCE_ONE_TEST_SOURCE}")
        message(FATAL_ERROR "HCIREDUCE_ONE_TEST source does not exist: ${HCIREDUCE_ONE_TEST_SOURCE}")
    endif()

    add_executable(hcireduceTestOne ${HCIREDUCE_TEST_MAIN_SOURCE} ${HCIREDUCE_ONE_TEST_SOURCE})
    target_include_directories(hcireduceTestOne PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(hcireduceTestOne PRIVATE hcireduceTestSupport)

    if(NOT HCIREDUCE_BUILD_TESTS)
        set_target_properties(
            hcireduceTestOne
            PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE
        )
    endif()

    add_custom_target(hcireduceTestOneRun
        COMMAND $<TARGET_FILE:hcireduceTestOne>
        DEPENDS hcireduceTestOne
        USES_TERMINAL
        VERBATIM
    )

    if(HCIREDUCE_BUILD_TESTS AND BUILD_TESTING)
        add_test(NAME hcireduceTestOne COMMAND hcireduceTestOne)
    endif()
endif()
