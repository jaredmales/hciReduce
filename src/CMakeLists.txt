


############################
# libhcireduce
############################

# list of source files
set(libheaders common/HCI.hpp
               common/HCIobservation.hpp
               common/ADIobservation.hpp
               common/ADIDerotator.hpp
               common/KLIPreduction.hpp )

set(libsource common/HCIobservation.cpp
              common/ADIobservation.cpp
              common/ADIDerotator.cpp
              common/KLIPreduction.cpp )

# this is the "object library" target: compiles the sources only once
add_library(objlib OBJECT ${libheaders} ${libsource})

# shared libraries need PIC
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

target_include_directories(objlib PUBLIC "${PROJECT_SOURCE_DIR}")

# shared and static libraries built from the same object files
add_library(hcireduce-shared SHARED $<TARGET_OBJECTS:objlib>)
add_library(hcireduce-static STATIC $<TARGET_OBJECTS:objlib>)

target_include_directories(hcireduce-shared PUBLIC "${PROJECT_SOURCE_DIR}")
target_include_directories(hcireduce-static PUBLIC "${PROJECT_SOURCE_DIR}")

#This makes it so we have libhcireduce.so and libhcireduce.a, without the -shared/-static
SET_TARGET_PROPERTIES(hcireduce-shared PROPERTIES OUTPUT_NAME hcireduce CLEAN_DIRECT_OUTPUT 1)
SET_TARGET_PROPERTIES(hcireduce-static PROPERTIES OUTPUT_NAME hcireduce CLEAN_DIRECT_OUTPUT 1)


# this is the "klipReduce" target
add_executable(klipReduce)
target_sources(klipReduce PRIVATE "apps/klipReduce.cpp")

target_include_directories(klipReduce PUBLIC "${PROJECT_SOURCE_DIR}")

target_link_libraries(objlib PUBLIC PkgConfig::MXLIB)
target_link_libraries(hcireduce-shared PUBLIC PkgConfig::MXLIB)
target_link_libraries(hcireduce-static PUBLIC PkgConfig::MXLIB)

target_link_libraries(klipReduce PUBLIC hcireduce-shared)

if(TARGET hcireduce-coverage-options)
    target_link_libraries(objlib PUBLIC hcireduce-coverage-options)
    target_link_libraries(hcireduce-shared PUBLIC hcireduce-coverage-options)
    target_link_libraries(hcireduce-static PUBLIC hcireduce-coverage-options)
    target_link_libraries(klipReduce PUBLIC hcireduce-coverage-options)
endif()

install(TARGETS hcireduce-shared DESTINATION lib)
install(TARGETS hcireduce-static DESTINATION lib)
install(FILES ${libheaders} DESTINATION include/hcireduce)
install(TARGETS klipReduce)
