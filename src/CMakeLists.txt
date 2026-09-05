


############################
# libhcireduce
############################

# list of source files
set(libheaders common/HCI.hpp
               common/ConfigUtils.hpp
               common/HCIobservation.hpp
               common/ADIobservation.hpp
               common/ADIDerotator.hpp
               common/P4PixelGrid.hpp
               common/P4RotatedGrid.hpp
               common/P4PCA.hpp
               common/P4TemporalPCA.hpp
               common/P4LocalProcessing.hpp
               common/P4NegativeOptimizer.hpp
               common/KLIPPSFModel.hpp
               common/RadialPSFModel.hpp
               common/P4PSFFilter.hpp
               common/P4PSFModel.hpp
               common/P4PSFReconstructor.hpp
               common/P4Reduction.hpp
               common/KLIPreduction.hpp )

set(libsource common/HCIobservation.cpp
              common/ADIobservation.cpp
              common/ADIDerotator.cpp
              common/P4PixelGrid.cpp
              common/P4RotatedGrid.cpp
              common/P4PCA.cpp
              common/P4TemporalPCA.cpp
              common/P4LocalProcessing.cpp
              common/P4NegativeOptimizer.cpp
              common/KLIPPSFModel.cpp
              common/RadialPSFModel.cpp
              common/P4PSFFilter.cpp
              common/P4PSFModel.cpp
              common/P4PSFReconstructor.cpp
              common/P4Reduction.cpp
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

# this is the "hciAnalyze" target
add_executable(hciAnalyze)
target_sources(hciAnalyze PRIVATE "apps/hciAnalyze.cpp")

target_include_directories(hciAnalyze PUBLIC "${PROJECT_SOURCE_DIR}")

# this is the "p4Reduce" target
add_executable(p4Reduce)
target_sources(p4Reduce PRIVATE "apps/p4Reduce.cpp")

target_include_directories(p4Reduce PUBLIC "${PROJECT_SOURCE_DIR}")

target_link_libraries(objlib PUBLIC PkgConfig::MXLIB)
target_link_libraries(hcireduce-shared PUBLIC PkgConfig::MXLIB)
target_link_libraries(hcireduce-static PUBLIC PkgConfig::MXLIB)

target_link_libraries(klipReduce PUBLIC hcireduce-shared)
target_link_libraries(hciAnalyze PUBLIC hcireduce-shared)
target_link_libraries(p4Reduce PUBLIC hcireduce-shared)

if(HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION)
    target_compile_definitions(objlib PRIVATE HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION=1)
    target_compile_definitions(hcireduce-shared INTERFACE HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION=1)
    target_compile_definitions(hcireduce-static INTERFACE HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION=1)
endif()

# Keep installed executables paired with the libhcireduce from the same prefix.
set_target_properties(klipReduce hciAnalyze p4Reduce PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")

if(TARGET hcireduce-coverage-options)
    target_link_libraries(objlib PUBLIC hcireduce-coverage-options)
    target_link_libraries(hcireduce-shared PUBLIC hcireduce-coverage-options)
    target_link_libraries(hcireduce-static PUBLIC hcireduce-coverage-options)
    target_link_libraries(klipReduce PUBLIC hcireduce-coverage-options)
    target_link_libraries(hciAnalyze PUBLIC hcireduce-coverage-options)
    target_link_libraries(p4Reduce PUBLIC hcireduce-coverage-options)
endif()

install(TARGETS hcireduce-shared DESTINATION lib)
install(TARGETS hcireduce-static DESTINATION lib)
install(FILES ${libheaders} DESTINATION include/hcireduce)
install(TARGETS klipReduce hciAnalyze p4Reduce)
