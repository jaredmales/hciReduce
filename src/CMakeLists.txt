
# this is the "klipReduce" target
add_executable(klipReduce)

target_sources(klipReduce PRIVATE "apps/klipReduce.cpp")

target_include_directories(klipReduce PUBLIC "${PROJECT_SOURCE_DIR}" "/usr/include/eigen3")
add_subdirectory("common")

target_link_libraries(klipReduce mxlib cfitsio)

install(TARGETS klipReduce)


