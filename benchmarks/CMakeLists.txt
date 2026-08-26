add_executable(p4DecompositionBenchmark p4DecompositionBenchmark.cpp)

target_include_directories(p4DecompositionBenchmark PRIVATE "${PROJECT_SOURCE_DIR}")

target_link_libraries(
    p4DecompositionBenchmark
    PRIVATE PkgConfig::MXLIB CUDA::cudart CUDA::cusolver
)
