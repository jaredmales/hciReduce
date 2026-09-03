if(HCIREDUCE_BUILD_CPU_BENCHMARKS)
    add_executable(klipPrecisionBenchmark klipPrecisionBenchmark.cpp)

    target_include_directories(klipPrecisionBenchmark PRIVATE "${PROJECT_SOURCE_DIR}")

    target_link_libraries(klipPrecisionBenchmark PRIVATE PkgConfig::MXLIB)

    add_executable(p4PrecisionBenchmark p4PrecisionBenchmark.cpp)

    target_link_libraries(p4PrecisionBenchmark PRIVATE PkgConfig::MXLIB)

    if(HCIREDUCE_ENABLE_EXPERIMENTAL_P4_PRECISION)
        add_executable(
            p4ReductionPrecisionBenchmark
            p4ReductionPrecisionBenchmark.cpp
        )

        target_include_directories(
            p4ReductionPrecisionBenchmark PRIVATE "${PROJECT_SOURCE_DIR}"
        )

        target_link_libraries(p4ReductionPrecisionBenchmark PRIVATE hcireduce-shared)
    endif()
endif()

if(HCIREDUCE_BUILD_CUDA_BENCHMARKS)
    add_executable(p4DecompositionBenchmark p4DecompositionBenchmark.cpp)

    target_include_directories(p4DecompositionBenchmark PRIVATE "${PROJECT_SOURCE_DIR}")

    target_link_libraries(
        p4DecompositionBenchmark
        PRIVATE PkgConfig::MXLIB CUDA::cudart CUDA::cusolver
    )
endif()
