#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.10)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    project(Htool_Benchmarks VERSION 1.0 LANGUAGES CXX)
else()
    project(Htool_Benchmarks VERSION 1.0 DESCRIPTION "Benchmarks for Htool" HOMEPAGE_URL "https://github.com/htool-ddm/htool" LANGUAGES CXX)
endif()

# To force c++11
if(${CMAKE_VERSION} VERSION_LESS 3.1)
    add_compile_options(-std=c++14)
elseif(${CMAKE_VERSION} VERSION_LESS 3.6.3 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    add_compile_options(-std=c++14)
else()
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # To set default CMAKE_BUILD_TYPE
    set(default_build_type "Release")
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
        set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    endif()

    if(${CMAKE_BUILD_TYPE} STREQUAL Release_w_vecto)
        message(STATUS "Setting build type to 'Release_w_vecto'.")
        set(CMAKE_BUILD_TYPE Release)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math")
    endif()

    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

    # Information about compilation exported
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.
    include(CTest)

    # # Docs only available if this is the main app
    # find_package(Doxygen)
    # if(Doxygen_FOUND)
    # add_subdirectory(docs)
    # else()
    # message(STATUS "Doxygen not found, not building docs")
    # endif()

    # Sanitizers
    include(cmake/sanitizers.cmake)

endif()

#=============================================================================#
#========================== External Libraries ===============================#
#=============================================================================#
# MPI
find_package(MPI REQUIRED)
message("-- MPI libraries found in " "${MPI_LIBRARIES}")
message("-- MPI include files found in " "${MPI_INCLUDE_PATH}")
separate_arguments(MPIEXEC_PREFLAGS) # to support multi flags

# OPENMP
find_package(OpenMP)

# BLAS
find_package(BLAS REQUIRED)
message("-- Found Blas implementation:" "${BLAS_LIBRARIES}")

# LAPACK
find_package(LAPACK)
message("-- Found Lapack:" "${LAPACK_LIBRARIES}")

#=============================================================================#
#========================== Libraries ========================================#
#=============================================================================#
#=== HTOOL as header only library
add_library(htool INTERFACE)
target_include_directories(htool INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/htool/include> $<INSTALL_INTERFACE:external/htool/include> ${MPI_INCLUDE_PATH} ${HPDDM_INCLUDE_DIRS} ${MKL_INC_DIR})
target_link_libraries(htool INTERFACE ${MPI_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} OpenMP::OpenMP_CXX)

target_compile_features(htool INTERFACE cxx_std_14)

if("${BLA_VENDOR}" STREQUAL "Intel10_64lp_seq" OR "${BLA_VENDOR}" STREQUAL "Intel10_64lp")
    target_compile_definitions(htool INTERFACE "-DHTOOL_MKL")
endif()

#=== benchmark
# set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)

# add_subdirectory(external/benchmark)
# add_library(benchmark INTERFACE)
# target_include_directories(
#     benchmark
#     INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/benchmark/include>
#               $<INSTALL_INTERFACE:external/benchmark/include>)

add_library(htool_benchmarks INTERFACE)
target_link_libraries(htool_benchmarks INTERFACE htool)
target_include_directories(htool_benchmarks INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/> $<INSTALL_INTERFACE:include/>)

#=============================================================================#
#========================== Benchmarks =======================================#
#=============================================================================#

add_custom_target(build-benchmarks)

# add_executable(bench_hmatrix_build_vs_pbl_size src/bench_hmatrix_build_vs_pbl_size.cpp)
# target_link_libraries(bench_hmatrix_build_vs_pbl_size htool benchmark::benchmark benchmark::benchmark_main)

# add_executable(bench_hmatrix_build_vs_thread_numbers src/bench_hmatrix_build_vs_thread_numbers.cpp)
# target_link_libraries(bench_hmatrix_build_vs_thread_numbers htool benchmark::benchmark benchmark::benchmark_main)

# add_executable(bench_hmatrix_matrix_product_vs_pbl_size src/bench_hmatrix_matrix_product_vs_pbl_size.cpp)
# target_link_libraries(bench_hmatrix_matrix_product_vs_pbl_size htool benchmark::benchmark benchmark::benchmark_main)

# add_executable(bench_hmatrix_matrix_product_vs_thread_numbers src/bench_hmatrix_matrix_product_vs_thread_numbers.cpp)
# target_link_libraries(bench_hmatrix_matrix_product_vs_thread_numbers htool benchmark::benchmark benchmark::benchmark_main)

# add_executable(bench_hmatrix_matrix_product src/bench_hmatrix_matrix_product.cpp)
# target_link_libraries(bench_hmatrix_matrix_product htool benchmark::benchmark benchmark::benchmark_main)

# add_executable(bench_hmatrix_factorization src/bench_hmatrix_factorization.cpp)
# target_link_libraries(bench_hmatrix_factorization htool benchmark::benchmark benchmark::benchmark_main)

# add_dependencies(build-benchmarks bench_hmatrix_build_vs_pbl_size bench_hmatrix_build_vs_thread_numbers bench_hmatrix_matrix_product_vs_pbl_size bench_hmatrix_matrix_product_vs_thread_numbers bench_hmatrix_factorization bench_hmatrix_matrix_product)

#=============================================================================#
#=====================VS pbl_size benchmarks==================================#
#=============================================================================#
add_executable(bench_hmatrix_build_vs_pbl_size src/Bench_WO_GBm/bench_hmatrix_build_vs_pbl_size.cpp)
target_link_libraries(bench_hmatrix_build_vs_pbl_size htool_benchmarks)

add_executable(bench_hmatrix_factorization_vs_pbl_size src/Bench_WO_GBm/bench_hmatrix_factorization_vs_pbl_size.cpp)
target_link_libraries(bench_hmatrix_factorization_vs_pbl_size htool_benchmarks)

add_executable(bench_hmatrix_matrix_product_vs_pbl_size src/Bench_WO_GBm/bench_hmatrix_matrix_product_vs_pbl_size.cpp)
target_link_libraries(bench_hmatrix_matrix_product_vs_pbl_size htool_benchmarks)

#=============================================================================#
#=====================VS thread benchmarks====================================#
#=============================================================================#
add_executable(bench_hmatrix_build_vs_thread src/Bench_WO_GBm/bench_hmatrix_build_vs_thread.cpp)
target_link_libraries(bench_hmatrix_build_vs_thread htool_benchmarks)

add_executable(bench_hmatrix_matrix_product_vs_thread src/Bench_WO_GBm/bench_hmatrix_matrix_product_vs_thread.cpp)
target_link_libraries(bench_hmatrix_matrix_product_vs_thread htool_benchmarks)

#=============================================================================#
#=====================VS ratio benchmarks====================================#
#=============================================================================#
add_executable(bench_hmatrix_build_vs_ratio src/Bench_WO_GBm/bench_hmatrix_build_vs_ratio.cpp)
target_link_libraries(bench_hmatrix_build_vs_ratio htool_benchmarks)

add_executable(bench_hmatrix_matrix_product_vs_ratio src/Bench_WO_GBm/bench_hmatrix_matrix_product_vs_ratio.cpp)
target_link_libraries(bench_hmatrix_matrix_product_vs_ratio htool_benchmarks)

add_dependencies(build-benchmarks bench_hmatrix_build_vs_pbl_size bench_hmatrix_factorization_vs_pbl_size bench_hmatrix_matrix_product_vs_pbl_size bench_hmatrix_build_vs_thread bench_hmatrix_matrix_product_vs_thread bench_hmatrix_build_vs_ratio bench_hmatrix_matrix_product_vs_ratio)
