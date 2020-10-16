#=============================================================================#
#============================== Project ======================================#
#=============================================================================#
cmake_minimum_required(VERSION 3.10)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    project(Htool_benchmarck VERSION 1.0 LANGUAGES CXX)
else()
    project(Htool_benchmarck VERSION 1.0 DESCRIPTION "A header only c++ library that provides Hierarchical matrices." HOMEPAGE_URL "https://github.com/htool-ddm/htool" LANGUAGES CXX)
endif()

# To force c++11
if(${CMAKE_VERSION} VERSION_LESS 3.1)
    add_compile_options(-std=c++11)
elseif(${CMAKE_VERSION} VERSION_LESS 3.6.3 AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    add_compile_options(-std=c++11)
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # To set default CMAKE_BUILD_TYPE
    set(default_build_type "Release")
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
        set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    endif()

    if(${CMAKE_BUILD_TYPE} STREQUAL Release_native)
        message(STATUS "Setting build type to 'Release_native'.")
        set(CMAKE_BUILD_TYPE Release)
        set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -march=native")
    endif()

    message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

    # Sanitizers
    include(cmake/sanitizers.cmake)

    # Formatting
    include(cmake/formatting.cmake)
    file(GLOB_RECURSE ALL_CODE_FILES ${PROJECT_SOURCE_DIR}/include/*.[h]pp ${PROJECT_SOURCE_DIR}/include/*.[h] ${PROJECT_SOURCE_DIR}/examples/*.[ch]pp ${PROJECT_SOURCE_DIR}/examples/*.[ch] ${PROJECT_SOURCE_DIR}/tests/*.[ch]pp ${PROJECT_SOURCE_DIR}/tests/*.[ch])

    clang_format(format ${ALL_CODE_FILES})
    file(GLOB_RECURSE CMAKE_FILES ${PROJECT_SOURCE_DIR}/CMakeLists.txt)
    cmake_format(cmake_format ${CMAKE_FILES})

    # Files to do find_package for some module
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake_modules")

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
message(${OpenMP_CXX_LIBRARIES})

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
add_library(htool_benchmark INTERFACE)
target_include_directories(htool_benchmark INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/external/htool/include> $<INSTALL_INTERFACE:external/htool/include> ${MPI_INCLUDE_PATH} ${HPDDM_INCLUDE_DIRS} ${MKL_INC_DIR})
target_link_libraries(htool_benchmark INTERFACE ${MPI_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

#=============================================================================#
#========================== Lists ======================================#
#=============================================================================#

#=== Compressions
set(compressions "SVD")
# list(APPEND compressions "partialACA")
# list(APPEND compressions "sympartialACA")
# list(APPEND compressions "fullACA")

#=== type
set(types "double")

#=== clustering
set(clusterings "GeometricClustering")

#=== Admissible condition
set(admissible_conditions "RjasanowSteinbach")

#=== Kernel
set(kernels "Gravitational")

#=============================================================================#
#========================== Repertories ======================================#
#=============================================================================#
add_subdirectory(src)
