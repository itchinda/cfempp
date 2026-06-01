cmake_minimum_required(VERSION 3.20)

project(cfem LANGUAGES C CXX)

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler CXX: ${CMAKE_CXX_COMPILER}")

# ==============================================================================
# IPO / LTO support detection
# ==============================================================================
include(CheckIPOSupported)

check_ipo_supported(RESULT ipo_supported)

if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    message(STATUS "IPO/LTO enabled for Release builds")
else()
    message(STATUS "IPO/LTO not supported by this toolchain")
endif()

# ==============================================================================
# HPC Optimization Target (Modern CMake)
# ==============================================================================
add_library(cfem_hpc_flags INTERFACE)

# Release flags
target_compile_options(cfem_hpc_flags INTERFACE 
    $<$<CONFIG:Release>:-O3 -flto -fopenmp-simd -march=native -ffp-contract=fast -fno-math-errno>
    $<$<CONFIG:Debug>:-O0 -g -fno-omit-frame-pointer>
)

# Debug macros management (PETSc, Teuchos, NDEBUG)
target_compile_definitions(cfem_hpc_flags INTERFACE
    $<$<CONFIG:Release>:NDEBUG TEUCHOS_DEBUG_OFF>
    $<$<CONFIG:Debug>:PETSC_USE_DEBUG>
)

# ==============================================================================
# Dependencies
# ==============================================================================

# PETSc
if(NOT DEFINED PETSC_DIR OR NOT DEFINED PETSC_ARCH)
    message(FATAL_ERROR "PETSC_DIR and PETSC_ARCH must be defined via CMakePresets.")
endif()

set(PETSC_PATH "${PETSC_DIR}/${PETSC_ARCH}")
set(ENV{PKG_CONFIG_PATH} "${PETSC_PATH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET GLOBAL petsc)

# MPI
find_library(PETSC_MPI_LIB NAMES mpi HINTS "${PETSC_PATH}/lib" NO_DEFAULT_PATH)
target_link_libraries(PkgConfig::PETSC INTERFACE ${PETSC_MPI_LIB})

# Math Libraries (GMP, MPFR, MPC)
find_library(GMP_LIB NAMES gmp REQUIRED)
find_library(MPFR_LIB NAMES mpfr REQUIRED)
find_library(MPC_LIB NAMES mpc REQUIRED)

# SymEngine
find_package(SymEngine REQUIRED CONFIG)
if(SymEngine_FOUND)
    message(STATUS "Found SymEngine version: ${SymEngine_VERSION}")
    
    if(TARGET symengine)
        target_link_libraries(symengine INTERFACE ${GMP_LIB} ${MPFR_LIB} ${MPC_LIB})
        add_library(SymEngine::libsymengine ALIAS symengine)
    elseif(TARGET libsymengine)
        target_link_libraries(libsymengine INTERFACE ${GMP_LIB} ${MPFR_LIB} ${MPC_LIB})
        add_library(SymEngine::libsymengine ALIAS libsymengine)
    endif()
endif()

# LLVM
find_package(LLVM REQUIRED CONFIG)
if(LLVM_FOUND)
    llvm_map_components_to_libnames(llvm_libs core mcjit orcjit native ipo bitreader target mc support)
    message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}")
endif()

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    add_library(cfem_openmp_interface INTERFACE)
    target_link_libraries(cfem_openmp_interface INTERFACE OpenMP::OpenMP_CXX)
    target_compile_definitions(cfem_openmp_interface INTERFACE WITH_OPENMP)
endif()

# ==============================================================================
# Subdirectories
# ==============================================================================
add_subdirectory(libs)
add_subdirectory(apps)

# ==============================================================================
# Testing (GoogleTest)
# ==============================================================================
include(CTest)
if(BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.14.0
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_subdirectory(tests)
endif()

# ==============================================================================
# Git Hooks Configuration
# ==============================================================================
if(EXISTS "${CMAKE_SOURCE_DIR}/.githooks")
    execute_process(
        COMMAND git config core.hooksPath .githooks
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        RESULT_VARIABLE GIT_HOOKS_RESULT
        OUTPUT_QUIET ERROR_QUIET
    )
    if(GIT_HOOKS_RESULT EQUAL 0)
        message(STATUS "Git hooks correctly configured to use .githooks/")
    endif()
endif()