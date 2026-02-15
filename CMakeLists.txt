# # ------------------------------------------------------------
# # Compilation modes:
# #   cmake .. -DCMAKE_BUILD_TYPE=Release
# #   cmake .. -DCMAKE_BUILD_TYPE=Debug
# # ------------------------------------------------------------

cmake_minimum_required(VERSION 3.15)

# ------------------------------------------------------------
# 1. Build type / PETSc configuration (Must happen before project())
# ------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(TARGET_PETSC_ARCH "arch-linux-c-debug")
    add_definitions(-DPETSC_USE_DEBUG)
else()
    set(TARGET_PETSC_ARCH "arch-linux-c-opt")
endif()

set(ENV{PETSC_ARCH} ${TARGET_PETSC_ARCH})

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Forced PETSc Architecture: $ENV{PETSC_ARCH}")

set(PETSC_PATH $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
set(ENV{PKG_CONFIG_PATH} "${PETSC_PATH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")

find_package(PkgConfig REQUIRED)

pkg_get_variable(PETSC_CXX PETSc cxxcompiler)
if(PETSC_CXX)
    set(CMAKE_CXX_COMPILER ${PETSC_CXX})
endif()

# ------------------------------------------------------------
# Project Initialization
# ------------------------------------------------------------
project(cfem LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# set(CMAKE_CXX_FLAGS_ASAN "-O1 -g -fsanitize=address -fno-omit-frame-pointer" CACHE STRING "") # For intense debug/segmentation fault/ illegal memory acces
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG -DTEUCHOS_DEBUG_OFF" CACHE STRING "")   # For pure performance (HPC)

# Define the App name here so it's available in the subdirectories
if(TARGET_PETSC_ARCH STREQUAL "arch-linux-c-debug")
    set(APP_PROJECT_NAME cfem_app_dev)
    set(SYMENGINE_BUILD_DIR "build_debug")
else()
    set(APP_PROJECT_NAME cfem_app_opt)
    set(SYMENGINE_BUILD_DIR "build_release")
endif()

# Global PETSc target available to all subdirectories
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET petsc)

# ------------------------------------------------------------
# Find SymEngine (Sibling directory lookup)
# ------------------------------------------------------------

set(SymEngine_DIR "${CMAKE_SOURCE_DIR}/../symengine/${SYMENGINE_BUILD_DIR}")
find_package(SymEngine REQUIRED CONFIG)

message(STATUS "Found SymEngine at: ${SymEngine_DIR}")

find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP")
endif()

# ------------------------------------------------------------
# Add Modular Libraries & Applications
# ------------------------------------------------------------
add_subdirectory(libs)
add_subdirectory(apps) # This will look for apps/CMakeLists.txt

# ------------------------------------------------------------
# Testing / GoogleTest
# ------------------------------------------------------------
include(CTest)

if(BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.17.0
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    include(GoogleTest)
    add_subdirectory(tests)
endif()