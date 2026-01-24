cmake_minimum_required(VERSION 3.15)

# ------------------------------------------------------------
# Build type / PETSc configuration
# ------------------------------------------------------------
# Compilation modes:
#   cmake .. -DCMAKE_BUILD_TYPE=Release
#   cmake .. -DCMAKE_BUILD_TYPE=Debug
# ------------------------------------------------------------
# 1. Build type / PETSc configuration (Must happen before project())
# ------------------------------------------------------------
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
# 2. Project Initialization
# ------------------------------------------------------------
project(cfem LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Define the App name here so it's available in the subdirectories
if(TARGET_PETSC_ARCH STREQUAL "arch-linux-c-debug")
    set(APP_PROJECT_NAME cfem_app_dev)
else()
    set(APP_PROJECT_NAME cfem_app_opt)
endif()

# Global PETSc target available to all subdirectories
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET petsc)

# ------------------------------------------------------------
# 3. Add Modular Libraries & Applications
# ------------------------------------------------------------
add_subdirectory(libs)
add_subdirectory(apps) # This will look for apps/CMakeLists.txt

# ------------------------------------------------------------
# 4. Testing / GoogleTest
# ------------------------------------------------------------
include(CTest)

if(BUILD_TESTING)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG         v1.17.0
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    include(GoogleTest)
    add_subdirectory(tests)
endif()