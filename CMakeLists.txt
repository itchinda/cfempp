cmake_minimum_required(VERSION 3.15)

# ------------------------------------------------------------
# Build type / PETSc configuration
# ------------------------------------------------------------
# Compilation modes:
#   cmake .. -DCMAKE_BUILD_TYPE=Release
#   cmake .. -DCMAKE_BUILD_TYPE=Debug
if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(TARGET_PETSC_ARCH "arch-linux-c-debug")
    add_definitions(-DPETSC_USE_DEBUG)
else()
    set(TARGET_PETSC_ARCH "arch-linux-c-opt")
endif()

# Force PETSC_ARCH to match CMake build type
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
# Project
# ------------------------------------------------------------
project(odepp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ------------------------------------------------------------
# Library
# ------------------------------------------------------------
add_library(odepp_lib
    src/PetscMatrix.cpp
    src/PetscVector.cpp
    src/PetscKSPSolver.cpp
)

target_include_directories(odepp_lib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

pkg_search_module(PETSC REQUIRED IMPORTED_TARGET petsc)

target_link_libraries(odepp_lib
    PUBLIC
        PkgConfig::PETSC
)

# ------------------------------------------------------------
# Application
# ------------------------------------------------------------
if(TARGET_PETSC_ARCH STREQUAL "arch-linux-c-debug")
    set(APP_PROJECT_NAME odepp_app_dev)
else()
    set(APP_PROJECT_NAME odepp_app_opt)
endif()

add_executable(${APP_PROJECT_NAME}
    apps/main.cpp
)

target_link_libraries(${APP_PROJECT_NAME}
    PRIVATE
        odepp_lib
)

# ------------------------------------------------------------
# Testing / GoogleTest
# ------------------------------------------------------------
include(CTest)   # enables BUILD_TESTING + enable_testing()

if(BUILD_TESTING)

    include(FetchContent)

    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.17.0
    )

    # Prevent overriding parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(googletest)

    add_executable(unit_tests
        tests/test_petsc_solver.cpp
    )

    target_link_libraries(unit_tests
        PRIVATE
            odepp_lib
            GTest::gtest_main
    )

    add_test(NAME lib_gtest COMMAND unit_tests)

endif()
