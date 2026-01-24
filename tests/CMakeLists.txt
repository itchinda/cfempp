add_executable(unit_tests
    test_main.cpp
    wrappers/test_petsc_wrappers.cpp
    solvers/test_solvers.cpp
    fem/test_fem.cpp
    mesh/test_mesh.cpp
)

target_link_libraries(unit_tests
    PRIVATE
        GTest::gtest
        cfem_solvers       # This pulls in fem, mesh, and wrappers automatically
)

# This command tells CMake to find all TEST() macros in the binary
gtest_discover_tests(unit_tests)