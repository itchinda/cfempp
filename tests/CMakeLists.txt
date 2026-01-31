add_executable(unit_tests
    test_main.cpp
    wrappers/test_petsc_wrappers.cpp
    solvers/test_solvers.cpp
    fem/test_fem.cpp
    mesh/test_mesh.cpp
    utils/test_vector3d.cpp
    utils/test_static_vector.cpp
    utils/test_static_matrix.cpp
    utils/test_sym_tensor_o2.cpp
)

target_link_libraries(unit_tests
    PRIVATE
        GTest::gtest
        cfem_solvers
        cfem_utils
)

# This command tells CMake to find all TEST() macros in the binary
gtest_discover_tests(unit_tests)