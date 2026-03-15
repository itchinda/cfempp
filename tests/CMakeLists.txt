add_executable(unit_tests
    test_main.cpp
    wrappers/test_petsc_wrappers.cpp
    solvers/test_solvers.cpp
    expressions/test_symbolic_simplification.cpp
    fem/test_fem.cpp
    fem/test_reference_element.cpp
    fem/test_quadrature.cpp
    fem/test_fields.cpp
    mesh/test_mesh.cpp
    mesh/test_gmsh_parser.cpp
    utils/test_dynamic_vector.cpp
    utils/test_vector3d.cpp
    utils/test_fixed_vector.cpp
    utils/test_fixed_matrix.cpp
    utils/test_sym_tensor_o2.cpp
)

target_link_libraries(unit_tests
    PRIVATE
        GTest::gtest
        cfem_solvers
        cfem_utils
        cfem_mesh
        cfem_fem_integration
        cfem_fields
        cfem_expressions
        SymEngine::libsymengine 
        ${llvm_libs}
        ${GMP_LIB}
        ${MPFR_LIB}
        ${LLVM_SYSTEM_LIBS}
)

# This command tells CMake to find all TEST() macros in the binary
gtest_discover_tests(unit_tests)