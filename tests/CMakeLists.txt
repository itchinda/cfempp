add_executable(unit_tests
    test_main.cpp
    wrappers/test_petsc_wrappers.cpp
    solvers/test_solvers.cpp
    expressions/test_symbolic_simplification.cpp
    fem/test_fem.cpp
    fem/test_domain_integrator.cpp
    fem/test_reference_element.cpp
    fem/test_quadrature.cpp
    fem/test_fields.cpp
    fem/test_fefield.cpp
    mesh/test_mesh.cpp
    mesh/test_gmsh_parser.cpp
    utils/test_dynamic_matrix.cpp
    utils/test_dynamic_vector.cpp
    utils/test_vector3d.cpp
    utils/test_fixed_vector.cpp
    utils/test_fixed_matrix.cpp
    utils/test_sym_tensor_o2.cpp
    utils/tensors/test_tensor_addition.cpp
    utils/tensors/test_tensor_product.cpp
    utils/tensors/test_tensor_substraction.cpp
)

target_link_libraries(unit_tests
    PRIVATE
        GTest::gtest
	    cfem_openmp_interface
        cfem_solvers
        cfem_utils
        cfem_kernel
        cfem_quadrature
        cfem_reference_element
        cfem_mesh
        cfem_fespace
        cfem_fefield
        cfem_expressions
        cfem_numerical_integration
)

target_compile_definitions(unit_tests PRIVATE CFEM_TESTING)

gtest_discover_tests(unit_tests)
