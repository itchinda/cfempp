add_library(cfem_numerical_integration
    src/FunctionIntegrator.cpp
)

target_include_directories(cfem_numerical_integration
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_numerical_integration
    PUBLIC cfem_utils
           cfem_reference_element
           cfem_quadrature
           cfem_mesh
           cfem_fespace
           cfem_expressions
)