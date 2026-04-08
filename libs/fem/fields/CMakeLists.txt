add_library(cfem_fields
    src/FESpace.cpp
    src/LagrangeSpace.cpp
    src/DGSpace.cpp
    src/FEField.cpp
    src/Integrator.cpp
)

target_include_directories(cfem_fields
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${SYMENGINE_INCLUDE_DIRS}
)

target_link_libraries(cfem_fields
    PUBLIC 
        cfem_utils
        cfem_mesh
        cfem_fem_integration
        cfem_expressions
        SymEngine::libsymengine 
        PkgConfig::PETSC
        OpenMP::OpenMP_CXX
)