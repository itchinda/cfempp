add_library(cfem_expressions
    src/AlgebraicSimplifier.cpp
    src/BinaryExpression.cpp
    src/GradExpression.cpp
    src/LambdaExpression.cpp
    src/UnaryExpression.cpp
    src/symengine/ScalarSymbolicExpression.cpp
    src/symengine/SymEngineSymbols.cpp
    src/symengine/TensorSymbolicExpression.cpp
    src/symengine/VectorSymbolicExpression.cpp
)

# Use the variables SymEngine provides explicitly
target_include_directories(cfem_expressions
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${SYMENGINE_INCLUDE_DIRS}  # <--- Explicitly add SymEngine headers        
)

target_link_libraries(cfem_expressions
    PUBLIC 
        cfem_utils 
        cfem_mesh 
        cfem_fem_integration
        # On utilise la cible IMPORTED, elle s'occupe de GMP, MPFR, et des includes !
        SymEngine::libsymengine 
        PkgConfig::PETSC
        OpenMP::OpenMP_CXX
)