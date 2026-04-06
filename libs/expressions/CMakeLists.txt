add_library(cfem_expressions
    src/AlgebraicSimplifier.cpp
    src/CfemExpression.cpp
    src/BinaryExpression.cpp
    src/ConstantExpression.cpp
    src/GradExpression.cpp
    src/LambdaExpressionBase.cpp
    src/MathOperators.cpp
    src/UnaryExpression.cpp
    src/VectorExpression.cpp
    src/TensorExpression.cpp
    src/ZeroExpression.cpp
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
        cfem_fields
        # On utilise la cible IMPORTED, elle s'occupe de GMP, MPFR, et des includes !
        SymEngine::libsymengine 
        PkgConfig::PETSC
        OpenMP::OpenMP_CXX
)