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
        ${SYMENGINE_INCLUDE_DIRS}
)

target_link_libraries(cfem_expressions
    PUBLIC 
        cfem_utils
        cfem_mesh
        cfem_reference_element
        SymEngine::libsymengine
        #${llvm_libs}
        ${GMP_LIB}
        ${MPFR_LIB}
        ${MPC_LIB}
        # ${LLVM_SYSTEM_LIBS}
	cfem_openmp_interface
)
