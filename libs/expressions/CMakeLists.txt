add_library(cfem_expressions
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

find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARY NAMES gmp)

target_link_libraries(cfem_expressions
    PUBLIC 
        cfem_utils 
        cfem_mesh 
        cfem_fem_integration 
        cfem_fields
        symengine
        mpfr
        mpc
        gmp
        OpenMP::OpenMP_CXX
        ${PETSC_LIBRARIES} # Use the variable if the imported target is acting up
        ${GMP_LIBRARY}
)