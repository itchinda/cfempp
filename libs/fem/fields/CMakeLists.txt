add_library(cfem_fields
    src/FESpace.cpp
    src/LagrangeSpace.cpp
    src/FEField.cpp
    src/DomainIntegrator.cpp
)

# Use the variables SymEngine provides explicitly
target_include_directories(cfem_fields
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${SYMENGINE_INCLUDE_DIRS}  # <--- Explicitly add SymEngine headers
)

# find_path(GMP_INCLUDE_DIR NAMES gmp.h)
# find_library(GMP_LIBRARY NAMES gmp)

# target_link_libraries(cfem_fields
#     PUBLIC 
#         cfem_utils 
#         cfem_mesh 
#         cfem_fem_integration
#         cfem_expressions
#         symengine
#         mpfr
#         mpc
#         gmp
#         OpenMP::OpenMP_CXX
#         ${PETSC_LIBRARIES} # Use the variable if the imported target is acting up
#         ${GMP_LIBRARY}
# )

target_link_libraries(cfem_fields
    PUBLIC 
        cfem_utils
        cfem_mesh
        cfem_fem_integration
        cfem_expressions
        # On utilise la cible IMPORTED, elle s'occupe de GMP, MPFR, et des includes !
        SymEngine::libsymengine 
        PkgConfig::PETSC
        OpenMP::OpenMP_CXX
)