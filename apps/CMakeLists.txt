# apps/CMakeLists.txt

# Define the main simulation executable
# add_executable(${APP_PROJECT_NAME} 
#     main.cpp
# )

# Link against the solver library
# This automatically pulls in fem, mesh, wrappers, and PETSc
# target_link_libraries(${APP_PROJECT_NAME}
#     PRIVATE
#         cfem_solvers
#         cfem_utils
#         cfem_mesh
#         cfem_fields
#         cfem_exporter
#         cfem_expressions
#         cfem_petsc_wrappers
#         # --- LA RÉSOLUTION DU PROBLÈME EST ICI ---
#         SymEngine::libsymengine
#         ${llvm_libs}
#         dl     # Souvent requis par LLVM sur Linux
#         z      # Souvent requis par LLVM sur Linux
# )

add_executable(${APP_PROJECT_NAME} main.cpp)

target_link_libraries(${APP_PROJECT_NAME}
    PRIVATE
        cfem_solvers
        cfem_utils
        cfem_mesh
        cfem_fefield
        cfem_expressions
        cfem_exporter
        cfem_numerical_integration
        # SymEngine::libsymengine   # 1. La lib statique qui demande des symboles
        # ${llvm_libs}              # 2. Les libs LLVM qui fournissent les symboles
        # ${GMP_LIB}                # 3. Les libs de calcul
        # ${MPFR_LIB}
        # ${MPC_LIB}
        # ${LLVM_SYSTEM_LIBS}       # 4. Les libs système de bas niveau
)