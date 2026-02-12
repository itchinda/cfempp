# apps/CMakeLists.txt

# Define the main simulation executable
add_executable(${APP_PROJECT_NAME} 
    main.cpp
)

# Link against the solver library
# This automatically pulls in fem, mesh, wrappers, and PETSc
target_link_libraries(${APP_PROJECT_NAME}
    PRIVATE
        cfem_solvers
        cfem_utils
        cfem_mesh
        cfem_fields
        cfem_exporter
        cfem_petsc_wrappers
)

# Optional: Add more apps as you go
# add_executable(wave_equation_2d wave_2d.cpp)
# target_link_libraries(wave_equation_2d PRIVATE cfem_solvers)