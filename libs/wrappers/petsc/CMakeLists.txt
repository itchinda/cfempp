# libs/wrappers/petsc/CMakeLists.txt
add_library(cfem_petsc_wrappers 
    src/PetscMatrix.cpp 
    src/PetscVector.cpp 
    src/PetscKSPSolver.cpp
)

target_include_directories(cfem_petsc_wrappers PUBLIC include)
#target_link_libraries(cfem_petsc_wrappers PUBLIC ${PETSC_LIBRARIES})
target_link_libraries(cfem_petsc_wrappers PUBLIC PkgConfig::PETSC)

# Alias for easier linking later
add_library(odepp::wrappers::petsc ALIAS cfem_petsc_wrappers)