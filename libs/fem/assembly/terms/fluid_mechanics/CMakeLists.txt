add_library(cfem_fluid_mechanics_wfterm
    src/FluidViscosityWFTerm.cpp
)

target_include_directories(cfem_fluid_mechanics_wfterm
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_fluid_mechanics_wfterm
    PUBLIC 
        cfem_utils
        cfem_petsc_wrappers
        cfem_reference_element
        cfem_expressions
        cfem_mesh
        cfem_fespace
)