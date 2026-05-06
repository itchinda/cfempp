add_library(cfem_fefield
    src/FEField.cpp

)

target_include_directories(cfem_fefield
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${SYMENGINE_INCLUDE_DIRS}
)

target_link_libraries(cfem_fefield
    PUBLIC 
        cfem_utils
        cfem_petsc_wrappers
        cfem_reference_element
        cfem_expressions
        cfem_mesh
        cfem_fespace
)