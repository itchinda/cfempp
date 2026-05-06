add_library(cfem_fespace
    src/FESpace.cpp
    src/LagrangeSpace.cpp
    src/DGSpace.cpp
)

target_include_directories(cfem_fespace
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${SYMENGINE_INCLUDE_DIRS}
)

target_link_libraries(cfem_fespace
    PUBLIC 
        cfem_utils
        cfem_reference_element
        cfem_mesh
        cfem_expressions
        cfem_fefield
)