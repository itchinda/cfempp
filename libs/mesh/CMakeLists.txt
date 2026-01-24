add_library(cfem_mesh
    src/Mesh.cpp
)

target_include_directories(cfem_mesh
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_mesh
    PUBLIC cfem_petsc_wrappers
)