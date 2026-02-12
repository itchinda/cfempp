add_library(cfem_mesh
    src/Mesh.cpp
    src/MeshEntity.cpp
    src/GmshParser.cpp
)

target_include_directories(cfem_mesh
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_mesh
    PUBLIC cfem_petsc_wrappers cfem_utils cfem_fields
)