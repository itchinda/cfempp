add_library(cfem_fem_integration
    src/ReferenceElementFactory.cpp
)

target_include_directories(cfem_fem_integration
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_fem_integration
    PUBLIC cfem_petsc_wrappers cfem_mesh
)