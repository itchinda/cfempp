add_library(cfem_fields
    src/FESpace.cpp
)

target_include_directories(cfem_fields
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_fields
    PUBLIC cfem_mesh
)