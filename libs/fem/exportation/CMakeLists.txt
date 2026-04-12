add_library(cfem_exporter
    src/MeshExporter.cpp
    src/VTKExporter.cpp
    src/VTUExporter.cpp
)

target_include_directories(cfem_exporter
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_exporter
    PUBLIC cfem_mesh
           cfem_expressions
           cfem_fefield
           cfem_reference_element
)