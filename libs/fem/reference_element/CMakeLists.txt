add_library(cfem_reference_element
    src/ReferenceElementFactory.cpp
)

target_include_directories(cfem_reference_element
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_reference_element
    PUBLIC cfem_utils cfem_fespace
)