add_library(cfem_fem_assembly
    src/Assembly.cpp
)

target_include_directories(cfem_fem_assembly
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_fem_assembly
    PUBLIC cfem_petsc_wrappers
)