add_library(cfem_fem_formulation
    src/Formulation.cpp
)

target_include_directories(cfem_fem_formulation
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_fem_formulation
    PUBLIC cfem_petsc_wrappers
)