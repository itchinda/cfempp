add_library(cfem_solvers
    src/Solvers.cpp
)

target_include_directories(cfem_solvers
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_solvers
    PUBLIC cfem_petsc_wrappers
)