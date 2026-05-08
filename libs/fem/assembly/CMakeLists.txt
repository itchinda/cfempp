add_subdirectory(terms)

add_library(cfem_assembly
    src/BilinearCoupling.cpp
    src/BoundaryWFTerm.cpp
    src/VariationalProblem.cpp
    src/VolumetricWFTerm.cpp
    src/WFTerm.cpp
)

target_include_directories(cfem_assembly
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_assembly
    PUBLIC cfem_petsc_wrappers
    cfem_utils
    cfem_openmp_interface
)
