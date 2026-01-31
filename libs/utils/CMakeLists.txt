add_library(cfem_utils
    src/cfemutils.cpp
    src/ContractException.cpp
)

target_include_directories(cfem_utils
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_utils
    PUBLIC cfem_petsc_wrappers
)