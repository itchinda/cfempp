add_library(cfem_fields
    src/FESpace.cpp
    src/LagrangeSpace.cpp
    src/FEField.cpp
)

find_package(OpenMP REQUIRED)

# 2. Add the flags to your executable or library
if(OpenMP_CXX_FOUND)
    target_link_libraries(cfem_fields PUBLIC OpenMP::OpenMP_CXX)
endif()


target_include_directories(cfem_fields
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(cfem_fields
    PUBLIC cfem_utils cfem_mesh cfem_fem_integration
)

