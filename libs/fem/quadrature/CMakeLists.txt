set(QUADRATURE_HEADERS
    include/Quadrature.h
    include/Quadrature.tpp
)

add_library(cfem_quadrature
    src/Quadrature.cpp
    ${QUADRATURE_HEADERS}
)

target_include_directories(cfem_quadrature
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(cfem_quadrature
    PUBLIC 
        cfem_utils
)