add_library(cfem_kernel
    src/EvalContext.cpp
)

target_include_directories(cfem_kernel
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(cfem_kernel
    PUBLIC cfem_utils
)