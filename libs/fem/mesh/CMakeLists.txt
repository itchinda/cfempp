add_library(cfem_mesh
    src/builders/AnnulusMeshBuilder.cpp
    src/builders/BlockMeshBuilder.cpp
    src/builders/DiskMeshBuilder.cpp
    src/builders/HollowCylinderMeshBuilder.cpp
    src/builders/LineMeshBuilder.cpp
    src/builders/QuadMeshBuilder.cpp
    src/builders/SolidSphereBuilder.cpp
    src/builders/SurfaceSphereBuilder.cpp
    src/builders/CylinderMeshBuilder.cpp
    src/Mesh.cpp
    src/MeshEntity.cpp
    src/MeshConformityChecker.cpp
    src/GmshParser.cpp
    src/MeshQualityTool.cpp

)

target_include_directories(cfem_mesh
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include 
           ${CMAKE_CURRENT_SOURCE_DIR}/include/builders
)

target_link_libraries(cfem_mesh
    PUBLIC  cfem_utils
            cfem_openmp_interface
)
