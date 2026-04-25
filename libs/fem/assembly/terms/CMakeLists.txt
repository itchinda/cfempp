option(BUILD_GENERAL_TERMS "Build General Terms module" ON)
option(BUILD_SOLID_MECHANICS_TERMS "Build Solid Mechanics module" ON)
option(BUILD_FLUID_MECHANICS_TERMS "Build Fluid Mechanics module" ON)
option(BUILD_ELECTROMAGNETISM_TERMS "Build Electromagnetism module" ON)

if(BUILD_GENERAL_TERMS)
    add_subdirectory(general)
endif()

if(BUILD_SOLID_MECHANICS_TERMS)
    add_subdirectory(solid_mechanics)
endif()

if(BUILD_FLUID_MECHANICS_TERMS)
    add_subdirectory(fluid_mechanics)
endif()

if(BUILD_ELECTROMAGNETISM_TERMS)
    add_subdirectory(electromagnetism)
endif()