# Source files
get_source_files(SOURCES)
append_source_files(SOURCES Adsorption)
append_source_files(SOURCES SolidModels)
append_source_files(SOURCES FractureModels)
append_source_files(SOURCES FractureModels/Permeability)

append_source_files(SOURCES Fluid)
append_source_files(SOURCES Fluid/Density)
append_source_files(SOURCES Fluid/Viscosity)
append_source_files(SOURCES Fluid/GibbsFreeEnergy)
append_source_files(SOURCES Fluid/FluidProperties)
append_source_files(SOURCES Fluid/SpecificHeatCapacity)
append_source_files(SOURCES Fluid/ThermalConductivity)
append_source_files(SOURCES Fluid/WaterVaporProperties)

append_source_files(SOURCES MPL)
append_source_files(SOURCES MPL/Properties)
append_source_files(SOURCES MPL/Components)
append_source_files(SOURCES MPL/Utils)

append_source_files(SOURCES PorousMedium)
append_source_files(SOURCES PorousMedium/Porosity)
append_source_files(SOURCES PorousMedium/Storage)
append_source_files(SOURCES PorousMedium/Permeability)
append_source_files(SOURCES PorousMedium/UnsaturatedProperty/CapillaryPressure)
append_source_files(SOURCES
                    PorousMedium/UnsaturatedProperty/RelativePermeability)
append_source_files(SOURCES TwoPhaseModels)

if(OGS_USE_MFRONT)
    add_subdirectory(SolidModels/MFront)
endif()

add_library(MaterialLib ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS MaterialLib LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

include(GenerateExportHeader)
generate_export_header(MaterialLib)
target_include_directories(MaterialLib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(MaterialLib PRIVATE MathLib MeshLib ParameterLib)

if(OGS_USE_MFRONT)
    target_link_libraries(MaterialLib PUBLIC MaterialLib_SolidModels_MFront)
endif()

if(OGS_USE_PCH)
    cotire(MaterialLib)
endif()
