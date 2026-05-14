cmake_minimum_required(VERSION 3.16)

# ----------
# Parameters
# ----------
set(PROJ "Heat2D_viz")
set(EXE  ${PROJ})
set(SRC Main.cc Read_dset.cc)


project(${PROJ} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel  # For vtkStructuredGrid
    RenderingCore
    FiltersGeometry
    InteractionStyle
    RenderingFreeType
    RenderingOpenGL2
##    RenderingAnnotation   # vtkCubeAxesActor, vtkScalarBarActor
##    InteractionWidgets    # vtkSliderWidget
##    FiltersGeometry       # vtkStructuredGridGeometryFilter
##    FiltersCore           # vtkPolyDataNormals
)

add_executable(${EXE} ${SRC})

target_link_libraries(${EXE} PRIVATE HDF5::HDF5)
target_link_libraries(${EXE} PRIVATE ${VTK_LIBRARIES})

vtk_module_autoinit(
    TARGETS ${EXE}
    MODULES ${VTK_LIBRARIES}
)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

install(TARGETS ${EXE}
        RUNTIME DESTINATION bin
)
