#include <cassert>
#include <vector>
#include <string>

#include <hdf5.h>

#include <vtkDoubleArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkWarpScalar.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleJoystickCamera.h>

#include "include/Declare_functions.hh"
#include "Parameters.hh"


using namespace std;


int main() {
    const auto file_id = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(file_id >= 0);

    /* -----------------------------------------
     * Load the coordinates and the initial data
     * ----------------------------------------- */
    auto [x, xdims] = read_dset(file_id, "Coordinates/x", 1);
    auto [y, ydims] = read_dset(file_id, "Coordinates/y", 1);
    auto [u, udims] = read_dset(file_id, "Solution/Iteration 0", 2);

    assert(xdims.size() == 1);
    assert(ydims.size() == 1);
    assert(udims.size() == 2);

    const auto nx = xdims[0];
    const auto ny = ydims[0];

    assert(udims[0] == nx);
    assert(udims[1] == ny);


    /* ------------------------------------------------
     * Build the grid and fill it with the initial data
     * ------------------------------------------------ */
    vtkNew<vtkRectilinearGrid> grid;
    grid->SetDimensions(static_cast<int>(nx), static_cast<int>(ny), 1);

    grid->SetXCoordinates(x);
    grid->SetYCoordinates(y);

    grid->GetPointData()->SetScalars(u);


    /* --------------------------------
     * Build the visualization pipeline
     * -------------------------------- */
    /* Warp the 2D data into a filter elevating the 2D grid cells based on
     * their scalar value                                                       */
    vtkNew<vtkWarpScalar> warp;
    warp->SetInputData(grid);
    warp->SetScaleFactor(1.0);   // adjust this to exaggerate or flatten the surface
    warp->UseNormalOn();
    warp->SetNormal(0.0, 0.0, 1.0);
    warp->Update();

    vtkNew<vtkDataSetMapper> mapper;
    //mapper->SetInputData(grid);  // Good for color map visualization
    mapper->SetInputConnection(warp->GetOutputPort());

    vtkNew<vtkActor> surface_actor;
    surface_actor->SetMapper(mapper);

    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(surface_actor);

    vtkNew<vtkRenderWindow> window;
    window->AddRenderer(renderer);

    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(window);

    vtkNew<vtkInteractorStyleJoystickCamera> istyle;

    interactor->Initialize();
    interactor->Start();


    /* --------
     * Close up
     * -------- */
    assert(H5Fclose(file_id) >= 0);

    return 0;
}
