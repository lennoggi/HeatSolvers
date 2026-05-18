#include <cassert>
#include <vector>
#include <string>
#include <iostream>

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
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCommand.h>

#include "include/Declare_functions.hh"
#include "include/SwitchIterationCallback.hh"
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
    warp->SetScaleFactor(0.01);
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

    // Toggle fullscreen view
    ////window->SetFullScreen(1);
    //const auto screen_size = window->GetScreenSize();
    //window->SetSize(screen_size[0], screen_size[1]);

    vtkNew<vtkRenderWindowInteractor> interactor;
    interactor->SetRenderWindow(window);

    vtkNew<vtkInteractorStyleTrackballCamera> istyle;
    interactor->SetInteractorStyle(istyle);


    // Add a custom callback to switch between iterations
    vtkNew<SwitchIterationCallback> switch_iter_cb;

    const auto sol_group_id = H5Gopen(file_id, "Solution", H5P_DEFAULT);
    assert(sol_group_id >= 0);

    H5G_info_t group_info;
    assert(H5Gget_info(sol_group_id, &group_info) >= 0);

    const auto niters = group_info.nlinks;  // NOTE: assumes the "Solution" group only contains datasets!
    assert(niters > 0);

    assert(H5Gclose(sol_group_id) >= 0);

    switch_iter_cb->file_id      = file_id;
    switch_iter_cb->current_iter = 0;
    switch_iter_cb->niters       = static_cast<size_t>(niters);
    switch_iter_cb->nx           = nx;
    switch_iter_cb->ny           = ny;
    switch_iter_cb->grid         = grid.Get();
    switch_iter_cb->warp         = warp.Get();
    switch_iter_cb->window       = window.Get();

    istyle->AddObserver(vtkCommand::KeyPressEvent, switch_iter_cb);

    interactor->Initialize();
    interactor->Start();


    /* --------
     * Close up
     * -------- */
    assert(H5Fclose(file_id) >= 0);

    return 0;
}
