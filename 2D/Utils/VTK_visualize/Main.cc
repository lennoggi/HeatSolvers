// NOTE: included here to have SCALE correctly imported
#include "Parameters.hh"

#include <cassert>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include <hdf5.h>

#include <vtkDoubleArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkWarpScalar.h>
#include <vtkDataSetMapper.h>

#include <vtkActor.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#if SCALE == LIN
#include <vtkLookupTable.h>
#elif SCALE == LOG
#include <vtkLogLookupTable.h>
#else
#error "Invalid color scale"
#endif

#include <vtkScalarBarActor.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include <vtkCommand.h>

#include "include/Declare_functions.hh"
#include "include/SwitchIterationCallback.hh"
#include "include/Check_parameters.hh"


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


    // Add informative text
    ostringstream text_ss;
    text_ss << "Iteration 0" << endl
            << "Left/Right: switch iteration" << endl
            << "<Q>/<E>: quit" << endl
            << "Left mouse button: rotate the camera around its focal point" << endl
            << "Shift + left mouse button: pan the camera" << endl
            << "Ctrl + left mouse button: roll the camera around its view axis" << endl
            << "Ctrl + shift + left mouse button / Right mouse button / Mouse wheel rolling: zoom in/out";

    vtkNew<vtkTextActor> text_actor;
    text_actor->SetInput(text_ss.str().c_str());
    text_actor->GetTextProperty()->SetFontSize(18);
    text_actor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);  // White
    text_actor->SetPosition(10, 10);  // Pixels from bottom-left


    // Add a color bar
    #if SCALE == LIN
    vtkNew<vtkLookupTable> lookup_table;
    #elif SCALE == LOG
    vtkNew<vtkLogLookupTable> lookup_table;
    #else
    #error "Invalid color scale"
    #endif

    lookup_table->SetTableRange(CBAR_MIN, CBAR_MAX);
    lookup_table->Build();

    vtkNew<vtkScalarBarActor> colorbar_actor;
    colorbar_actor->SetNumberOfLabels(7);  // XXX XXX XXX XXX XXX XXX
    colorbar_actor->SetLabelFormat("%.3e");
    colorbar_actor->SetLookupTable(lookup_table);

    // Make sure the data also use the same color bar
    mapper->SetLookupTable(lookup_table);
    mapper->SetScalarRange(CBAR_MIN, CBAR_MAX);
    mapper->SetColorModeToMapScalars();
    mapper->ScalarVisibilityOn();


    // Add actors to the renderer and set up the rendering window
    vtkNew<vtkRenderer> renderer;
    renderer->AddActor(surface_actor);
    renderer->AddActor2D(text_actor);
    renderer->AddActor2D(colorbar_actor);

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
    switch_iter_cb->text_actor   = text_actor.Get();
    switch_iter_cb->grid         = grid.Get();
    switch_iter_cb->warp         = warp.Get();
    switch_iter_cb->window       = window.Get();

    istyle->AddObserver(vtkCommand::KeyPressEvent, switch_iter_cb);


    // Start
    interactor->Initialize();
    interactor->Start();


    /* --------
     * Close up
     * -------- */
    assert(H5Fclose(file_id) >= 0);

    return 0;
}
