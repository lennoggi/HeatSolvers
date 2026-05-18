#include <sstream>
#include <stdexcept>

#include <hdf5.h>

#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkWarpScalar.h>
#include <vtkRenderWindow.h>

#include "include/Declare_functions.hh"

using namespace std;


/* =========================================================================
 * Routine loading a new iteration and updating the VTK pipeline accordingly
 * ========================================================================= */
void load_iteration(const size_t &it,
                    const size_t &nx,
                    const size_t &ny,
                    const hid_t  &file_id,
                    vtkRectilinearGrid *grid,
                    vtkWarpScalar      *warp,
                    vtkRenderWindow    *window) {
    ostringstream it_ss;
    it_ss << "Solution/Iteration " << it;

    if (file_id < 0) {
        throw runtime_error("load_iteration(): got invalid file_id");
    }

    auto [u, udims] = read_dset(file_id, it_ss.str(), 2);

    if (udims[0] != nx or udims[1] != ny) {
        ostringstream err_ss;
        err_ss  << "load_iteration(): unexpected solution dataset dimesions ("
                << udims[0] << ", " << udims[1] << ") at iteration " << it
                << "(expected dimensions: (" << nx << ", " << ny << ")";
        throw runtime_error(err_ss.str());
    }

    if (grid == nullptr) {
        throw runtime_error("load_iteration(): got invalid VTK grid");
    }

    grid->GetPointData()->SetScalars(u);
    grid->Modified();
    warp->Update();
    window->Render();

    return;
}
