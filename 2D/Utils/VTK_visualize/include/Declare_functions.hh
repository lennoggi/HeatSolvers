#ifndef DECLARE_FUNCTIONS_HH
#define DECLARE_FUNCTIONS_HH

#include <vector>
#include <string>
#include <utility>

#include <hdf5.h>

#include <vtkDoubleArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkWarpScalar.h>
#include <vtkRenderWindow.h>


std::pair<vtkSmartPointer<vtkDoubleArray>, std::vector<hsize_t>>
read_dset(const hid_t       &loc_id,
          const std::string &dset_name,
          const size_t      &dset_ndims);

void load_iteration(const size_t &it,
                    const size_t &nx,
                    const size_t &ny,
                    const hid_t  &file_id,
                    vtkRectilinearGrid *grid,
                    vtkWarpScalar      *warp,
                    vtkRenderWindow    *window);


#endif  // DECLARE_FUNCTIONS_HH
