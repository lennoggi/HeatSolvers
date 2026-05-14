#include <cassert>
#include <vector>
#include <string>
#include <utility>

#include <hdf5.h>

#include <vtkDoubleArray.h>

#include "include/Declare_functions.hh"

using namespace std;


/* ============================================================================
 * Routine reading a double-precision dataset of dimensionality dset_ndims into
 * a vtkDoubleArray
 * ============================================================================ */
pair<vtkSmartPointer<vtkDoubleArray>, vector<hsize_t>>
read_dset(const hid_t  &loc_id,
          const string &dset_name,
          const size_t &dset_ndims)
{
    const auto dset_id  = H5Dopen(loc_id, dset_name.c_str(), H5P_DEFAULT);
    assert(dset_id >= 0); 

    const auto space_id = H5Dget_space(dset_id);
    assert(space_id >= 0); 

    vector<hsize_t> size(dset_ndims);
    assert(H5Sget_simple_extent_dims(space_id, size.data(), nullptr) == dset_ndims);

    for (size_t i = 0; i < dset_ndims; ++i) {
        assert(size[i] > 0); 
    }

    auto arr = vtkSmartPointer<vtkDoubleArray>::New();
    arr->SetNumberOfComponents(1);  // Scalar field

    size_t size_flattened = size[0];

    for (size_t i = 1; i < dset_ndims; ++i) {
        size_flattened *= size[i];
    }   

    arr->SetNumberOfTuples(size_flattened);  // Allocates memory

    assert(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr->GetPointer(0)) >= 0); 
    assert(H5Dclose(dset_id) >= 0); 

    return {arr, size};
}
