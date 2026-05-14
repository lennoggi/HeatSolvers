#ifndef DECLARE_FUNCTIONS_HH
#define DECLARE_FUNCTIONS_HH

#include <vector>
#include <string>
#include <utility>

#include <hdf5.h>

#include <vtkDoubleArray.h>


std::pair<vtkSmartPointer<vtkDoubleArray>, std::vector<hsize_t>>
read_dset(const hid_t       &loc_id,
          const std::string &dset_name,
          const size_t      &dset_ndims);


#endif  // DECLARE_FUNCTIONS_HH
