#ifndef DECLARE_FUNCTIONS_HH
#define DECLARE_FUNCTIONS_HH

#include <vector>


// Routines from Initial_data.cc
double gaussian(const double &x,
                const double &y);


// Routines from Utils.cc
std::vector<double> vecvec_lincomb(const std::vector<double> &v,
                                   const std::vector<double> &w,
                                   const double &a,
                                   const double &b);

double vecvec_dot(const std::vector<double> &v,
                  const std::vector<double> &w);


#endif  // DECLARE_FUNCTIONS_HH
