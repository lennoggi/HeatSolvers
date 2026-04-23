#include <cmath>
#include "include/Declare_functions.hh"
#include "Parameters.hh"

using namespace std;


// Gaussian
double gaussian(const double &x,
                const double &y) {
    constexpr double norm = 1.0/(2.0*M_PI*SIGMAX*SIGMAY);

    constexpr double expfacx = 1.0/(2.0*SIGMAX*SIGMAX);
    constexpr double expfacy = 1.0/(2.0*SIGMAY*SIGMAY);

    const double xs = x - X0;
    const double ys = y - Y0;

    //auto gauss = norm*exp(-expfacx*xs*xs - expfacy*ys*ys);

    //if (gauss < 1.0e-16) { gauss = 0.0; }

    return norm*exp(-expfacx*xs*xs - expfacy*ys*ys);
}
