#ifndef PARAMETERS_HH
#define PARAMETERS_HH


// ***** DON'T TOUCH *****
#define GAUSSIAN 0
// ***********************

// Initial data
#define INITIAL_DATA GAUSSIAN


// Diffusion constant
constexpr inline double ALPHA = 1.0;

// Left and right time-independent function values (Dirichlet BCs)
constexpr inline double UX0 = 0.0;
constexpr inline double UXL = 0.0;
constexpr inline double UY0 = 0.0;
constexpr inline double UYL = 0.0;


// Number of spatial grid points
constexpr inline size_t NX = 1000;
constexpr inline size_t NY = 1000;

// Spatial bounds are [0, L - L/N] along both x and y
constexpr inline double LX = 1.0;
constexpr inline double LY = 1.0;


// Number of time steps
constexpr inline size_t NT = 200;

// Time step
constexpr inline double DT = 1.0e-04;


// Initial function center
constexpr inline double X0 = 0.3*LX;
constexpr inline double Y0 = 0.3*LY;

// Initial wavefunction half-width
constexpr inline double SIGMAX = 0.03*LX;
constexpr inline double SIGMAY = 0.03*LY;


/* The conjugate gradient method stops after the relative residual
 *   r = || Au - u_old || / || u_old || ,
 * where u is the solution and A is the 2D heat equation kernel, falls below
 * this threshold                                                               */
constexpr inline double CGTOL = 1.e-06;

// Output frequency
constexpr inline size_t OUT_EVERY = 1;

// Output filename
#define FILENAME "Data.h5"

// Verbosity
constexpr inline bool     VERBOSE = true;
constexpr inline bool VERYVERBOSE = false;


#endif  // PARAMETERS_HH
