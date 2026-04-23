#ifndef CHECK_PARAMETERS_HH
#define CHECK_PARAMETERS_HH

#include "../Parameters.hh"


static_assert(INITIAL_DATA == GAUSSIAN);

static_assert(ALPHA > 0.0);

static_assert(NX > 1);
static_assert(NY > 1);

static_assert(LX > 0.0);
static_assert(LY > 0.0);

static_assert(NT > 1, "Need at least two time steps for correct output handling");
static_assert(DT > 0.0);

static_assert(SIGMAX > 0.0);
static_assert(SIGMAY > 0.0);

static_assert(X0 - 3.0*SIGMAX > 0.0 and X0 + 3.0*SIGMAX < LX - LX/static_cast<double>(NX),
              "Initial function is quite out of bounds along x");
static_assert(Y0 - 3.0*SIGMAY > 0.0 and Y0 + 3.0*SIGMAY < LY - LY/static_cast<double>(NY),
              "Initial function is quite out of bounds along y");

static_assert(CGTOL > 0.0 and CGTOL < 1.0);

static_assert(OUT_EVERY > 0);


static_assert(VERBOSE     or not     VERBOSE);
static_assert(VERYVERBOSE or not VERYVERBOSE);
static_assert(not VERYVERBOSE or VERBOSE);  // VERYVERBOSE allowed only if VERBOSE is true


#endif  // CHECK_PARAMETERS_HH
