#ifndef CHECK_PARAMETERS_HH
#define CHECK_PARAMETERS_HH

#include "../Parameters.hh"


static_assert(SCALE == LIN or SCALE == LOG);
static_assert(CBAR_MAX > CBAR_MIN);

#if SCALE == LOG
static_assert(CBAR_MIN > 0.0);
#endif


#endif  // CHECK_PARAMETERS_HH
