#include "input.h"
#include "../MyTime.h"

void gauge_wilson(ptBoson_fld&);

/* void fermion_wilson(); */

void zero_modes_subtraction(ptBoson_fld&);

void stochastic_gauge_fixing(ptBoson_fld&);

// Quenched evolution
void NsptEvolve(ptBoson_fld&);

/* // Unquenched evolution */
/* void NsptEvolve(ptGluon_fld&, ptSpinColor_fld&, SpinColor_fld&, ptGluon_fld&); */

void plaquette_measure(ptBoson_fld&, act_params_t&);

int plaquette_check(Cplx*, Cplx*);

SU2 SU2rand(MyRand&);

#ifdef __WILLOOP_2x2__
void WL2x2(ptBoson_fld&);
#endif
