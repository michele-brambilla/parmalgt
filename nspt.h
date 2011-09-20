#include "input.h"
#include "MyTime.h"

void gauge_wilson(ptGluon_fld&);

void fermion_wilson(ptGluon_fld&, ptSpinColor_fld&, SpinColor_fld&);

void zero_modes_subtraction(ptGluon_fld&);

void stochastic_gauge_fixing(ptGluon_fld&);
void FAstochastic_gauge_fixing(ptGluon_fld&);

// Quenched evolution
void NsptEvolve(ptGluon_fld&);

// Unquenched evolution
void NsptEvolve(ptGluon_fld&, ptSpinColor_fld&, SpinColor_fld&);

void plaquette_measure(ptGluon_fld&, act_params_t&);

int plaquette_check(Cplx*, Cplx*);

void WL2x2(ptGluon_fld&);
void ComputeLoops(ptGluon_fld&);


