AM_CPPFLAGS = -DSF -DSFBC -D__PARALLEL_OMP__ -DDSFMT_MEXP=19937 \
  -DHAVE_TR1

noinst_HEADERS = Background.h nspt.h MassC.h input.h MyTime.h \
  QCDenvNODEpt.h QCDenvNODE.h Point.hpp newQCDpt.h lattice.h \
  PtTypes.hpp MyMath.h MyQCD.h MyRand.h Types.h dSFMT/dSFMT.h \
  choices.h PtAlgo.hpp GammaStuff.h

dist_data_DATA = damocle.dag Quench.cfg run.sh

bin_PROGRAMS = Quenched

Quenched_SOURCES = Quenched.cc MyMath.cc MyRand.cc nspt.cc \
  newQCDpt.cc MyQCD.cc Point.cc Background.cc
Quenched_LDADD = -lfftw3
Quenched_CXXFLAGS = -fopenmp
