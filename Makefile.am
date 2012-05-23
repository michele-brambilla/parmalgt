AM_CPPFLAGS = -DSF -DSFBC -D__PARALLEL_OMP__ -DDSFMT_MEXP=19937 \
  -DHAVE_TR1
ACLOCAL_AMFLAGS = -I m4
noinst_HEADERS = Background.h nspt.h MassC.h input.h MyTime.h \
  QCDenvNODEpt.h QCDenvNODE.h Point.hpp newQCDpt.h lattice.h \
  PtTypes.hpp MyMath.h MyQCD.h MyRand.h Types.h dSFMT/dSFMT.h \
  choices.h PtAlgo.hpp GammaStuff.h

dist_data_DATA = damocle.dag Quench.cfg run.sh

bin_PROGRAMS = LocalQuench

#ParQuenched_SOURCES = ParQuenched.cc MyMath.cc MyRand.cc \
  newQCDpt.cc Point.cc Background.cc

LocalQuench_SOURCES = LocalQuench.cc MyMath.cc MyRand.cc \
  newQCDpt.cc Point.cc Background.cc

.PHONY: docs

docs:
	mkdir -p ./docs/doxygen
	doxygen Doxyfile

