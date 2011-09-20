include ../makefile.local

DEF = -D_VERBOSE_
#DEF = -D_NO_VERBOSE_

LIBs = ../rand55.o ../MyMath.o ../MyQCD.o ../QCDpt.o
INCs = ../rand55.h ../MyMath.h ../MyQCD.h ../QCDpt.h ../lattice.h ../QCDenvNODEpt.h dirac_spec.h \
	dirac_par.h

diractest: $(LIBs) $(INCs) diractest.o dirac_spec.o
	$(CC) $(MyO) $(DEF) -lgsl -lgslcblas -o diractest.exe diractest.o dirac_spec.o $(LIBs) -I./ -I../  -I$(FFTWdir)/include/  -I$(GSLdir) -L$(FFTWdir)/lib/ -L$(GSLdir) -lfftw3

diractest.o: diractest.cc
	$(CC) $(MyO) $(DEF) -c diractest.cc -I./ -I../ -I$(FFTWdir)/include/ -I$(GSLdir)

dirac_spec.o: dirac_spec.cc dirac_spec.h
	$(CC) $(MyO) $(DEF) -c dirac_spec.cc -I./ -I../ -I$(FFTWdir)/include/ -I$(GSLdir)

test_lus: $(LIBs) $(INCs) test_lus.cc
	$(CC) $(MyO) $(DEF) -lfftw3 -o test_lus.exe test_lus.cc $(LIBs) -I./ -I../

test_dd: $(LIBs) $(INCs) test_dd.cc
	$(CC) $(MyO) $(DEF) -lfftw3 -o test_dd.exe test_dd.cc $(LIBs) -I./ -I../

.PHONY : clean

clean:
	\rm *.exe *.o
