include ../makefile.local

LIBs = ../MyRand.o ../MyMath.o ../MyQCD.o ../QCDpt.o 
INCs = ../MyRand.h ../MyMath.h ../MyQCD.h ../QCDpt.h ../lattice.h ../QCDenvNODEpt.h


test: ptjnk accessi verifica

ptjnk:  $(LIBs) $(INCs) ptjnk.cc
	$(CC) $(MyO) -o ptjnk.exe PTjnk.cc $(LIBs) -I../

accessi: $(LIBs) $(INCs) accessi.cc
#	$(CC) $(MyO) -o accessi.exe accessi.cc $(LIBs) -I../
	$(CC) $(MyO) -o accessi3.exe accessi3.cc $(LIBs) -I../

verifica: verifica.cc $(LIBs) $(INCs)
	$(CC) $(MyO) -o verifica.exe verifica.cc $(LIBs) -I../

clean: 
	rm *.exe; rm *.o