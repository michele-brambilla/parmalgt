SHELL = /bin/sh

.SUFFIXES: 
.SUFFIXES: .cc .o

include ../makefile.local

objects := $(filter-out ../nspt.o ../Quenched.o ../Unquenched.o ../MyTime.o, $(patsubst %.cc,%.o,$(wildcard ../*.cc RenLib.cc)))
mains := $(patsubst %.cc,%.o,NewRenorm.cc GaugeFixing.cc)
exe := $(patsubst %.o,%,$(mains))


all: $(exe) 

$(objects): %.o : %.cc
	$(CC) $(MyO) -c $< -I. -I../ -I$(FFTWdir)include $(PARALLEL_WLOOP) 

$(mains): %.o : %.cc
	$(CC) $(MyO) -c $< -I. -I../ -I$(FFTWdir)include $(PARALLEL_WLOOP) 

$(exe): $(mains) $(objects)
	$(CC) $(MyO) -o $@.exe $@.o $(objects)  -I. -I../ -L$(FFTWdir)lib -lfftw3 -L. $(PARALLEL_WLOOP) 

include $(objects:.o=.d)
include $(mains:.o=.d)

%.d: %.cc 
	@set -e; rm -f $@; \
	$(CC) -MM $(PARALLEL_WLOOP) $(MyO) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

.PHONY : clean

clean:
	\rm *.exe *.o 
