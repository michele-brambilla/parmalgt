SHELL = /bin/sh
VPATH = $(FFTWdir)lib

.SUFFIXES: 
.SUFFIXES: .cc .o

include ../makefile.local

mains   := $(patsubst %.cc,%.o,GaugeSU2.cc)
exe     := $(patsubst %.o,%,$(mains))
objects := $(filter-out ../nspt.o ../Quenched.o ../Unquenched.o $(mains), $(patsubst %.cc,%.o,$(wildcard *.cc ../*.cc)))

all: $(exe)

$(objects): %.o : %.cc
	$(CC) $(MyO) -c $< -I. -I../ $(PARALLEL_GEN) 

$(mains): %.o : %.cc
	$(CC) $(MyO) -c $< -I. -I../ $(PARALLEL_GEN) 

$(exe): $(objects) $(mains) -lfftw3
	$(CC) $(MyO) -o $@.exe $@.o $(objects) -L$(FFTWdir)lib $(PARALLEL_GEN) 

include $(objects:.o=.d)
include $(mains:.o=.d)

%.d: %.cc 
	@set -e; rm -f $@; \
	$(CC) -MM $(PARALLEL_GEN) $(MyO) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

.PHONY : clean

clean:
	\rm *.exe *.o 
