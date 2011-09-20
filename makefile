SHELL = /bin/sh
AR=ar rcs
VPATH = $(FFTWdir)lib

.SUFFIXES: 
.SUFFIXES: .cc .o

LIBS += $(shell pkg-config --libs fftw3)

include makefile.local

mains   := $(patsubst %.cc,%.o,Quenched.cc Unquenched.cc)
exe   := $(patsubst %.o,%,$(mains))
objects := $(filter-out $(mains), $(patsubst %.cc,%.o,$(wildcard *.cc)))

all: $(exe)
	$(MAKE) -C renorm
	$(MAKE) -C tools
	$(MAKE) -C SU2

$(objects): %.o : %.cc
	$(CC) $(MyO) -c $< -I. -I$(FFTWdir)include $(PARALLEL_GEN) 

$(mains): %.o : %.cc
	$(CC) $(MyO) -c $< -I. -I$(FFTWdir)include $(PARALLEL_GEN) 

$(exe): $(mains) $(objects)
	$(CC) $(MyO) -L$(FFTWdir)lib $(LIBS) -o $@.exe $@.o $(objects) ./dSFMT/dSFMT.c $(PARALLEL_GEN) 

include $(objects:.o=.d)
include $(mains:.o=.d)

%.d: %.cc 
	@set -e; rm -f $@; \
	$(CC) -MM $(PARALLEL_GEN) $(MyO) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

.PHONY : clean

clean:
	\rm *.exe *.o renorm/*.o renorm/*.exe tools/*.o tools/*.exe SU2/*.o SU2/*.exe
