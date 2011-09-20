SHELL = /bin/sh

.SUFFIXES: 
.SUFFIXES: .cc .o

include ../makefile.local

objects := $(filter-out ../nspt.o ../Quenched.o ../Unquenched.o ../MyTime.o, $(patsubst %.cc,%.o,$(wildcard ../*.cc)))
mains := $(patsubst %.cc,%.o,$(wildcard *.cc))
exe := $(patsubst %.o,%,$(mains))


all: $(exe) 

$(objects): %.o : %.cc
	$(CC) -c $(PARALLEL_WLOOP) $(MyO) -I. -I../ $< -o $@ -I$(FFTWdir)include 

$(mains): %.o : %.cc
	$(CC) -c $(PARALLEL_WLOOP) $(MyO) -I. -I../ $< -o $@ -I$(FFTWdir)include 

$(exe): $(mains) $(objects)
	$(CC) $(PARALLEL_WLOOP) $(MyO) -o $@.exe $@.o $(objects) -I. -I../  -L$(FFTWdir)lib $(LDFLAGS)

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
