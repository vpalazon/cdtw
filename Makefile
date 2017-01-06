
OPTIMIZE=1

CXX=g++

ifeq ($(OPTIMIZE),1)
OPTIMIZEFLAGS=-O2
else
OPTIMIZEFLAGS=
endif

CXXFLAGS=-Wall -Wno-deprecated $(OPTIMIZEFLAGS)
CFLAGS=-Wall $(OPTIMIZEFLAGS)

ALL=test_cdtw test_dtw

all: $(ALL)

test_cdtw: test_cdtw.o chronometer.o Samples.o Useful.o LocalDistances.o Dtw.o
	$(CXX) $(CXXFLAGS) -o test_cdtw $^

test_dtw: test_dtw.o chronometer.o Samples.o Useful.o LocalDistances.o Dtw.o
	$(CXX) $(CXXFLAGS) -o test_dtw $^

Dtw.o: Dtw.cc Dtw.hh
	$(CXX) $(CXXFLAGS) -c $<

Samples.o: Samples.cc Samples.hh
	$(CXX) $(CXXFLAGS) -c $<

LocalDistances.o: LocalDistances.cc LocalDistances.hh
	$(CXX) $(CXXFLAGS) -c $<

chronometer.o: chronometer.c chronometer.h
	gcc $(CFLAGS) -c $<

clean: 
	rm -rf *.o $(ALL)
