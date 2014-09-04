CXX=g++
CXXARGS=-Wall -O3 -g
all: bamQC

bamQC: bamQC.o StringRef.o miscFunctions.o
	$(CXX) $(CXXARGS) -o bamQC bamQC.o StringRef.o miscFunctions.o

bamQC.o: bamQC.cpp
	$(CXX) $(CXXARGS) -c -o bamQC.o bamQC.cpp

StringRef.o: StringRef.cpp StringRef.h
	$(CXX) $(CXXARGS) -c -o StringRef.o StringRef.cpp

miscFunctions.o: miscFunctions.cpp miscFunctions.h
	$(CXX) $(CXXARGS) -c -o miscFunctions.o miscFunctions.cpp

clean:
	rm bamQC