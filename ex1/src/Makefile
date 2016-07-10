.PHONY: all clean
	
FFLAGS=-g -O3 -Wall
MODULES=random.f90 routines.f90
OBJECTS=$(MODULES:.f90=.o)

all: mcnvt

mcnvt: $(OBJECTS) mcnvt.o
	gfortran $(FFLAGS) -o mcnvt $(OBJECTS) mcnvt.o
	
%.o: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm *.o *.mod 
