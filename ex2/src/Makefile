.PHONY: all clean
	
FFLAGS=-g -O3
MODULES=random.f90 routines.f90 glecn.f90
OBJECTS=$(MODULES:.f90=.o)

all: mdcode

mdcode: $(OBJECTS) mdcode.o
	gfortran $(FFLAGS) -o mdcode $(OBJECTS) mdcode.o
	
%.o: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm *.o *.mod 
