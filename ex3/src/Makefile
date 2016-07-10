.PHONY: all clean
	
FFLAGS=-g -O3
MODULES=random.f90 routines.f90 structure.f90
OBJECTS=$(MODULES:.f90=.o)

all: mdcode mdcommit

mdcode: $(OBJECTS) mdcode.o
	gfortran $(FFLAGS) -o mdcode $(OBJECTS) mdcode.o

mdcommit: $(OBJECTS) mdcommit.o
	gfortran $(FFLAGS) -o mdcommit $(OBJECTS) mdcommit.o
	
%.o: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm *.o *.mod 
