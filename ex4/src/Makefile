.PHONY: all clean
	
FFLAGS=-g -O1
MODULES=random.f90 metadynamics.f90 routines.f90 structure.f90
OBJECTS=$(MODULES:.f90=.o)

all: mdcode sumhills 

mdcode: $(OBJECTS) mdcode.o
	gfortran $(FFLAGS) -o mdcode $(OBJECTS) mdcode.o

sumhills: $(OBJECTS) sumhills.o
	gfortran $(FFLAGS) -o sumhills $(OBJECTS) sumhills.o

%.o: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm *.o *.mod 
