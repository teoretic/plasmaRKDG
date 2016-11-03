# Makefile created by mkmf $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

include mktemplate


.DEFAULT:
	-touch $@
all: plasma3DR
basisRKDG.o: ./basisRKDG.f90 names.o servicefunc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./basisRKDG.f90
boundc.o: ./boundc.f90 names.o servicefunc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./boundc.f90
flux.o: ./flux.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./flux.f90
initials.o: ./initials.f90 names.o servicefunc.o solution_export.o solvers.o time.o basisRKDG.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./initials.f90
mesh3D.o: ./mesh3D.f90 names.o servicefunc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./mesh3D.f90
names.o: ./names.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./names.f90
plasma3D.o: ./plasma3D.f90 names.o mesh3D.o solvers.o time.o initials.o basisRKDG.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./plasma3D.f90
servicefunc.o: ./servicefunc.f90 names.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./servicefunc.f90
solution_export.o: ./solution_export.f90 names.o servicefunc.o basisRKDG.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./solution_export.f90
solvers.o: ./solvers.f90 names.o servicefunc.o boundc.o flux.o basisRKDG.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./solvers.f90
time.o: ./time.f90 names.o servicefunc.o solution_export.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	./time.f90
SRC = ./plasma3D.f90 ./time.f90 ./solution_export.f90 ./servicefunc.f90 ./names.f90 ./flux.f90 ./boundc.f90 ./initials.f90 ./mesh3D.f90 ./solvers.f90 ./basisRKDG.f90
OBJ = plasma3D.o time.o solution_export.o servicefunc.o names.o flux.o boundc.o initials.o mesh3D.o solvers.o basisRKDG.o
clean: neat
	-rm -f .cppdefs $(OBJ) plasma3DR
neat:
	-rm -f $(TMPFILES)
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
plasma3DR: $(OBJ) 
	$(LD) $(OBJ) -o plasma3DR  $(LDFLAGS)
