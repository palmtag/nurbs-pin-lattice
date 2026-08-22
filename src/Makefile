# Makefile
# @version CVS $Id: Makefile,v 1.3 2026/08/21 03:57:13 palmtag Exp $
#--------------------------------------------

FC = gfortran
FFLAGS = -O2 -fbacktrace -Wall -fbounds-check

.SUFFIXES :
.SUFFIXES : .f90 .o

EXEC= driver.exe

default: $(EXEC)

driver.exe: Mod_input.o Mod_nurbsmesh.o driver.o pincell8.o meshutil.o
	$(FC) $(FFLAGS) $^ -o $@

.f90.o :
	$(FC) $(FFLAGS) -c $*.f90

clean:
	rm -f core *.o *.mod $(EXEC)

driver.o: Mod_input.o Mod_nurbsmesh.o
pincell8.o: Mod_input.o Mod_nurbsmesh.o
meshutil.o: Mod_input.o Mod_nurbsmesh.o


