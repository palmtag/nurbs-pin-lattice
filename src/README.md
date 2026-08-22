# NURBS Mesh creator for pincell arrays

This program will generate an MFEM NURBS mesh for an array of pincells,
such as those typically found in a PWR reactor.
Each pincell can be unique.


## Compiling the Code

To build the code, a fortran compiler must be installed.
The command "gfortran" must be part of your path.

```
  make
```

This will create an executable file called "driver.exe".

The program should compile with any standard Fortran compiler, such as ifort.
You will just need to modify the compile options in the Makefile.


