# **fortran\_fftw**

This Fortran project illustrates the usage of FFTW to compute derivatives of periodic functions.

-----------------------------------------------------------------------------

## Requirements
* The GNU Make tool https://www.gnu.org/software/make/
* The GNU gfortran compiler https://gcc.gnu.org/wiki/GFortran
* FFTW https://github.com/FFTW/fftw3.git

-----------------------------------------------------------------------------

## To build the project

Type the following command line arguments
```
git clone https://github.com/jlokimlin/fortran_fftw.git

cd fortran_fftw; make all
```
-----------------------------------------------------------------------------

## Result

```
For sample size n =   8

For f(x) = cos(2x) + sin(3x)

First derivative from periodic periodic_array
 Discretization error =    3.1086244689504383E-015

Second derivative from the periodic periodic_array
 Discretization error =    6.2172489379008766E-015

****************************************


For sample size n =  16

For f(x) = cos(2x) + sin(3x)

First derivative from periodic periodic_array
 Discretization error =    9.7699626167013776E-015

Second derivative from the periodic periodic_array
 Discretization error =    2.1316282072803006E-014

****************************************


For sample size n =  32

For f(x) = cos(2x) + sin(3x)

First derivative from periodic periodic_array
 Discretization error =    1.8207657603852567E-014

Second derivative from the periodic periodic_array
 Discretization error =    1.2745360322696797E-013

****************************************


For sample size n =  64

For f(x) = cos(2x) + sin(3x)

First derivative from periodic periodic_array
 Discretization error =    4.4408920985006262E-014

Second derivative from the periodic periodic_array
 Discretization error =    7.0476957603204937E-013

****************************************


For sample size n = 128

For f(x) = cos(2x) + sin(3x)

First derivative from periodic periodic_array
 Discretization error =    8.8817841970012523E-014

Second derivative from the periodic periodic_array
 Discretization error =    3.7623237858497305E-012

****************************************

 
This result was compiled by GCC version 5.3.1
```