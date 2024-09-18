Faddeyeva: A Fortran90 Package for multiple_precision (single, double, quad) and multi-accuracy computation of 
the set of Error functions (FADDEYEVA_V3_RK(z), ERFC_Z(z), ERF_Z(z), ERFI_Z(z), ERFCX_Z(z)) 
and the Dawson Integral DAWSON_Z(z) of a complex argument (z=x+i*y)
 
Installation
========
Save all the accompanying files in a single directory with a name that is compatible with your Fortran compiler.

Compilation
=========
1- Using Makefile
Use the provided makefile to compile the package:
make -f makefile

2- Manual Compilation
Alternatively, you can compile the package manually using the terminal. For the gfortran compiler, use:
> gfortran -o3 set_rk.f90 Faddeyeva_v2_mod_rk.f90 Faddeyeva_v3_mod_rk.f90 Faddeyeva_v3_driver_rk.f90  -o Faddeyeva_v3_driver_rk

For the Intel Fortran compiler (ifort), replace gfortran with ifort.

Running the Driver Code
==================
Run the driver code by typing the following in the terminal:
>faddeyeva_v3_driver_rk
The code will execute automatically and display the results on the screen.

Selecting Precision
----------------------
To select or change the precision, set the integer "rk" in the set_rk.f90 file to the desired value (sp or dp or qp)

File Descriptions
=============
Makefile
makefile: A script for compiling the package.

Source Files (.f90)
---------------
set_rk.f90                                               : An auxiliary module to select the precision by setting the value of the integer rk.
Faddeyeva_v3_parameters.f90              : Contains numerical constants and parameters used in calculating the set of Error
                                                                 functions and Dawson integral of a complex variable.

cheb_t_erfcx_parameters_sdq.f90          : Contains Chebyshev subinterval polynomial coefficients used in the present
                                                                  erfcx or scaled complementary error function of a real variable used in the 
                                                                  present package

cheb_t100_daw_parameters_sdq.f90     : Contains Chebyshev subinterval polynomial coefficients used in the present
                                                                  Daw(x) or Dawson integral of a real variable used in the present package

Faddeyeva_v2_mod_rk.f90                     : A Fortran90 module including a Fortran translation of Algorithm 916 
                                                                  used herein for efficiency comparison  (for single & double precision only)
Faddeyeva_v3_mod_rk.f90                     : A Fortran90 module including implementation of the  present algorithm  
                                                                  for calculating the set of error functions and the Dawson integral of complex 
                                                                  variables   (works for single, double & quad precision)

Faddeyeva_v3_driver_rk.f90                  : An example of a Fortran driver code or main program.
constants.f90                                           : Contains some constants used in both of Faddeyeva_v2_mod_rk &
                                                                  Faddeyeva_v3_mod_rk
Text Files
-----------
Readme.txt                             : The present file describing the contents of the package.
Disclaimer_and_License.txt  : A file containing a disclaimer and the license agreement for using the package.
Faddeyeva_ref_Maple.txt      : Externally generated data using Maple for accuracy checking using double and quad precision.
Faddeyeva_ref_Maplesp.txt   : Externally generated data using Maple for accuracy checking suitable for single precision.

Other externally generated data files for accuracy check
-----------
Reference_data_erfcz_25_25.txt
Reference_data_erfz_25_25.txt
Reference_data_erfiz_25_25.txt
Reference_data_dawz_12_12.txt
Reference_data_erfcxz_25_25.txt

