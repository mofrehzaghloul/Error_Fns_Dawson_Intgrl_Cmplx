  ! MODULE SET_RK PROVIDES THE KIND TYPE PARAMETER NEEDED TO DEFINE
  ! THE PRECISION OF A COMPLETE PACKAGE ALONG WITH VALUES
  ! FOR COMMONLY USED STANDARD PRECISIONS

  Module set_rk
    Use iso_fortran_env, Only: real32, real64, real128
    Implicit None
    Public :: rk, sp, dp, qp

    ! PRECISIONS
    Integer, Parameter :: sp = real32 ! SINGLE PRECISION
    Integer, Parameter :: dp = real64 ! DOUBLE PRECISION
    Integer, Parameter :: qp = real128 ! QUADRUPLE PRECISION

    ! SET THE PRECISION FOR THE WHOLE PACKAGE
    Integer, Parameter :: rk = dp

    ! TO CHANGE THE DEFAULT PACKAGE-PRECISION CHANGE
    ! THE PARAMETER ASSIGNED TO RK ABOVE TO THE REQUIRED
    ! PRECISION AND RECOMPILE THE COMPLETE PACKAGE.

  End Module
