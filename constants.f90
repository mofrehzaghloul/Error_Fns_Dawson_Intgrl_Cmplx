  ! Common constants used in Faddeyeva_v2 and Faddeyeva_v3
    Real (rk), Parameter :: zero = 0.0E0_rk
    Real (rk), Parameter :: huntw = 1.0E0_rk/120.0E0_rk
    Real (rk), Parameter :: sixth = 1.0E0_rk/6.0E0_rk
    Real (rk), Parameter :: quarter = 0.25E0_rk
    Real (rk), Parameter :: half = 0.5E0_rk
    Real (rk), Parameter :: one = 1.0E0_rk
    Real (rk), Parameter :: three_halfs = 1.5E0_rk
    Real (rk), Parameter :: two = 2.0E0_rk
    Real (rk), Parameter :: three = 3.0E0_rk
    Real (rk), Parameter :: four = 4.0E0_rk
    Real (rk), Parameter :: five = 5.0E0_rk
    Real (rk), Parameter :: ten = 10.0E0_rk
    Real (rk), Parameter :: rmin = tiny(0.0E0_rk) ! UNDERFLOW LIMIT
    Real (rk), Parameter :: log_rmin = log(rmin)
    Real (rk), Parameter :: sqrt_rmin = sqrt(tiny(0.0E0_rk))
    Real (rk), Parameter :: rmax = huge(1.0E0_rk) ! OVERFLOW LIMIT
    Real (rk), Parameter :: abs_log_rmin = abs(log(rmin))
    Real (rk), Parameter :: log_rmax = log(rmax)
    Real (rk), Parameter :: sqrt_abs_log_rmin = sqrt(abs_log_rmin)
    Real (rk), Parameter :: sqrt_log_rmax = sqrt(log_rmax)

    Real (rk), Parameter :: pi = 3.1415926535897932384626433832795028841971693_rk
    Real (rk), Parameter :: one_sqrt_pi = 0.564189583547756286948079451560772585844_rk ! 1/SQRT(PI)
    Real (rk), Parameter :: two_sqrt_pi = 1.128379167095512573896158903121545171688_rk ! 2/SQRT(PI)
   
    Complex (rk), Parameter :: j1 = (zero, one)
   
    Real (rk), Parameter :: eps = epsilon(1.0E0_rk), epsinv = 1.0E0_rk/eps, &
      xsmall = sqrt(eps*1.5_rk), xlarge = sqrt(0.5E0_rk/eps), xhuge = huge(1.0E0_rk)
      
 Real (rk), Parameter :: one3 = 1.0E0_rk/3.0E0_rk, two_p_one = 2.1_rk, &
      sqrt_lg_xhuge = sqrt(log(xhuge)), zeps0 = eps/one_sqrt_pi/4.0E0_rk, &
      zeps1 = sqrt(eps/2.0E0_rk), zeps2 = eps**0.25_rk, zeps3 = (3.0E0_rk*eps)**(1.0E0_rk/6.0E0_rk)
