  Module faddeyeva_v2_mod_rk
    Use set_rk, Only: rk, sp
  !  Use rk_erfcx_cody, Only: erfc_scaled => rk_erfcx ! Comment for compilers providing Erfc_scaled
    Implicit None

    Include 'constants.f90'

    Real (rk), Parameter :: a = 0.5_rk + (sp/rk)*.2_rk ! A=0.5 FOR R8 & 0.7 FOR SP
    Real (rk), Parameter :: a_pi = a/3.141592653589793238462643383279502884197_rk ! 0.25_RK*A/ATAN(1.0E0_RK)
                                                                                  ! !A/PI
    Real (rk), Parameter :: a_sqr = a**2
    Real (rk), Parameter :: inv_a_sqr = 1.0E0_rk/a_sqr
    Real (rk), Parameter :: two_a_pi = 2.0E0_rk*a_pi
    Real (rk), Parameter :: two_a_pi_half_a = 0.5_rk*two_a_pi*a
    Real (rk), Parameter :: two_a_sqr = 2.0E0_rk*a_sqr

    Integer, Parameter :: ncycles_4 = 2*(rk/sp+1) + 1
    Integer, Parameter :: ncycles_5 = 4*(rk/sp) + (sp/rk)
    Integer, Parameter :: ncycles_6 = ncycles_5 + 1
    Integer, Parameter :: ndigits_max = 13 - 7*(sp/rk)
    Integer, Parameter :: bigz_border(10) = (/ 107, 110, 180, 380, 810, 1750, 3800, 8200, 17500, &
      38000 /)
    Integer, Parameter :: ncycles(10) = (/ ncycles_4, ncycles_5, ncycles_6, 9, 10, 10, 11, 11, 12, &
      13 /)


    ! ..EXP(-TWO_A_SQR*(N-1))
    Real (rk), Parameter :: exp_2a_sqr_n_1(13) = (/ (1.000000000000000E+000_rk*(1+(sp/rk)*0.0E0_rk)) &
      , (6.065306597126334E-001_rk*(1-(sp/rk)*3.812166081938592E-001_rk)), ( &
      3.678794411714423E-001_rk*(1-(sp/rk)*6.171071140248879E-001_rk)), (2.231301601484298E-001_rk*( &
      1-(sp/rk)*7.630722413178782E-001_rk)), (1.353352832366127E-001_rk*(1- &
      (sp/rk)*8.533930378696498E-001_rk)), (8.208499862389880E-002_rk*(1- &
      (sp/rk)*9.092820467105875E-001_rk)), (4.978706836786394E-002_rk*(1- &
      (sp/rk)*9.438652371658662E-001_rk)), (3.019738342231850E-002_rk*(1- &
      (sp/rk)*9.652647410552614E-001_rk)), (1.831563888873418E-002_rk*(1- &
      (sp/rk)*9.785063986549101E-001_rk)), (1.110899653824231E-002_rk*(1- &
      (sp/rk)*9.867001164575563E-001_rk)), (6.737946999085467E-003_rk*(1- &
      (sp/rk)*9.917702529509800E-001_rk)), (4.086771438464067E-003_rk*(1- &
      (sp/rk)*9.949075692073008E-001_rk)), (2.478752176666359E-003_rk*(1- &
      (sp/rk)*9.968488884015556E-001_rk)) /)

    ! REAL(RK), PARAMETER :: INV_ASQR_EXP_ASQR_NSQR(NDIGITS_MAX) = &
    ! INV_A_SQR*(/ ( EXP(-A_SQR * NN**2), NN = 1, NDIGITS_MAX ) /)
    ! ..EXP(-A_SQR * NN**2)/A**2
    Real (rk), Parameter :: inv_asqr_exp_asqr_nsqr(13) = (/ inv_a_sqr*7.788007830714049E-001_rk*(1-( &
      sp/rk)*2.133721389334466E-001_rk), inv_a_sqr*3.678794411714423E-001_rk*(1-(sp/ &
      rk)*6.171071140248879E-001_rk), inv_a_sqr*1.053992245618643E-001_rk*(1-(sp/ &
      rk)*8.846748789619375E-001_rk), inv_a_sqr*1.831563888873418E-002_rk*(1-(sp/ &
      rk)*9.785063986549101E-001_rk), inv_a_sqr*1.930454136227709E-003_rk*(1-(sp/ &
      rk)*9.975212478233336E-001_rk), inv_a_sqr*1.234098040866796E-004_rk*(1-(sp/ &
      rk)*9.998231130977574E-001_rk), inv_a_sqr*4.785117392129010E-006_rk*(1-(sp/ &
      rk)*9.999921891752664E-001_rk), inv_a_sqr*1.125351747192591E-007_rk*(1-(sp/ &
      rk)*9.999997865791929E-001_rk), inv_a_sqr*1.605228055185612E-009_rk*(1-(sp/ &
      rk)*9.999999963915950E-001_rk), inv_a_sqr*1.388794386496402E-011_rk*(1-(sp/ &
      rk)*9.867001164575563E-001_rk), inv_a_sqr*7.287724095819692E-014_rk*(1-(sp/ &
      rk)*9.999999999997556E-001_rk), inv_a_sqr*2.319522830243570E-016_rk*(1-(sp/ &
      rk)*9.999999999999990E-001_rk), inv_a_sqr*4.477732441718302E-019_rk*(1-(sp/ &
      rk)*1.000000000000000E+000_rk) /)


    Private
    Public :: faddeyeva_v2_rk

  Contains

    ! ..
    ! IF DERIVATIVES ARE NOT REQUIRED, ONLY POSITIVE VALUES OF THE IMAGINARY
    ! PART OF THE INPUT COMPLEX NUMBER Z ARE CONSIDERED AND WHEN THE  OPTIONAL
    ! "STAT" ARGUMENT IS NOT NEEDED, EFFICIENCY CAN BE IMPROVED BY
    ! COMMENTING THE FIRST LINE DEFINING THE SUBROUTINE AND UNCOMMENTING
    ! THE FOLLOWING LINE (LINE 64) TOGETHER WITH COMMENTING THE LINES
    ! (128,130,138,140,143, 157-160,231,234)
    ! ..


    ! ELEMENTAL SUBROUTINE FADDEYEVA_V2_RK ( Z, W, SDGTS, DVDX, DVDY, STAT )
    Elemental Subroutine faddeyeva_v2_rk(z, w, sdgts, stat)

      ! ----------
      ! FADDEYEVA_V2_RK IS AN ELEMENTAL FORTRAN SUBROUTINE THAT RECEIVES, AS INPUT, A
      ! SINGLE DUMMY SCALAR COMPLEX NUMBER (Z=X+IY) AND RETURNS AS OUTPUT
      ! THE FADDEYEVA FUNCTION,W, DEFINED AS W(Z)=EXP(-Z^2)*ERFC(-I*Z) WHERE ERFC(Z)
      ! IS THE COMPLEX COMPLEMENTARY ERROR FUNCTION. AN OPTIONAL INTEGER INPUT (SDGTS),
      ! MAY BE USED REPRESENTING THE DESIRED NUMBER OF SIGNIFICANT FIGURES IN
      ! THE CALCULATED FADDEYEVA FUNCTION.
      ! IN ADDITION, THE SUBROUTINE RETURNS AS OPTIONAL THE DERIVATIVES
      ! OF THE REAL PART OF THE FADDEYEVA FUNCTION (DVDX AND DVDY) WITH RESPECT TO
      ! THE REAL AND IMAGINARY PARTS OF Z, RESPECTIVELY.
      ! ------------
      ! THIS VERSION USES A REFORMED FORM OF HUMLICEK'S W4 RATIONAL APPROXIMATION
      ! [ZAGHLOUL, M. R. 2015. A SIMPLE REFORM FOR TREATING THE LOSS OF ACCURACY
      ! OF HUMLICECK’S W4 ALGORITHM NEAR THE REAL AXIS. ARXIV:1505.05596V1 [ASTRO-PH.IM]]
      ! FOR THE CASE OF 4 SIGNIFICANT FIGURES
      ! ------------
      ! THE ROUTINE CAN BE RUN IN DEFAULT REAL (SINGLE PRECISION) OR IN DOUBLE
      ! PRECISION DEPENDING ON THE CHOICE OF THE INTEGER PARAMETER RK IN THE
      ! SUBSIDIARY MODULE SET_RK
      ! -------------
      ! THE SUBROUTINE IS SET TO REPRODUCE THE HIGHEST POSSIBLE ACCURACY OBTAINABLE
      ! FROM ALGORITHM 916 [ZAGHLOUL AND ALI,TOMS, VOL. 38, NO. 2, ARTICLE 15:1-22
      ! (2011)] THROUGH REQUESTING A NUMBER OF SIGNIFICANT FIGURES (SDGTS) EQUAL
      ! TO 13 WHEN RUN IN DOUBLE PRECISION.
      ! THE DESIRED NUMBER OF SIGNIFICANT FIGURES CAN BE REDUCED FOR MARGINAL
      ! IMPROVEMENTS OF THE EFFICIENCY AT THE EXPENSE OF ACCURACY.
      ! THE RECOMMENDED RANGE FOR "SDGTS" IS BETWEEN 4 AND 13 FOR DOUBLE PRECISION
      ! AND BETWEEN 4 AND 6 FOR SINGLE PRECISION.
      ! A NUMBER OF SIGNIFICANT FIGURES SMALLER THAN 4 IS NOT RECOMMENDED FOR ACCURACY
      ! CONCERNS, PARTICULARLY REGARDING THE COMPUTATIONS OF THE DERIVATIVES, ANT IT
      ! WILL BE DIRECTLY CHANGED TO 4.

      ! VALUES OF SDGTS >13 FOR DOUBLE PRECISION ARE NOT RECOMMENDED FOR PERFORMANCE
      ! CONCERNS AND WILL BE AUTOMATICALLY CHANGED TO 13.
      ! SIMILARLY, VALUES OF SDGTS > 6 FOR SINGLE PRECISION ARE NOT RECOMMENDED
      ! FOR PERFORMANCE CONCERNS AND WILL BE AUTOMATICALLY CHANGED TO 6.

      ! AN OPTIONAL "STAT" ARGUMENT IS USED BECAUSE "WRITE" STATEMENTS ARE NOT PERMITTED
      ! IN ELEMENTAL SUBROUTINES. THE "STAT" RETURNS INTEGER VALUES 0,1,2 AND -1
      ! CORRESPONDING TO NORMAL, TOO FEW SIGNIFICANT FIGURES <4, TOO MANY SIGNIFICANT
      ! FIGURES (>13 FOR DOUBLE PRECISION AND >6 FOR DEFAULT REAL OR SINGLE PRECISION)
      ! AND OVERFLOW, RESPECTIVELY.
      ! ----------
      ! THE ACCOMPANYING DRIVER CODE “FADDEYEVA_DRIVER_RK.F90” CAN BE RUN FOR
      ! COMPUTATION OF THE FADDEYEVA FUNCTION W(Z)=V+IL AND THE PARTIAL
      ! DERIVATIVES OF ITS REAL PART, V(X,Y) (OPTIONAL), ON A SCALAR OR AN ARRAY OF
      ! THE COMPLEX VARIABLE Z. THE PARTIAL DERIVATIVES OF THE IMAGINARY PART, L(X,Y),
      ! ARE SIMPLY GIVEN BY EQ. (23) IN THE ORIGINAL ARTICLE OF ALGORITHM 916,
      ! TOMS. VOL. 38, NO. 2, ARTICLE 15:1-22 (2011) AND AS COMMENTED BY THE
      ! END OF THIS SUBROUTINE.
      ! AN EXAMPLE OF GENERATING AN ARRAY OF Z IS INCLUDED IN THE DRIVER CODE.
      ! ----------
      ! AUTHOR: MOFREH R. ZAGHLOUL
      ! UNITED ARAB EMIRATES UNIVERSITY, JUNE 17, 2015
      ! ----------


      ! PRIVATE MATHEMATICAL AND REPEATEDLY USED CONSTANTS

      Complex (rk), Intent (In) :: z
      Integer, Intent (In), Optional :: sdgts
      Complex (rk), Intent (Out) :: w
      ! REAL(RK), INTENT(OUT), OPTIONAL :: DVDX, DVDY
      Real (rk) :: xsqr_plus_ysqr, ysqr_minus_xsqr, x_sqr, y_sqr, s, y_abs
      Integer, Intent (Out), Optional :: stat ! 0 => NORMAL,
      Real (rk) :: xx, yy

      ! 1 => TOO FEW (<4)
      ! 2 => TOO MANY SDGTS (>6) FOR SP &
      ! (>13) FOR R8,
      ! -1 => OVERFLOW
      Integer :: ndgts

      If (.Not. present(sdgts)) ndgts = 4
      If (present(sdgts)) ndgts = sdgts
      If (present(stat)) stat = 0
      If (ndgts<4) Then
        If (present(stat)) stat = 1
        ndgts = 4
      Else If (ndgts>ndigits_max) Then
        If (present(stat)) stat = 2
        ndgts = ndigits_max
      End If

      xx = real(z)
      yy = aimag(z)
      x_sqr = xx*xx
      y_sqr = yy*yy
      
      xsqr_plus_ysqr = (x_sqr+y_sqr)
      ysqr_minus_xsqr = y_sqr - x_sqr


      If (yy<zero .And. (ysqr_minus_xsqr)>=abs_log_rmin) Then
        If (present(stat)) stat = -1
        Return
      Else If (abs(xx)<rmin) Then
        w = erfc_scaled(yy)
        Return
      End If


      If (ndgts>4) Then
        ! -------------
        ! ASYMPTOTIC EXPRESSION FOR LARGE |Z| (|Z|^2>BIGZ_BORDER)
        ! EXAPNDED FORM OF (J1*ONE_SQRT_PI)*TWO*(Z*Z - ONE)/(Z*(TWO*Z*Z - 3.0E0_RK))
        If (xsqr_plus_ysqr>=bigz_border(ndgts-3)) Then
          ! W=(J1*ONE_SQRT_PI)*TWO*(Z*Z - ONE)/(Z*(TWO*Z*Z - 3.0E0_RK))
          y_abs = abs(yy)
          w = one_sqrt_pi*((y_abs*((ysqr_minus_xsqr+one)*(ysqr_minus_xsqr+three_halfs)+ &
            four*x_sqr*y_sqr+x_sqr))+j1*(xx*((ysqr_minus_xsqr+one)*(ysqr_minus_xsqr+three_halfs)+ &
            four*x_sqr*y_sqr-y_sqr)))/(xsqr_plus_ysqr*((ysqr_minus_xsqr+ &
            three_halfs)*(ysqr_minus_xsqr+three_halfs)+four*x_sqr*y_sqr))


          If (yy<zero) Then
            w = two*exp(ysqr_minus_xsqr+j1*two*y_abs*xx) - (real(w)-j1*aimag(w))
          End If

          ! --------------
          ! FOR X=0, USE THE ASYMPTOTIC EXPRESSIONS FOR X--->0 FROM EQ. (6) IN THE ORIGINAL
          ! ARTICLE OF THE ALGORITHM 916, TOMS. VOL. 38, NO. 2, ARTICLE 15:1-22 (2011).
        Else If (abs(xx)<rmin .And. xsqr_plus_ysqr<bigz_border(ndgts-3)) Then
          w = cmplx(erfc_scaled((yy)), 0.0E0_rk, kind=rk)


          ! -------CALCULATING FADDEYEVA FN FOR VALUES OF 0<|X|<SQRT(-LOG(RMIN))
          ! -------WHILE |Z|^2<BIGZ_BORDER
        Else If (abs(xx)>zero .And. abs(xx)<sqrt_abs_log_rmin .And. xsqr_plus_ysqr<bigz_border(ndgts-3)) &
            Then
          w = use_cycles_1(xx, yy, ncycles(ndgts-3))
        Else

          ! -------CALCULATING FADDEYEVA FOR VALUES OF X>=SQRT(-LOG(RMIN))
          ! -------WHILE |Z|^2<BIGZ_BORDER
          w = use_cycles_2(xx, yy, ncycles(ndgts-3))
        End If
      Else
        ! .. REFORMED HUMLICEK ROUTINE
        y_abs = abs(yy)
        s = abs(xx) + y_abs
        If (s>=15.0E0_rk) Then
          w = ((y_abs*(half+xsqr_plus_ysqr))+j1*(xx*(xsqr_plus_ysqr-half)))*(one_sqrt_pi/(( &
            xsqr_plus_ysqr*xsqr_plus_ysqr+ysqr_minus_xsqr)+quarter))

        Else If (s<15.0E0_rk .And. s>=5.5_rk .And. y_sqr>1.000E-12_rk) Then
          w = -ysqr_minus_xsqr + j1*two*xx*y_abs
          w = (one_sqrt_pi*(-y_abs+j1*xx))*(w-2.5_rk)/(w*(w-3.0E0_rk)+0.75_rk)

        Else If (yy>=0.195_rk*abs(xx)-0.176_rk) Then
          w = y_abs - j1*xx
          w = (16.4955_rk+w*(20.20933_rk+w*(11.96482_rk+w*(3.778987_rk+0.5642236_rk*w))))/( &
            16.4955_rk+w*(38.82363_rk+w*(39.27121_rk+w*(21.69274_rk+w*(6.699398_rk+w)))))
        Else
          w = ysqr_minus_xsqr - j1*two*xx*y_abs
          w = exp(w) - (-j1*(xx+j1*y_abs)*(36183.31_rk-w*(3321.99_rk-w*(1540.787_rk-w*(219.031_rk-w* &
            (35.7668_rk-w*(1.320522_rk-w*0.56419_rk))))))/(32066.6_rk-w*(24322.84_rk- &
            w*(9022.228_rk-w*(2186.181_rk-w*(364.2191_rk-w*(61.57037_rk-w*(1.841439_rk-w))))))))
        End If

        If (yy<zero) Then
          w = two*exp(ysqr_minus_xsqr+j1*two*y_abs*xx) - (real(w)-j1*aimag(w))
        End If

      End If

      ! --------------- CALCULATION OF THE DERIVATIVES
      ! PARTIAL DERIVATIVE OF V W.R.T. X       =DLDY
      ! IF ( PRESENT(DVDX) ) DVDX=-TWO*REAL(Z*W)

      ! PARTIAL DERIVATIVE OF V W.R.T. Y       =-DLDX
      ! IF ( PRESENT(DVDY) ) DVDY=TWO*AIMAG(Z*W)-TWO_SQRT_PI

    End Subroutine


    ! ==========================================
    Elemental Function use_cycles_1(xx, yy, ncycles) Result (w)

      Real (rk), Intent (In) :: xx, yy
      Integer, Intent (In) :: ncycles
      Complex (rk) :: w
      Real (rk) :: cos_2yx, del2_tmp, del3_tmp, del3_3_tmp, den1, erfcx_y, exp_x_sqr, exp1, exp2, &
        exp3, exp3_den, exp3_3_den, l_old, n3, n3_3, sigma1, sigma2_3, sigma4_5, two_a_pi_y, &
        two_a_sqr_n3, two_x, two_yx, two_a_x, v_old, x, x_sqr, y, y_sqr_a_sqr

      Integer :: n

      x = abs(xx)
      y = max(rmin, abs(yy))
      erfcx_y = erfc_scaled(y)
      x_sqr = x*x
      two_x = two*x
      two_a_x = a*two_x
      two_yx = y*two_x
      cos_2yx = cos(two_yx)
      exp_x_sqr = exp(-x_sqr)
      two_a_pi_y = two_a_pi/y
      n3 = real(ceiling(x/a), kind=rk)
      n3_3 = n3 - one
      two_a_sqr_n3 = two_a_sqr*n3
      y_sqr_a_sqr = inv_a_sqr*y*y

      sigma1 = zero
      sigma2_3 = zero
      sigma4_5 = zero
      v_old = exp_x_sqr*(erfcx_y*cos_2yx+two_a_pi_y*sin(two_yx*half)**2)
      l_old = exp_x_sqr*(-erfcx_y+half*two_a_pi_y)

      exp1 = exp(-two_a_x)
      exp3 = exp(-two_a_sqr_n3+two_a_x+two_a_sqr)
      exp2 = exp(two_a_sqr)/(exp3*exp3)

      del2_tmp = one
      del3_tmp = exp(-((a_sqr*n3*n3-two_a_x*n3-two_a_sqr_n3)+x_sqr+two_a_x+a_sqr))
      del3_3_tmp = exp(-a_sqr)*exp3

      Do n = 1, ncycles
        den1 = inv_asqr_exp_asqr_nsqr(n)*exp_x_sqr/(n*n+y_sqr_a_sqr)
        del2_tmp = del2_tmp*exp1
        del3_tmp = del3_tmp*exp3
        exp3_den = del3_tmp*inv_asqr_exp_asqr_nsqr(n)/((n3_3+n)**2+y_sqr_a_sqr)
        sigma1 = sigma1 + den1
        If (n3_3>=n) Then
          del3_3_tmp = del3_3_tmp*exp2
          exp3_3_den = del3_3_tmp*del3_tmp*inv_asqr_exp_asqr_nsqr(n)/((n3-n)**2+y_sqr_a_sqr)
          sigma2_3 = sigma2_3 + del2_tmp*den1 + exp3_3_den + exp3_den
          sigma4_5 = sigma4_5 + (n3-n)*exp3_3_den + (n3_3+n)*exp3_den - n*del2_tmp*den1
        Else
          sigma2_3 = sigma2_3 + del2_tmp*den1 + exp3_den
          If (x>=5.0E-3_rk) Then
            sigma4_5 = sigma4_5 + (n3_3+n)*exp3_den - n*del2_tmp*den1
          Else
            sigma4_5 = sigma4_5 + two*n*n*two_a_x*den1*(one+sixth*n*n*two_a_x*two_a_x+huntw*(n*n* &
              two_a_x*two_a_x)**2)
          End If
        End If
      End Do

      If ((y<five) .And. two_yx>rmin) Then
        w = (v_old+y*two_a_pi*(-cos_2yx*sigma1+half*sigma2_3)+j1*sign(one,xx)*(sin(two_yx)*(l_old+ &
          two_a_pi*y*sigma1)+two_a_pi_half_a*sigma4_5))
      Else If ((y<five) .And. two_yx<=rmin) Then
        w = (v_old+y*two_a_pi*(-cos_2yx*sigma1+half*sigma2_3)+j1*sign(one,xx)*(y*(two_x*l_old+ &
          two_x*two_a_pi*y*sigma1)+two_a_pi_half_a*sigma4_5))
      Else
        w = (v_old+y*two_a_pi*(-cos_2yx*sigma1+half*sigma2_3)+j1*sign(one,xx)*two_a_pi_half_a* &
          sigma4_5)
      End If
      If (yy<zero) Then
        w = two*exp(y*y-x*x+j1*two*y*xx) - (real(w)-j1*aimag(w))
      End If

    End Function

    ! =========================================
    Elemental Function use_cycles_2(xx, yy, ncycles) Result (w)

      Real (rk), Intent (In) :: xx, yy
      Integer, Intent (In) :: ncycles
      Complex (rk) :: w
      Real (rk) :: del3_3_tmp, den1, exp2, exp3_den, exp3_3_den, factor, n3, n3_3, sigma3, sigma5, &
        x, y, y_sqr_a_sqr
      Integer :: n

      x = abs(xx)
      y = max(rmin, abs(yy))
      n3 = real(ceiling(x/a), kind=rk) ! CEILING OF (X/A)
      n3_3 = n3 - 1
      y_sqr_a_sqr = inv_a_sqr*y*y
      sigma3 = zero
      sigma5 = zero
  
      del3_3_tmp = exp((two*a*x-two_a_sqr*n3)+a_sqr)
      exp2 = 1/(del3_3_tmp*del3_3_tmp)
      factor = del3_3_tmp
      exp3_den = inv_a_sqr*exp(-(a*n3_3-x)**2)

      Do n = 1, ncycles
        del3_3_tmp = del3_3_tmp*exp2
        exp3_den = exp3_den*factor*exp_2a_sqr_n_1(n)
        exp3_3_den = exp3_den*del3_3_tmp/((n3-n)**2+y_sqr_a_sqr)
        den1 = exp3_den/((n3_3+n)**2+y_sqr_a_sqr)
        sigma3 = sigma3 + exp3_3_den + den1
        sigma5 = sigma5 + (n3-n)*exp3_3_den + (n3_3+n)*den1
      End Do

      w = y*a_pi*sigma3 + j1*sign(one, xx)*two_a_pi_half_a*sigma5

      If (yy<zero) Then
        w = two*exp(y*y-x*x+j1*two*y*xx) - (real(w)-j1*aimag(w))
      End If

    End Function
  End Module
