  Module faddeyeva_v3_mod_rk
    Use set_rk
     USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN,&
    ieee_POSITIVE_INF,ieee_Negative_inf,ieee_POSITIVE_zero,ieee_Negative_zero
    
    Implicit None

    Include 'Faddeyeva_v3_parameters.f90'
    Private
    Public :: faddeyeva_v3_rk, erfc_z, erf_z, erfi_z, erfcx_z, Dawson_z
    
  Contains


  !!::---------
  ELEMENTAL FUNCTION erfc_z( z ,sdgts) RESULT(erfcz)
    ! erfc_z : Complementary error function for a complex
    ! input z in terms of the Faddeyeva function
    ! erfc(z)=exp(-z^2)*w(i*z)

    COMPLEX(rk), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: sdgts
    REAL(rk) ::xx, yy
    COMPLEX(rk) :: erfcz
    yy=AIMAG(z)
    xx=REAL(z)
    IF (yy == zero)THEN
      erfcz=erfc(Real(z))
    ELSEIF(xx == zero .and. (yy*yy)<=log_Rmax)THEN
      !!: 1-erf(iy)=1-(1-exp(y^2)w(y))=1-(1-exp(y^2)*[exp(-y^2)+i*(2/sqrt(pi))*Daw(y)]
      !!: erfc(iy)=1-iexp(y^2)*Imag(w(y))

      !      CALL Faddeyeva_v3_rk(-j1*z, erfcz, sdgts)
      !      erfcz=one-j1*EXP(yy*yy)*AIMAG(erfcz)
      erfcz=one-j1*EXP(yy*yy)*two_sqrt_pi*dawson_rk(yy)

    ELSEIF(xx == zero .and. yy>sqrt_log_Rmax)THEN
      erfcz=one-j1*IEEE_VALUE(1.0_rk, IEEE_Positive_Inf)
    ELSEIF(xx == zero .and. yy<-sqrt_abs_log_rmin)THEN
      erfcz=one-j1*IEEE_VALUE(1.0_rk, IEEE_negative_Inf)
      
     ELSEIF ((xx*xx+yy*yy)<=one)THEN
!      !1- series for erf(z) about z=0
!      !erfz=(2/sqrt(pi))*(z-(1/3)z^3+(1/10)z^5-(1/42)z^7+...+a(n)*z^(2n-1)
!      !where for  & a(n) = (-1)^(n-1)/(2n-1)*(n-1)!  n=1,2,..
      erfcz=one-erf_z(z,sdgts)

      
    ELSE
 !!     
   !  IF (xx >= zero .and. (yy*yy-xx*xx)<=log_Rmax .and. (yy*yy+xx*xx)>127.0_rk)THEN
    !   elseIF (xx >= zero .and. (yy*yy-xx*xx)<=log_Rmax )THEN
       
       IF (xx >= zero .and. (yy*yy-xx*xx)<=log_Rmax )THEN
       
      CALL Faddeyeva_v3_rk(j1*z, erfcz, sdgts)
     
        erfcz=exp(-z*z)*erfcz
        
       
        elseIF ( xx<zero .and. (yy*yy-xx*xx)<=log_Rmax )THEN
       CALL Faddeyeva_v3_rk(-j1*z, erfcz, sdgts)
      
        erfcz=two-EXP(-z*z)*erfcz
     
      ELSEIF((yy*yy-xx*xx)<-abs(log_Rmin) )THEN
        erfcz=zero+j1*zero
      ELSE
        erfcz=EXP(-j1*two*xx*yy)*erfcz*IEEE_VALUE(1.0_rk, IEEE_positive_Inf)
      ENDIF
    ENDIF
  END FUNCTION erfc_z
  !!::--------------------------------------

 !!::---------
  ELEMENTAL FUNCTION erf_z( z ,sdgts) RESULT(erfz)
    ! erf_z : Error function for a complex input z
    ! in terms of the Faddeyeva function
    ! erf(z)= 1- exp(-z^2)w(i*z)

    COMPLEX(rk), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: sdgts
    REAL(rk):: xx, yy, eps_trc
    COMPLEX(rk) :: zsqr, erfz, f_n
    Integer :: n
     Real (rk), Parameter :: cff(70) = (/(-(two*n-1)/(n*(two*n+1)),n=1,70) /)
    
    xx=REAL(z, kind=rk)
    yy=AIMAG(z)
   
    IF (yy == zero .and. xx .ne.zero )THEN
      erfz=erf(xx)
    ELSEIF(xx == zero .and. yy*yy<=log_Rmin)THEN

      erfz=j1*EXP(yy*yy)*two_sqrt_pi*dawson_rk(yy)

    ELSEIF(xx == zero .and. yy>sqrt_log_Rmax)THEN
      erfz=j1*IEEE_VALUE(1.0_rk, IEEE_Positive_Inf)
    ELSEIF(xx == zero .and. yy<-sqrt_abs_log_rmin) THEN
      erfz=j1*IEEE_VALUE(1.0_rk, IEEE_negative_Inf)
    ELSEIF ((xx*xx+yy*yy)<=one)THEN
      !series about z=0
      !erfz=(2/sqrt(pi))*(z-(1/3)z^3+(1/10)z^5-(1/42)z^7+...+a(n)*z^(2n-1)
      !where for  & a(n) = (-1)^(n-1)/(2n-1)*(n-1)!  n=1,2,..

     eps_trc=tolerance(sdgts)
      zsqr=z*z
      f_n = (one, zero)
      erfz = f_n
      Do n = 1, 70
        f_n = f_n*cff(n)*zsqr
        If (abs(real(f_n,kind=rk))< abs(real(erfz,kind=rk))*eps_trc .And. &
        abs(aimag(f_n))<abs(aimag(erfz))*eps_trc) Exit
      erfz = erfz + f_n
      End Do
      erfz=two_sqrt_pi*z*erfz


    ELSEIF (xx>=zero .and. (yy*yy-xx*xx)<=log_Rmax )THEN
       CALL Faddeyeva_v3_rk(j1*z,erfz,sdgts)
        erfz=one-EXP(-z*z)*erfz
    
    ElseIF (xx<zero .and. (yy*yy-xx*xx)<=log_Rmax )THEN
       CALL Faddeyeva_v3_rk(-j1*z,erfz,sdgts)
        erfz=EXP(-z*z)*erfz-one
      ELSEIF( (yy*yy-xx*xx)<-log_Rmin )THEN
        erfz=one-(zero+j1*zero)
      ELSE
        !      erfz=(one+j1)*IEEE_VALUE(1.0_rk, IEEE_QUIET_NAN)
        erfz=one-EXP(-j1*two*xx*yy)*erfz*IEEE_VALUE(1.0_rk, IEEE_positive_Inf)
      

    ENDIF
  END FUNCTION erf_z

!:-
!:-

!!::---------
  ELEMENTAL FUNCTION erfi_z(z, sdgts) RESULT(erfiz)
    ! erfi_z : Imaginary error function for a complex input z
    ! in terms of the Faddeyeva function
    ! erfi(z)= i*(1- exp(z^2)w(z))

    COMPLEX(rk), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: sdgts
    COMPLEX(rk):: erfiz, zsqr, f_n
    REAL(rk):: xx, yy, eps_trc
    integer :: n
     Real (rk), Parameter :: cff(70) = (/((two*n-1)/(n*(two*n+1)),n=1,70) /)
    
    zsqr=z*z
    xx=REAL(z)
    yy=AIMAG(z)
    
      !series about z=0
      !erfiz=(1/sqrt(pi))*(2z+(2/3)z^3+(1/5)z^5+...+a(n)*z^(2n-1)
      !where for n>2 & a(n) = 2/(2n-1)*(n-1)!

     

    IF(xx==zero )THEN
      erfiz=j1*ERF(yy)

    ELSEIF(yy==zero .and. (xx*xx)<log_Rmax)THEN
      erfiz=EXP(xx*xx)*two_sqrt_pi*dawson_rk(xx)


   ELSEIF ((xx*xx+yy*yy)<=one)THEN
      !series about z=0
      !erfiz=(1/sqrt(pi))*(2z+(2/3)z^3+(1/5)z^5+...+a(n)*z^(2n-1)
      !where for n>2 & a(n) = 2/(2n-1)*(n-1)!

     eps_trc=tolerance(sdgts)
      f_n = (one, zero)
      erfiz = f_n
      Do n = 1, 70
        f_n = f_n*cff(n)*zsqr
        If (abs(real(f_n,kind=rk))< abs(real(erfiz,kind=rk))*eps_trc .And. &
        abs(aimag(f_n))<abs(aimag(erfiz))*eps_trc) Exit
      erfiz = erfiz + f_n
      End Do
      erfiz=two_sqrt_pi*z*erfiz
   
   
    ELSE
      
      IF ((xx*xx-yy*yy)<=log_Rmax )THEN
        erfiz=-j1*erf_z(j1*z,sdgts)
      ELSEIF ((xx*xx-yy*yy)<-log_Rmin )THEN
        erfiz=j1*(one-(zero+j1*zero))
      ELSE
        erfiz=j1*(one-erfiz*EXP(j1*two*xx*yy)*IEEE_VALUE(1.0_rk, IEEE_positive_inf))
      ENDIF

    ENDIF
  END FUNCTION erfi_z
  !!::--------------------------------------


 !!::---------
  ELEMENTAL FUNCTION erfcx_z( z , sdgts) RESULT(erfcxz)
    ! erfcx_z : Scaled complementary error function for complex
    ! input z in terms of the Faddeyeva function
    ! erfcxz(z) =  w(iz)

    COMPLEX(rk), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: sdgts
    COMPLEX(rk) ::  erfcxz
    Real(rk) :: xx, yy
    xx=real(z, kind=rk)
    yy=aimag(z)
    IF (yy==zero)THEN
      erfcxz=ERFCx_rk(xx)
    ELSE
      CALL Faddeyeva_v3_rk(j1*z,erfcxz,sdgts)
      erfcxz=erfcxz
    ENDIF
  END FUNCTION erfcx_z
  !!::---------------------------------------


 !!::---------
  ELEMENTAL FUNCTION Dawson_z( z,sdgts) RESULT(dwz)
    ! Dawson_z : Dawson function for a complex input z
    ! in terms of the Faddeyeva function
    ! dwz(z)= i*sqrt(pi)/2 * (exp(-z^2)-w(z))

    COMPLEX(rk), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: sdgts
    COMPLEX(rk) :: dwz
    REAL(rk) :: xx, yy, y_sqr, dw_1, dw_2, dw_3, daw
    xx=REAL(z, kind=rk)
    yy=AIMAG(z)
    y_sqr=yy*yy
    IF (yy == zero)THEN
      dwz=dawson_rk(xx)+j1*zero

!    ELSEIF ((xx*xx+y_sqr)<28.5_rk .and. (y_sqr)<=1.0E-10_rk ) THEN
!      daw=dawson_rk(xx)
!      dw_1=(one-two*xx*daw)
!      dw_2=-(xx*dw_1+daw)
!      dw_3=-0.666666666666666666666666666666666666667_rk*(xx*dw_2+dw_1)
!      dwz=daw+dw_1*j1*yy-dw_2*y_sqr-j1*dw_3*y_sqr*yy

!    ELSE
     elseif((y_sqr + xx*xx)<=log_Rmax )THEN
     dwz= 0.886226925452758013649083741670572591398754234_rk*EXP(-z*z)*erfi_z(z,sdgts) 
!      IF ( xx>=zero)then ! .and. (y_sqr - xx*xx)<=log_Rmax )THEN
!      CALL Faddeyeva_v3_rk(-z,dwz,sdgts)
!        dwz=j1*0.886226925452758013649083741670572591398754234_rk*(dwz-EXP(-z*z))
!      elseIF ( xx<zero)then ! .and. (y_sqr - xx*xx)<=log_Rmax )THEN
!       CALL Faddeyeva_v3_rk(z,dwz,sdgts)
!        dwz=j1*0.886226925452758013649083741670572591398754234_rk*(EXP(-z*z)-dwz)

      ELSEIF((y_sqr - xx*xx)<-log_Rmin )THEN
        dwz=j1*0.886226925452758013649083741670572591398754234_rk*(zero-dwz)
      ELSEIF( (y_sqr - xx*xx)>log_Rmax )THEN
        !       dwz=(one+j1)*IEEE_VALUE(1.0_rk, IEEE_QUIET_NAN)
        dwz=j1*0.886226925452758013649083741670572591398754234_rk*(EXP(-j1*two*xx*yy)-two*EXP(j1*two*xx*ABS(yy)))*&
          IEEE_VALUE(1.0_rk, IEEE_positive_inf)
      ELSE
        dwz=j1*0.886226925452758013649083741670572591398754234_rk*(EXP(-j1*two*xx*yy))*&
          IEEE_VALUE(1.0_rk, IEEE_positive_inf)

      ENDIF
    
  END FUNCTION Dawson_z
  !!::--------------------------------------

!:--

 !!::---------
  ELEMENTAL FUNCTION zeta_z( z ,sdgts) RESULT(zetaz)
    ! zeta_z : Zeta function for a complex input z
    ! in terms of the Faddeyeva function
    ! zetaz(z)= i*sqrt(pi)*w(z)

    COMPLEX(rk), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: sdgts
    COMPLEX(rk) :: zetaz
  
    CALL Faddeyeva_v3_rk(j1*z,zetaz,sdgts)
    zetaz=j1*1.772453850905516027298167483341145182798E+000_rk*zetaz
  END FUNCTION zeta_z
  !!::---------------------------------------


  !!::---------
  ELEMENTAL FUNCTION Voigt_real( xx ,yy,sdgts) RESULT(Voigt_v)
    ! Voigt_v: Voigt Profile from Faddeyeva Function
    ! x=sqrt(ln2)*(nu-nu_0)/gamma_d & y=sqrt(ln2)gamma_l/gamma_d
    ! nu_0   : wave length @ line center
    ! nu     : wave length @ point of calculation
    ! gamma_d: Doppler half-width
    ! gamma_l: Lorentz half-width

    REAL(rk), INTENT(IN) :: xx, yy
    INTEGER, INTENT(IN) :: sdgts
    COMPLEX(rk) :: z, w
    REAL (rk)   :: Voigt_v
    z=(xx+j1*yy)
    CALL Faddeyeva_v3_rk(z,w,sdgts)
    Voigt_v=REAL(w)

  END FUNCTION Voigt_real
  !!::---------------------------------------


  !!::---------
  ELEMENTAL SUBROUTINE Voigt( xx , yy, sdgts, Voigt_v, Voigt_l, dvdx, dvdy)
    ! Voigt_v: Real Voigt function from Faddeyeva function
    ! Voigt_l: Imaginary Voigt function from Faddeyeva function
    ! dvdx   : Partial derivative of Voigt_V w.r.t. x      =dLdy
    ! dvdy   : Partial derivative of Voigt_V w.r.t. y      =-dLdx
    ! x=sqrt(ln2)*(nu-nu_0)/gamma_d & y=sqrt(ln2)gamma_l/gamma_d
    ! nu_0   : wave length @ line center
    ! nu     : wave length @ point of calculation
    ! gamma_d: Doppler half-width
    ! gamma_l: Lorentz half-width

    REAL(rk), INTENT(IN) :: xx, yy
    INTEGER, INTENT(IN) :: sdgts
    REAL(rk), INTENT (OUT):: Voigt_v
    REAL(rk), INTENT(OUT), OPTIONAL:: Voigt_l, dvdx, dvdy
    COMPLEX(rk) :: z, w
    z=(xx+j1*yy)
    CALL Faddeyeva_v3_rk(z, w, sdgts)
    Voigt_v=REAL(w)
    IF (PRESENT(Voigt_l)) Voigt_l=AIMAG(w)

    ! Partial derivative of V w.r.t. x       =dLdy
    IF ( PRESENT(dvdx) ) dvdx=-two*REAL(z*w)

    ! Partial derivative of V w.r.t. y       =-dLdx
    IF ( PRESENT(dvdy) ) dvdy=two*AIMAG(z*w)-two_sqrt_pi

  END SUBROUTINE Voigt
  !!::---------------------------------------















    ! !::************************************************************************************
    ! !::************************************************************************************
    ! ..
    Elemental Subroutine faddeyeva_v3_rk(z, w, sdgts, stat)

      ! ----------
      ! FADDEYEVA_V3_RK IS AN ELEMENTAL FORTRAN SUBROUTINE THAT RECEIVES, AS INPUT,
      ! A SCALAR OR AN ARRAY Z OF COMPLEX NUMBERS AND RETURNS AS OUTPUT A SCALER OR
      ! AN ARRAY OF THE CORRESPONDING FADDEYEVA FUNCTION,W, DEFINED AS
      ! W(Z)=EXP(-Z^2)*ERFC(-I*Z) WHERE ERFC(Z) IS THE COMPLEX COMPLEMENTARY
      ! ERROR FUNCTION. AN OPTIONAL INTEGER INPUT (SDGTS), MAY BE USED REPRESENTING
      ! THE DESIRED NUMBER OF SIGNIFICANT DIGITS IN THE CALCULATED FADDEYEVA FUNCTION.
      ! ------------
      ! ------------
      ! THE ROUTINE CAN BE RUN IN DEFAULT REAL (SINGLE PRECISION), IN DOUBLE
      ! PRECISION OR IN QUADRUPLE PRECISION DEPENDING ON THE CHOICE OF THE INTEGER
      ! PARAMETER "RK" IN THE SUBSIDIARY MODULE "SET_RK"
      ! -------------
      ! THE SUBROUTINE IS SET TO REPRODUCE THE HIGHEST POSSIBLE ACCURACY OBTAINABLE
      ! FROM ALGORITHM 916 [ZAGHLOUL AND ALI,TOMS, VOL. 38, NO. 2, ARTICLE 15:1-22
      ! (2011) & ZAGHLOUL M. ,ACM TOMS, VOL. 42, NO. 3, ARTICLE 26:1-9 (2016) ]
      ! THROUGH REQUESTING A NUMBER OF SIGNIFICANT DIGITS (SDGTS) EQUAL
      ! TO 14 WHEN RUN IN DOUBLE PRECISION. HOWEVER, THIS VERSION (WHEN RUN USING
      ! DOUBLE PRECISION ARITHMETIC EMBODIES SIGNIFICANT EFFICIENCY IMPROVEMENTS
      ! COMPARED TO THE PREVIOUS FORTRAN VERSION "FADDEYEVA_V2.F90" [ZAGHLOUL (2016)]

      ! THE DESIRED NUMBER OF SIGNIFICANT DIGITS CAN BE REDUCED FOR MARGINAL
      ! IMPROVEMENTS OF THE EFFICIENCY AT THE EXPENSE OF ACCURACY.
      ! THE RECOMMENDED RANGE FOR "SDGTS" IS BETWEEN 4 AND 6 FOR SINGLE PRECISION,
      ! 4 AND 14 FOR DOUBLE PRECISION AND BETWEEN 4 AND 31 FOR QUAD PRECISION.
      ! A NUMBER OF SIGNIFICANT DIGITS SMALLER THAN 4 IS NOT RECOMMENDED FOR ACCURACY
      ! CONCERNS, PARTICULARLY REGARDING THE COMPUTATIONS OF THE DERIVATIVES (IF NEEDED),
      ! ANT IT WILL BE AUTOMATICALLY CHANGED TO 4.

      ! VALUES OF SDGTS >14 FOR DOUBLE PRECISION ARE NOT RECOMMENDED FOR PERFORMANCE
      ! CONCERNS AND WILL BE AUTOMATICALLY CHANGED TO 14.
      ! SIMILARLY, VALUES OF SDGTS > 6 FOR SINGLE PRECISION ARE NOT RECOMMENDED
      ! FOR PERFORMANCE CONCERNS AND WILL BE AUTOMATICALLY CHANGED TO 6.
      ! FOR TARGETED ACCURACY BETWEEN 15 AND 31 SIGNFICANT DIGITS, ONE HAS TO
      ! USE QUAD PRECISIN THRROUGH SETTING RK=QP IN THE AUXILARY MODULE "SET_RK"

      ! ALSO WHEN USING QUAD PRECISION VALUES OF SDGTS > 31 ARE CHANGED AUTOMATICALLY
      ! TO 31

      ! AN OPTIONAL "STAT" ARGUMENT IS USED BECAUSE "WRITE" STATEMENTS ARE NOT PERMITTED
      ! IN ELEMENTAL SUBROUTINES. THE "STAT" RETURNS INTEGER VALUES 0,1,2 AND -1
      ! CORRESPONDING TO NORMAL, TOO FEW SIGNIFICANT FIGURES <4, TOO MANY SIGNIFICANT
      ! FIGURES (>31 FOR QUAD PRECISION, >14 FOR DOUBLE PRECISION AND >6 FOR DEFAULT REAL
      ! OR SINGLE PRECISION) AND OVERFLOW, RESPECTIVELY.
      ! ----------
      ! THE ACCOMPANYING DRIVER CODE “FADDEYEVA_V3_DRIVER_RK.F90” CAN BE RUN FOR
      ! COMPUTATION OF THE FADDEYEVA FUNCTION W(Z)=V+IL, ON A SCALAR OR AN ARRAY
      ! OF THE COMPLEX VARIABLE Z.
      ! ----------
      ! AUTHOR: MOFREH R. ZAGHLOUL
      ! UNITED ARAB EMIRATES UNIVERSITY, DECEMBER 3, 2023
      ! ----------

      Complex (rk), Intent (In) :: z
      Integer, Intent (In), Optional :: sdgts
      Complex (rk), Intent (Out) :: w
      Integer, Intent (Out), Optional :: stat
      ! 0 => NORMAL,
      ! 1 => TOO FEW SDGTS (<4)
      ! 2 => TOO MANY SDGTS(>6) FOR SP &
      ! (>14) FOR DP,
      ! (>31) FOR QP
      ! -1 => OVERFLOW
      Real (rk) :: xsqr_plus_ysqr, ysqr_minus_xsqr, x_sqr, y_sqr, xx, yy, x_abs, y_abs

      Integer :: ndgts

      If (.Not. present(sdgts)) ndgts = 4
      If (present(sdgts)) ndgts = sdgts
      If (present(stat)) stat = 0
      If (ndgts<4) Then
        If (present(stat)) stat = 1
        ndgts = 4
      Else If (ndgts>ndgts_max) Then
        If (present(stat)) stat = 2
        ndgts = ndgts_max
      End If

      xx = real(z, kind=rk)
      yy = aimag(z)

      x_abs = abs(xx)
      y_abs = abs(yy)
      If (x_abs>=sqrt_rmin) Then
        x_sqr = xx*xx
      Else
        x_sqr = zero
      End If
      If (y_abs>=sqrt_rmin) Then
        y_sqr = yy*yy
      Else
        y_sqr = zero
      End If

      ysqr_minus_xsqr = y_sqr - x_sqr
      xsqr_plus_ysqr = x_sqr + y_sqr


      If (yy<zero .And. (ysqr_minus_xsqr)>=abs_log_rmin) Then
        If (present(stat)) stat = -1
        Return
      Else If (y_abs<rmin .And. x_abs<sqrt_abs_log_rmin) Then
        w = exp(-x_sqr) + j1*two_sqrt_pi*dawson_rk(xx)
      Else If (y_abs<rmin .And. x_abs>sqrt_abs_log_rmin) Then
        w = zero + j1*two_sqrt_pi*dawson_rk(xx)
        Return
      Else If (x_abs<rmin) Then
        w = erfcx_rk(yy)
        Return
      End If


      If (ndgts>14) Then

        ! -------OUTER REGION ----> ASYMPTOTIC SERIES EXPANSION FOR BIG |Z|
        If (xsqr_plus_ysqr>=127.0E0_rk) Then      ! The asymptotic series needs revision
          w = use_bigz_series(xx+j1*y_abs, tolerance(sdgts+1))
         !return

          ! -----TAYLOR SERIES FOR SMALL |Z| <= 1.0
        Else If (xsqr_plus_ysqr<=1.0E0_rk) Then
          w = use_smallz_series(z, tolerance(sdgts))
          Return

          ! -----VERY SMALL Y REGION ----> TAYLOR SERIES OF DAWSON'S INTEGRAL
        Else If (y_sqr<=1.0E-15_rk) Then
          w = use_taylor_daw(xx, y_abs, 3)

        Else If (y_sqr<=1.0E-6_rk) Then
          w = use_taylor_daw(xx, y_abs, ntaylor2_qp(ndgts-3))

        Else If (y_sqr<=1.0E-3_rk) Then
          w = use_taylor_daw(xx, y_abs, ntaylor3_qp(ndgts-3))

          ! ------ REMAINING REGION ---> LOOP FOR|X|<SQRT(-LOG(R_MIN)) FROM ALGORITHM 916-V2
        Else
          w = use_cycles_1(xx, yy, ncycles_qp(ndgts-3))
        End If


      Else If (ndgts>4) Then


        ! --------OUTER REGION ----> ASYMPTOTIC SERIES EXPANSION FOR LARGE |Z|
        If (xsqr_plus_ysqr>=127.0E0_rk-40.0E0_rk*(sp/rk)) Then

          w = use_bigz_series(xx+j1*y_abs, tolerance(sdgts))
         ! return
          
          ! --------INNER REGION ----> TAYLOR SERIES FOR SMALL |Z|
        Else If (xsqr_plus_ysqr<=2.24E0_rk) Then
          w = use_smallz_series(z, tolerance(sdgts))
          Return

          ! ------VERY SMALL Y REGION ----> TAYLOR SERIES OF DAWSON'S INTEGRAL
        Else If (y_sqr<=1.0E-15_rk) Then
          w = use_taylor_daw(xx, y_abs, 1)

        Else If (y_sqr<=1.0E-6_rk) Then
          w = use_taylor_daw(xx, y_abs, ntaylor2_qp(ndgts-3))

        Else If (y_sqr<=1.0E-3_rk) Then
          w = use_taylor_daw(xx, y_abs, ntaylor3_qp(ndgts-3))

          ! ------ REMAINING REGION ----> LOOP FOR|X|<SQRT(-LOG(R_MIN)) SIMILAR TO
          ! ALGORITHM 916-V2 WITH PRECISION DEPENDENT VALUE OF THE PARAMETER A
        Else
          w = use_cycles_1(xx, yy, ncycles_qp(ndgts-3))
        End If

      Else
        ! HIGHLY EFFICIENT VERSION FOR 4 SD ACCURACY BASED ON ALGORITHM 985
        ! AND CORRECTED HUMLICEK W4 ALGORITHM

        If (xsqr_plus_ysqr>=16.0E3_rk) Then
          ! :-  W=((J1* ONE_SQRT_PI)/Z)
          w = (yy+j1*xx)*(one_sqrt_pi/xsqr_plus_ysqr)

        Else If (xsqr_plus_ysqr>=16.0E1_rk) Then
          w = ((y_abs*(half+xsqr_plus_ysqr))+j1*(xx*(xsqr_plus_ysqr-half)))*(one_sqrt_pi/(( &
            xsqr_plus_ysqr*xsqr_plus_ysqr+ysqr_minus_xsqr)+quarter))

        Else If (xsqr_plus_ysqr>=25.0E0_rk .And. y_sqr>1.0E-13_rk) Then
          w = -ysqr_minus_xsqr + j1*two*xx*y_abs
          w = (one_sqrt_pi*(-y_abs+j1*xx))*(w-2.5_rk)/(w*(w-3.0E0_rk)+0.75E0_rk)

        Else If (yy>=0.195E0_rk*abs(xx)-0.176E0_rk) Then
          w = y_abs - j1*xx
          w = (16.4955_rk+w*(20.20933_rk+w*(11.96482_rk+w*(3.778987_rk+0.5642236_rk*w))))/( &
            16.4955_rk+w*(38.82363_rk+w*(39.27121_rk+w*(21.69274_rk+w*(6.699398_rk+w)))))

        Else If (yy*xx<1.0E-5_rk .And. ysqr_minus_xsqr>=-sqrt_abs_log_rmin) Then
          w = ysqr_minus_xsqr - j1*two*xx*y_abs
          w = exp(ysqr_minus_xsqr) - (-j1*(xx+j1*y_abs)*(36183.31_rk-w*(3321.99_rk-w*(1540.787_rk-w* &
            (219.031_rk-w*(35.7668_rk-w*(1.320522_rk-w*0.56419_rk))))))/(32066.6_rk-w*(24322.84_rk- &
            w*(9022.228_rk-w*(2186.181_rk-w*(364.2191_rk-w*(61.57037_rk-w*(1.841439_rk-w))))))))

        Else
          w = ysqr_minus_xsqr - j1*two*xx*y_abs
          w = exp(w) - (-j1*(xx+j1*y_abs)*(36183.31_rk-w*(3321.99_rk-w*(1540.787_rk-w*(219.031_rk-w* &
            (35.7668_rk-w*(1.320522_rk-w*0.56419_rk))))))/(32066.6_rk-w*(24322.84_rk- &
            w*(9022.228_rk-w*(2186.181_rk-w*(364.2191_rk-w*(61.57037_rk-w*(1.841439_rk-w))))))))

        End If
      End If

      If (yy<zero) Then
        w = two*exp(ysqr_minus_xsqr+j1*two*y_abs*xx) - (real(w, kind=rk)-j1*aimag(w))
      End If

    End Subroutine
    ! !*************************************


    ! :- EVALUATION OF THE ASYMPTOTIC SERIES |Z|^2>127 USING
    ! DYNAMIC CHOICE OF THE NUMBER OF TERMS NEEDED.

    Elemental Function use_bigz_series(z, eps_trc) Result (w)
      Complex (rk), Intent (In) :: z
      Real (rk), Intent (In) :: eps_trc
      Complex (rk) :: alpha, w, f_n
      Integer :: n

      alpha = half/(z*z)
      f_n = (one, zero)
      w = f_n
      Do n = 1, 74, 2
        f_n = n*f_n*alpha
        If (abs(real(f_n,kind=rk))<abs(real(w,kind=rk))*eps_trc .And. abs(aimag(f_n))<abs(aimag(w))* &
          eps_trc) Exit
          w = w + f_n
      End Do
      w = j1*w*one_sqrt_pi/z
    End Function
    ! :--

    ! :- EVALUATION OF THE TAYLOR SERIES FOR SMALL |Z|
    ! (ABRAMOWITZ & STEGUN 7.1.8)USING
    ! DYNAMIC CHOICE OF THE NUMBER OF TERMS NEEDED.

    Elemental Function use_smallz_series(z, eps_trc) Result (w)
      Complex (rk), Intent (In) :: z
      Real (rk), Intent (In) :: eps_trc
      Real (rk) :: zabs
      Complex (rk) :: alpha, w, f_n
      Integer :: n

      alpha = j1*z
      f_n = (one, zero)
      w = f_n
      Do n = 1, 70
        f_n = f_n*gamma_ratio(n)*alpha
        w = w + f_n
        If (abs(real(f_n,kind=rk))<abs(real(w,kind=rk))*eps_trc .And. abs(aimag(f_n))<abs(aimag(w))* &
          eps_trc) Exit
      End Do

    End Function
    ! :--
    ! :*************************


    ! :--
    ! USE  OF A FEW TERMS OF TAYLOR SERIES OF DAWSON'S FUNCTION NEAR
    ! THE REAL AXIS
    Elemental Function use_taylor_daw(xx, yy, nterms) Result (w)
      Real (rk), Intent (In) :: xx, yy
      Integer, Intent (In) :: nterms
      Integer :: n
      Complex (rk) :: w, del_w
      Real (rk) :: x_sqr, y_sqr, ysqr_minus_xsqr, two_xy, dw, x_abs
      Real (rk), Dimension (20) :: d

      x_sqr = xx*xx
      y_sqr = yy*yy

      ysqr_minus_xsqr = y_sqr - x_sqr
      two_xy = two*xx*yy
      dw = dawson_rk(xx)
      d(1) = (one-two*xx*dw)
      w = exp(ysqr_minus_xsqr)*(cos(two_xy)-j1*sin(two_xy)) + (two_sqrt_pi*j1)*(dw+d(1)*j1*yy)

      If (nterms>1) Then
        d(2) = -(xx*d(1)+dw)
        del_w = (two_sqrt_pi*j1)*(-y_sqr)
        w = w + del_w*d(2)
        Do n = 3, nterms + 1
          d(n) = -(2.0E0_rk/n)*(xx*d(n-1)+d(n-2))
          del_w = del_w*j1*yy
          w = w + del_w*d(n)
        End Do
      End If

    End Function


    ! !:- REWORKED ALGORITHM 916 TO INCREASE ACCURACY TO 31 SD
    ! WHEN USING QUADRUPLE PRECISION AIRTHMETIC

    Elemental Function use_cycles_1(xx, yy, ncycles) Result (w)
      Real (rk), Intent (In) :: xx, yy
      Integer, Intent (In) :: ncycles
      Complex (rk) :: w
      Real (rk) :: cos_2yx, del2_tmp, del3_tmp, del3_3_tmp, den1, erfcx_y, exp_x_sqr, exp1, exp2, &
        exp3, exp3_den, exp3_3_den, n3, n3_3, sigma1, sigma2_3, sigma4_5, two_x, two_yx, two_a_x, x, &
        y, ysqr, y_sqr_a_sqr, two_nn_ax
      Integer :: n

      y = abs(yy)
      x = abs(xx)

      ysqr = y*y
      exp_x_sqr = exp(-x*x)
      erfcx_y = erfcx_rk(y)
      two_x = two*x
      two_a_x = a*two_x
      two_yx = y*two_x
      cos_2yx = cos(two_yx)
      n3 = real(ceiling(x/a), kind=rk)
      n3_3 = n3 - one
      y_sqr_a_sqr = inv_a_sqr*ysqr

      sigma1 = zero
      sigma2_3 = zero
      sigma4_5 = zero

      exp1 = exp(-two_a_x)
      exp3 = exp(-two_a_sqr*n3_3+two_a_x)
      exp2 = exp_2a_sqr/(exp3*exp3)

      del2_tmp = one
      del3_tmp = exp(-((a_sqr*n3*n3-two_a_x*n3-two_a_sqr*n3)+x*x+two_a_x+a_sqr))
      del3_3_tmp = exp_a_sqr*exp3

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
          If (x>=1.00E-3_rk*cnst1) Then
            sigma4_5 = sigma4_5 + (n3_3+n)*exp3_den - n*del2_tmp*den1
          Else
            two_nn_ax = n*n*two_a_x
            sigma4_5 = sigma4_5 + two*two_nn_ax*den1*(one+sixth*two_nn_ax*two_a_x+huntw*(two_nn_ax* &
              two_a_x)**2+(1.984126984126984126984126984126984126984E-4_rk)*(two_nn_ax*two_a_x)**3+ &
              (2.755731922398589065255731922398589065256E-6_rk)*(two_nn_ax*two_a_x)**4) ! +&
            ! (HUNTW/42.0E0_RK/72.0E0_RK/111.0E0_RK)*(N*N*TWO_A_X*TWO_A_X)**5)
          End If
        End If
      End Do

      If (y<cnst2 .And. two_yx>rmin) Then
        w = (exp_x_sqr*(erfcx_y*cos_2yx+(two_a_pi/y)*sin(two_yx*half)**2)+y*two_a_pi*(-cos_2yx* &
          sigma1+half*sigma2_3)+j1*sign(one,xx)*(sin(two_yx)*(exp_x_sqr*(-erfcx_y+half*two_a_pi/y)+ &
          two_a_pi*y*sigma1)+two_a_pi_half_a*sigma4_5))

      Else If ((y<cnst2) .And. two_yx<=rmin) Then
        w = (exp_x_sqr*(erfcx_y*cos_2yx+(two_a_pi/y)*sin(two_yx*half)**2)+y*two_a_pi*(-cos_2yx* &
          sigma1+half*sigma2_3)+j1*sign(one,xx)*(y*(two_x*exp_x_sqr*(-erfcx_y+half*two_a_pi/y)+ &
          two_x*two_a_pi*y*sigma1)+two_a_pi_half_a*sigma4_5))
      Else
        w = (exp_x_sqr*(erfcx_y*cos_2yx+(two_a_pi/y)*sin(two_yx*half)**2)+y*two_a_pi*(-cos_2yx* &
          sigma1+half*sigma2_3)+j1*sign(one,xx)*two_a_pi_half_a*sigma4_5)
      End If

    End Function

    ! ==============


    Elemental Function dawson_rk(x)
      ! DAWSON_RK EVALUATES THE DAWSON INTEGRAL OF A REAL ARGUMENT,X;
      ! DAWSON(X)=EXP(-X^2)*INT(0 <= T <= X) EXP(T^2 ) DT

      ! THE CALLING SEQUENCE FOR THIS FUNCTION IS
      ! Y=DAWSON_RK(X)
      ! THE FUNCTION CAN BE EVALUATED USING SINGLE, DOUBLE AND QUADRUPLE
      ! PRECISION DEPENDING ON THE CHOICE OF THE INTEGER "RK" IN THE
      ! ASSOCIATED MODULE "SET_RK"
      ! !---------------

      Implicit None
      Real (rk), Intent (In) :: x
      Real (rk) :: ax, dawson_rk
      Real (rk) :: y100, t
      Integer :: j, ycase, jmin
      Include 'cheb_t100_daw_parameters_sdq.f90'

      ax = abs(x)
      If (big_border_12<ax) Then
        dawson_rk = half/x
      Else If (ax<=xsmall) Then
        dawson_rk = x
      Else If (ax<=small_border_1) Then
        dawson_rk = priv_ncont_frac0(1, x)
      Else If (ax<=small_border_2) Then
        dawson_rk = priv_ncont_frac0(2, x)
      Else If (ax<=small_border_3) Then
        dawson_rk = priv_ncont_frac0(3, x)
      Else If (ax<=small_border_4) Then
        dawson_rk = priv_ncont_frac0(4, x)
      Else If (ax<=small_border_5) Then
        dawson_rk = priv_ncont_frac0(5, x)
      Else If (ax<=small_border_6) Then
        dawson_rk = priv_ncont_frac0(6, x)
      Else If (ax<=small_border_7) Then
        dawson_rk = priv_ncont_frac0(7, x)
      Else If (ax<=cheb_ax) Then
        y100 = one80/(ax+one_pt_8)
        ycase = (y100)
        t = two*y100 - (two*ycase+one)
        jmin = ycase*np_pls_1 + 1
        dawson_rk = t*cffs(jmin)
        Do j = jmin + 1, jmin + np_min_1
          dawson_rk = t*(dawson_rk+(cffs(j)))
        End Do
        dawson_rk = (dawson_rk+cffs(jmin+np))*((1-c_ax)*ax+c_ax) & ! C_AX=(RK/QP)+(SP/RK))
          *sign(one, x)

      Else If (ax<=big_border_1) Then
        dawson_rk = priv_ncont_frac(12, x)
      Else If (ax<=big_border_2) Then
        dawson_rk = priv_ncont_frac(11, x)
      Else If (ax<=big_border_3) Then
        dawson_rk = priv_ncont_frac(10, x)
      Else If (ax<=big_border_4) Then
        dawson_rk = priv_ncont_frac(9, x)
      Else If (ax<=big_border_5) Then
        dawson_rk = priv_ncont_frac(8, x)
      Else If (ax<=big_border_6) Then
        dawson_rk = priv_ncont_frac(7, x)
      Else If (ax<=big_border_7) Then
        dawson_rk = priv_ncont_frac(6, x)
      Else If (ax<=big_border_8) Then
        dawson_rk = priv_ncont_frac(5, x)
      Else If (ax<=big_border_9) Then
        dawson_rk = priv_ncont_frac(4, x)
      Else If (ax<=big_border_10) Then
        dawson_rk = priv_ncont_frac(3, x)
      Else If (ax<=big_border_11) Then
        dawson_rk = priv_ncont_frac(2, x)
      Else If (ax<=big_border_12) Then
        dawson_rk = x/(two*x*x-one)
      End If
      Return
    End Function

    ! :---------------------------------------------------
    ! LAPLACE CONTINUED FRACTION FOR LARGE VALUES OF X
    ! :---------------------------------------------------
    Elemental Function priv_ncont_frac(m, x) Result (y)
      Implicit None
      Real (rk), Intent (In) :: x
      Integer, Intent (In) :: m
      Real (rk) :: y
      Integer :: k

      y = real(m, kind=rk)/x
      Do k = m - 1, 1, -1
        y = real(k, kind=rk)/(x-half*y)
      End Do
      y = one/(two*x-y)
      Return
    End Function

    ! ----------------------------------------------------------
    ! CONTINUED FRACTION EXPANSION FOR |X| NEAR THE ORIGIN
    ! ----------------------------------------------------------
    Elemental Function priv_ncont_frac0(m, x) Result (y)

      Implicit None
      Real (rk), Parameter :: lc_0_cff(15) = (/ 0.66666666666666666666666666666666666666667_rk, &
        -0.26666666666666666666666666666666666666667_rk, &
        0.17142857142857142857142857142857142857142_rk, &
        -0.12698412698412698412698412698412698412698_rk, &
        0.10101010101010101010101010101010101010101_rk, &
        -0.08391608391608391608391608391608391608392_rk, &
        0.07179487179487179487179487179487179487179_rk, &
        -0.06274509803921568627450980392156862745098_rk, &
        0.05572755417956656346749226006191950464396_rk, &
        -0.05012531328320802005012531328320802005013_rk, &
        0.04554865424430641821946169772256728778468_rk, &
        -0.04173913043478260869565217391304347826087_rk, &
        0.03851851851851851851851851851851851851852_rk, &
        -0.03575989782886334610472541507024265644955_rk, &
        0.03337041156840934371523915461624026696329_rk /)

      Real (rk), Intent (In) :: x
      Integer, Intent (In) :: m
      Real (rk) :: y
      Real (rk) :: x_sqr
      Integer :: k

      x_sqr = x*x
      y = lc_0_cff(m)*x_sqr
      Do k = m - 1, 1, -1
        y = lc_0_cff(k)*x_sqr/(one+y)
      End Do
      y = x/(one+y)
      Return
    End Function




    ! ===============================

    Elemental Function erfcx_rk(x)

      ! ERFCX_RK IS AN ELEMENTAL FUNCTION THAT EVALUATES THE
      ! SCALED COMPLEMENTARY ERROR FUNCTION OF A REAL ARGUMENT X,
      ! I.E. "EXP(X*X) * ERFC(X)".

      ! THE FUNCTION ERFCX_RK(X) RECEIVES "X" AS AN INPUT AND RETURNS "ERFCX_RK'
      ! FOR THE SCALED COMPLEMENTARY ERROR FUNCTION USING THE PRECISION
      ! ARITHMETIC DETERMINED BY THE INTEGER "RK" SET IN THE SUBSIDIARY MODULE
      ! "SET_RK".
      ! AN OPTIONAL "STAT" ARGUMENT IS USED BECAUSE "WRITE" STATEMENTS ARE NOT PERMITTED
      ! IN ELEMENTAL SUBROUTINES. THE "STAT" RETURNS INTEGER VALUES 0,1,2 AND -1
      ! CORRESPONDING TO NORMAL, TOO FEW SIGNIFICANT FIGURES <4, TOO MANY SIGNIFICANT
      ! FIGURES (>13 FOR DOUBLE PRECISION AND >6 FOR DEFAULT REAL OR SINGLE PRECISION)
      ! AND OVERFLOW, RESPECTIVELY.
      ! !---------------

      Implicit None
      Real (rk), Intent (In) :: x
      Real (rk) :: erfcx_rk, ax, t, dumy
      Integer :: jmin, j, i
      Include 'cheb_t_erfcx_parameters_sdq.f90'

      ax = abs(x)

      If (x<-9.0E0_rk) Then
        t = x*x
        erfcx_rk = two*exp(t)
        Return
      Else

        If (ax<=zeps(7)) Then
          If (ax<=zeps(1)) Then
            erfcx_rk = one - two_by_sqrt_pi*x
            Return
          Else
            t = x*x
            erfcx_rk = zero
            Do j = 2, 7
              If (ax<=zeps(j)) Then
                Do i = j, 2, -1
                  erfcx_rk = t*(erfcx_rk+c_even(i)+c_odd(i)*x)
                End Do
                erfcx_rk = erfcx_rk + c_even(1) + x*c_odd(1)
                Return
              End If
            End Do

          End If

          ! :- CHEBYSHEV SUBINTERVAL POLYNOMIAL APPROXIMATION
        Else If (ax<=cheb_ax) Then
          t = two_p_one/(ax+two_p_one)
          jmin = n_dvs*t
          jmin = jmin*np_pls_1 + 1
          erfcx_rk = zero
          Do j = 1, (np/4)
            erfcx_rk = t*(t*(t*(t*(erfcx_rk+cffs(jmin))+cffs(jmin+1))+cffs(jmin+2))+cffs(jmin+3))
            jmin = jmin + 4
          End Do
          erfcx_rk = erfcx_rk + cffs(jmin)

        Else

          ! :- LAPLACE CONTINUED FRACTION IN A RATIONAL FUNCTION FORM
          If (x_vbig<=ax) Then
            erfcx_rk = one_by_sqrt_pi/ax

          Else If (x_big<=ax) Then
            erfcx_rk = one_by_sqrt_pi/(ax+half/x)

          Else

            t = ax*ax
            erfcx_rk = one
            dumy = one
outer:      Do j = 1, 4
              i = (j*(1+j))/2
              If (xbig_border(j)<=ax) Then
                Do jmin = i, i + j
                  erfcx_rk = t*erfcx_rk + cf_num(jmin)
                  dumy = t*dumy + cf_dnm(jmin)
                End Do
                erfcx_rk = (one_by_sqrt_pi/ax)*erfcx_rk/dumy
                Exit outer
              End If
            End Do outer
          End If
        End If

        If (x<zero) Then
          t = x*x
          erfcx_rk = two*exp(t) - erfcx_rk
        End If
      End If
      Return
    End Function

  End Module


