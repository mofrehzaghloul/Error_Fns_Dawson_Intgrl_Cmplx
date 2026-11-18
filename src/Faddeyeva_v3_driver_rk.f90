  ! last change:  Mofreh Zaghloul   April 24, 2023
  ! .. for compilation
  ! >(gfortran -o3 set_rk.f90 Faddeyeva_v2_mod_rk.f90 Faddeyeva_v3_mod_rk.f90
  ! Faddeyeva_v3_driver_rk.f90 -o Faddeyeva_v3_driver_rk)
  ! ! -wall -fcheck=bounds
  ! !! or
  ! >(ifort -o3 set_rk.f90 faddeyeva_v2_mod_rk.f90 faddeyeva_v3_mod_rk.f90
  ! faddeyeva_v3_driver_rk.f90 -o faddeyeva_v3_driver_rk)
  ! !  /fpe:0 /traceback
  ! ! /warn:all
  ! ! /warn:declarations
  ! ! /stand:f90
  ! ! /cu

  Program FADDEYEVA_V3_DRIVER_RK
    ! .. USE STATEMENTS ..
    Use, Intrinsic :: IEEE_ARITHMETIC
    Use, Intrinsic :: ISO_FORTRAN_ENV, Only: COMPILER_VERSION, COMPILER_OPTIONS
    Use SET_RK
    Use FADDEYEVA_V2_MOD_RK, Only: FADDEYEVA_V2_RK
    Use FADDEYEVA_V3_MOD_RK, Only: FADDEYEVA_V3_RK, ERFC_Z, ERF_Z, ERFI_Z, ERFCX_Z, DAWSON_Z

    Implicit None
    ! ..
    ! .. PARAMETERS ..
    Integer, Parameter :: NMAX = 250000, NREPEATS = 1 + 4*(1-(QP/RK)*(RK/QP))
    Complex (RK), Parameter :: J1 = (0.0E0_RK, 1.0E0_RK)

    ! .. LOCAL SCALARS ..
    Real (RK) :: TIME_BEGIN, TIME_END, ELAPSED_TIME, DUMY1, DUMY2
    Integer :: IM, JN, NDGTS, JN_ERRMAX, NDGTS_MAX
    Real (RK) :: ERR_RE_W_MAX, ERR_IM_W_MAX, ERR_W_MAX_OLD

    ! ..
    ! .. LOCAL ARRAYS ..
    Complex (RK), Dimension (NMAX) :: Z1, W1, Z2
    Complex (RK) :: WTMP
    Real (RK), Dimension (NMAX) :: XTEST, YTEST, RE_W, IM_W, RE_REF, IM_REF

    ! .. INTRINSIC FUNCTIONS ..
    Intrinsic CMPLX, CPU_TIME, REAL
  
!   Call FADDEYEVA_V3_RK((6.800000000000000e+00_rk,3.000000000000000e-02_rk), WTMP, 14)
!   write(*,*) WTMP
!   stop
  
    Open (Unit=12, File='Faddeyeva_ref_Maple.txt', Status='OLD')
!    Open (Unit=20, File='Re_Erfcx_err_3d_dp.txt', Status='Unknown')
!    Open (Unit=21, File='Im_Erfcx_err_3d_dp.txt', Status='Unknown')
    Do JN = 1, NMAX
      Read (12, *) XTEST(JN), YTEST(JN), RE_W(JN), IM_W(JN)
      Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
    End Do



    ! :
    Write (*, *) '*************************'
    Write (*, *) 'SAMPLE RESULTS FROM FADDEYEVA_V3_RK'
    Write (*, *) '*************************'
    ! GENERATE BEST ACCURACY (PRECISION-DEPENDENT), USING PRESENT
    ! VERSION FADDEYEVA_V3_RK WITH
    NDGTS_MAX = 6 + 8*(DP/RK)*(RK/DP) + 24*(RK/QP)*(QP/RK) ! 6, 14, 30

    Call FADDEYEVA_V3_RK(Z1, W1, NDGTS_MAX)

    ERR_RE_W_MAX = 0.0E0_RK
    ERR_IM_W_MAX = 0.0E0_RK
    JN_ERRMAX = 1
    Do JN = 1, NMAX
      If (ABS(REAL(RE_W(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_W(JN))>0.0E0_RK) Then
        ERR_W_MAX_OLD = MAX(ERR_RE_W_MAX, ERR_IM_W_MAX)
        ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_W(JN))/RE_W(JN)))
        ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_W(JN))/IM_W(JN)))
        ! !        write(20,'(ES37.30E3,2x,ES37.30E3,2x,ES37.30E3)')xtest(jn),
        ! ytest(jn),abs((real(w1(jn))-re_w(jn))/re_w(jn))
        ! !        write(21,'(ES37.30E3,2x,ES37.30E3,2x,ES37.30E3)')xtest(jn),
        ! ytest(jn),abs((aimag(w1(jn))-im_w(jn))/im_w(jn))

        ! write(20,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
        ! re_w(jn),real(w1(jn)), abs((real(w1(jn))-re_w(jn))/re_w(jn))
        ! write(21,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
        ! im_w(jn), aimag(w1(jn)),abs((aimag(w1(jn))-im_w(jn))/im_w(jn))


        If (ERR_IM_W_MAX>ERR_W_MAX_OLD .Or. ERR_RE_W_MAX>ERR_W_MAX_OLD) Then
          JN_ERRMAX = JN
        End If
      End If
    End Do

    Write (*, '(A24,I7)') 'JN_ERRMAX           =', JN_ERRMAX

    If (RK/=QP) Then
      Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
        AIMAG(Z1(JN_ERRMAX))
      Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') '  W(Z(JN_ERRMAX))     =', REAL(W1(JN_ERRMAX)), &
        AIMAG(W1(JN_ERRMAX))
    Else
      Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
        AIMAG(Z1(JN_ERRMAX))
      Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  W(Z(JN_ERRMAX))     =', REAL(W1(JN_ERRMAX)), &
        AIMAG(W1(JN_ERRMAX))
    End If

    Write (*, '(A24,ES11.3E3, 17X,A24,ES11.3E3)') 'MAX_ERR_RE_W PRESENT=', ERR_RE_W_MAX, &
      'MAX_ERR_IM_W PRESENT=', ERR_IM_W_MAX



    If (RK/=QP) Then
      Call FADDEYEVA_V2_RK(Z1, W1, NDGTS_MAX)

      ERR_RE_W_MAX = 0.0E0_RK
      ERR_IM_W_MAX = 0.0E0_RK
      Do JN = 1, NMAX
        If (ABS(REAL(RE_W(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_W(JN))>0.0E0_RK) Then
          ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_W(JN))/RE_W(JN)))
          ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_W(JN))/IM_W(JN)))
        End If
      End Do

      Write (*, '(A24,ES15.7E3, 17X,A24,ES15.7E3)') 'MAX_ERR_RE_W ALG916 =', ERR_RE_W_MAX, &
        'MAX_ERR_IM_W ALG916 =', ERR_IM_W_MAX

    End If


 
    Close (Unit=12)


    ! ***********************************************************************
    ! TIMING USING THE MAPLE TABULATED VALUES
    ! ***********************************************************************
    Write (*, *) ' '
    Write (*, '(/A,I5)') 'RESULTS: THE CASE OF MAPLE DATA'

    Do NDGTS = NDGTS_MAX, NDGTS_MAX ! 4 + 11*(1-(sp/rk)*(rk/sp))*(1-(dp/rk)*(dp/rk)), ndgts_max
      Write (*, *) '  '
      Write (*, *) '  '

      Open (Unit=12, File='Faddeyeva_ref_Maple.txt', Status='OLD')

      Do JN = 1, NMAX
        Read (12, *) XTEST(JN), YTEST(JN), RE_W(JN), IM_W(JN)
        Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
      End Do




      W1 = (0.0E0_RK, 0.0E0_RK)

      Call CPU_TIME(TIME_BEGIN)
      Do IM = 1, NREPEATS
        Call FADDEYEVA_V3_RK(Z1, W1, NDGTS)
      End Do
      Call CPU_TIME(TIME_END)
      ELAPSED_TIME = TIME_END - TIME_BEGIN
      Write (*, '(/A,2I5,A,G13.6,A)') 'ZAGHLOUL PRESENT', NDGTS, RK, '  ELAPSED TIME =', &
        ELAPSED_TIME, ' PROCESSOR DEPENDENT UNITS'
      Write (*, '(A,G11.4,A)') 'AVERAGE TIME PER CALL (NS) =', 1.0E9_RK*ELAPSED_TIME/(REAL(NMAX,KIND &
        =RK)*REAL(NREPEATS,KIND=RK)), ' PROCESSOR DEPENDENT UNITS'
      ERR_RE_W_MAX = 0.0E0_RK
      ERR_IM_W_MAX = 0.0E0_RK
      Do JN = 1, NMAX
        If (ABS(REAL(RE_W(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_W(JN))>0.0E0_RK) Then
          ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_W(JN))/RE_W(JN)))
          ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_W(JN))/IM_W(JN)))
        End If
      End Do

      Write (*, '(A16,ES11.3E3, 5X,A16,ES11.3E3)') 'MAX_ERR_RE_W PRESENT = ', ERR_RE_W_MAX, &
        'MAX_ERR_IM_W PRESENT = ', ERR_IM_W_MAX


      Write (*, *) '  '
      If (RK/=QP) Then
        W1 = (0.0E0_RK, 0.0E0_RK)

        Call CPU_TIME(TIME_BEGIN)
        Do IM = 1, NREPEATS
          Call FADDEYEVA_V2_RK(Z1, W1, NDGTS)
        End Do
        Call CPU_TIME(TIME_END)
        ELAPSED_TIME = TIME_END - TIME_BEGIN
        Write (*, '(/A,2I5,A,G13.6,A)') 'ZAGHLOUL 916    ', NDGTS, RK, '  ELAPSED TIME (S) =', &
          ELAPSED_TIME, ' PROCESSOR DEPENDENT UNITS'
        Write (*, '(A,G11.4,A)') 'AVERAGE TIME PER CALL (NS) =', 1.0E9_RK*ELAPSED_TIME/(REAL(NMAX, &
          KIND=RK)*REAL(NREPEATS,KIND=RK)), ' PROCESSOR DEPENDENT UNITS'
        ERR_RE_W_MAX = 0.0E0_RK
        ERR_IM_W_MAX = 0.0E0_RK
        Do JN = 1, NMAX
          If (ABS(REAL(RE_W(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_W(JN))>0.0E0_RK) Then
            ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_W(JN))/RE_W(JN)))
            ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_W(JN))/IM_W(JN)))
          End If
        End Do

        Write (*, '(A16,ES15.7E3, 5X,A16,ES15.7E3)') 'MAX_ERR_RE_W 916=', ERR_RE_W_MAX, &
          'MAX_ERR_IM_W 916=', ERR_IM_W_MAX

      End If
      Write (*, *) '  '
      Write (*, *) '  '




      ! :-
      ! :------------------------------- Accuracy Verification of Erfc(z)---------

      Open (Unit=13, File='Reference_data_erfcz_25_25.txt', Status='old')
      Do JN = 1, NMAX
        Read (13, *) XTEST(JN), YTEST(JN), RE_REF(JN), IM_REF(JN)
        Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
      End Do


      W1 = (0.0E0_RK, 0.0E0_RK)
      W1 = ERFC_Z(Z1, NDGTS)
      ERR_RE_W_MAX = 0.0E0_RK
      ERR_IM_W_MAX = 0.0E0_RK
      JN_ERRMAX = 1
      Do JN = 1, NMAX
        If (ABS(REAL(RE_REF(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_REF(JN))>0.0E0_RK) Then
          ERR_W_MAX_OLD = MAX(ERR_RE_W_MAX, ERR_IM_W_MAX)
          ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_REF(JN))/RE_REF(JN)))
          ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_REF(JN))/IM_REF(JN)))

          ! write(20,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn),re_ref(jn),&
          ! real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn), im_ref(jn), &
          ! aimag(w1(jn)), abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))

          ! write(20,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! re_ref(jn),real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! im_ref(jn), aimag(w1(jn)),abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))

          If (ERR_IM_W_MAX>ERR_W_MAX_OLD .Or. ERR_RE_W_MAX>ERR_W_MAX_OLD) Then
            JN_ERRMAX = JN
          End If
        End If
      End Do



      Write (*, '(A24,I7)') 'JN_ERRMAX           =', JN_ERRMAX

      If (RK/=QP) Then
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') '  erfcz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))
      Else
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfcz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))

        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfcz_ref(Z(JN_ERRMAX))     =', &
          RE_REF(JN_ERRMAX), IM_REF(JN_ERRMAX)
      End If

      Write (*, '(A24,ES11.3E3, 17X,A24,ES11.3E3)') 'MAX_ERR_RE_erfcz PRESENT=', ERR_RE_W_MAX, &
        'MAX_ERR_IM_erfcz PRESENT=', ERR_IM_W_MAX

      ! :---------------------------
      ! :-

      Write (*, *) ' '
      Write (*, *) ' '

      ! :------------------------------- Accuracy Verification of Erf(z)---------

      Open (Unit=14, File='Reference_data_erfz_25_25.txt', Status='old')
      Do JN = 1, NMAX
        Read (14, *) XTEST(JN), YTEST(JN), RE_REF(JN), IM_REF(JN)
        Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
      End Do

      W1 = (0.0E0_RK, 0.0E0_RK)
      W1 = ERF_Z(Z1, NDGTS)
      ERR_RE_W_MAX = 0.0E0_RK
      ERR_IM_W_MAX = 0.0E0_RK
      JN_ERRMAX = 1
      Do JN = 1, NMAX
        If (ABS(REAL(RE_REF(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_REF(JN))>0.0E0_RK) Then
          ERR_W_MAX_OLD = MAX(ERR_RE_W_MAX, ERR_IM_W_MAX)
          ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_REF(JN))/RE_REF(JN)))
          ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_REF(JN))/IM_REF(JN)))

          ! !        write(20,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn),re_ref(jn),&
          ! !        real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! !        write(21,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn), im_ref(jn), &
          ! !        aimag(w1(jn)), abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))
          ! write(20,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! re_ref(jn),real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! im_ref(jn), aimag(w1(jn)),abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))

          If (ERR_IM_W_MAX>ERR_W_MAX_OLD .Or. ERR_RE_W_MAX>ERR_W_MAX_OLD) Then
            JN_ERRMAX = JN
          End If
        End If
      End Do


      Write (*, '(A24,I7)') 'JN_ERRMAX           =', JN_ERRMAX

      If (RK/=QP) Then
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') '  erfz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))
      Else
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))

        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfz_ref(Z(JN_ERRMAX))     =', &
          RE_REF(JN_ERRMAX), IM_REF(JN_ERRMAX)
      End If

      Write (*, '(A24,ES11.3E3, 17X,A24,ES11.3E3)') 'MAX_ERR_RE_erfz PRESENT=', ERR_RE_W_MAX, &
        'MAX_ERR_IM_erfz PRESENT=', ERR_IM_W_MAX

      ! :---------------------------
      ! :-

      Write (*, *) ' '
      Write (*, *) ' '

      ! :------------------------------- Accuracy Verification of Erfi(z)---------

      Open (Unit=15, File='Reference_data_erfiz_25_25.txt', Status='old')
      Do JN = 1, NMAX
        Read (15, *) XTEST(JN), YTEST(JN), RE_REF(JN), IM_REF(JN)
        Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
      End Do

      W1 = (0.0E0_RK, 0.0E0_RK)
      W1 = ERFI_Z(Z1, NDGTS)
      ERR_RE_W_MAX = 0.0E0_RK
      ERR_IM_W_MAX = 0.0E0_RK
      JN_ERRMAX = 1
      Do JN = 1, NMAX
        If (ABS(REAL(RE_REF(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_REF(JN))>0.0E0_RK) Then
          ERR_W_MAX_OLD = MAX(ERR_RE_W_MAX, ERR_IM_W_MAX)
          ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_REF(JN))/RE_REF(JN)))
          ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_REF(JN))/IM_REF(JN)))
          ! write(20,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn),re_ref(jn),&
          ! real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn), im_ref(jn), &
          ! aimag(w1(jn)), abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))
          ! write(20,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! re_ref(jn),real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! im_ref(jn), aimag(w1(jn)),abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))

          If (ERR_IM_W_MAX>ERR_W_MAX_OLD .Or. ERR_RE_W_MAX>ERR_W_MAX_OLD) Then
            JN_ERRMAX = JN
          End If
        End If
      End Do


      Write (*, '(A24,I7)') 'JN_ERRMAX           =', JN_ERRMAX
      If (RK/=QP) Then
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') '  erfiz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))
      Else
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfiz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))

        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfiz_ref(Z(JN_ERRMAX))     =', &
          RE_REF(JN_ERRMAX), IM_REF(JN_ERRMAX)
      End If

      Write (*, '(A24,ES11.3E3, 17X,A24,ES11.3E3)') 'MAX_ERR_RE_erfi_z PRESENT=', ERR_RE_W_MAX, &
        'MAX_ERR_IM_erfi_z PRESENT=', ERR_IM_W_MAX

      ! :---------------------------
      ! :-


      Write (*, *) ' '
      Write (*, *) ' '

      ! :------------------------------- Accuracy Verification of Daw(z)---------

      Open (Unit=16, File='Reference_data_dawz_12_12.txt', Status='old')
      Do JN = 1, NMAX
        Read (16, *) XTEST(JN), YTEST(JN), RE_REF(JN), IM_REF(JN)
        Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
      End Do

      W1 = (0.0E0_RK, 0.0E0_RK)
      W1 = DAWSON_Z(Z1, NDGTS)
      ERR_RE_W_MAX = 0.0E0_RK
      ERR_IM_W_MAX = 0.0E0_RK
      JN_ERRMAX = 1
      Do JN = 1, NMAX
        If (ABS(REAL(RE_REF(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_REF(JN))>0.0E0_RK) Then
          ERR_W_MAX_OLD = MAX(ERR_RE_W_MAX, ERR_IM_W_MAX)
          ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_REF(JN))/RE_REF(JN)))
          ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_REF(JN))/IM_REF(JN)))

          ! write(20,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn),re_ref(jn),&
          ! real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn), im_ref(jn), &
          ! aimag(w1(jn)), abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))
          ! write(20,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! re_ref(jn),real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(1x,ES23.15E3,1x))')xtest(jn), ytest(jn), &
          ! im_ref(jn), aimag(w1(jn)),abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))

          If (ERR_IM_W_MAX>ERR_W_MAX_OLD .Or. ERR_RE_W_MAX>ERR_W_MAX_OLD) Then
            JN_ERRMAX = JN
          End If
        End If
      End Do


      Write (*, '(A24,I7)') 'JN_ERRMAX           =', JN_ERRMAX
      If (RK/=QP) Then
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') '  dawsonz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))
      Else
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  dawsonz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))

        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  dawsonz_ref(Z(JN_ERRMAX))     =', &
          RE_REF(JN_ERRMAX), IM_REF(JN_ERRMAX)
      End If

      Write (*, '(A24,ES11.3E3, 17X,A24,ES11.3E3)') 'MAX_ERR_RE_dawsonz PRESENT=', ERR_RE_W_MAX, &
        'MAX_ERR_IM_dawsonz PRESENT=', ERR_IM_W_MAX

      ! :---------------------------
      ! :-

      Write (*, *) ' '
      Write (*, *) ' '

      ! :------------------------------- Accuracy Verification of Erfcx(z)---------

      Open (Unit=17, File='Reference_data_erfcxz_12_12.txt', Status='old')
      Do JN = 1, NMAX
        Read (17, *) XTEST(JN), YTEST(JN), RE_REF(JN), IM_REF(JN)
        Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
      End Do

      W1 = (0.0E0_RK, 0.0E0_RK)
      W1 = ERFCX_Z(Z1, NDGTS)
      ERR_RE_W_MAX = 0.0E0_RK
      ERR_IM_W_MAX = 0.0E0_RK
      JN_ERRMAX = 1
      Do JN = 1, NMAX
        If (ABS(REAL(RE_REF(JN),KIND=RK))>0.0E0_RK .And. ABS(IM_REF(JN))>0.0E0_RK) Then
          ERR_W_MAX_OLD = MAX(ERR_RE_W_MAX, ERR_IM_W_MAX)
          ERR_RE_W_MAX = MAX(ERR_RE_W_MAX, ABS((REAL(W1(JN))-RE_REF(JN))/RE_REF(JN)))
          ERR_IM_W_MAX = MAX(ERR_IM_W_MAX, ABS((AIMAG(W1(JN))-IM_REF(JN))/IM_REF(JN)))
          ! write(20,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn),re_ref(jn),&
          ! real(w1(jn)), abs((real(w1(jn))-re_ref(jn))/re_ref(jn))
          ! write(21,'(5(ES39.30E3,2x))')xtest(jn), ytest(jn), im_ref(jn), &
          ! aimag(w1(jn)), abs((aimag(w1(jn))-im_ref(jn))/im_ref(jn))
          Write (20, '(5(1x,ES23.15E3,1x))') XTEST(JN), YTEST(JN), RE_REF(JN), REAL(W1(JN)), &
            ABS((REAL(W1(JN))-RE_REF(JN))/RE_REF(JN))
          Write (21, '(5(1x,ES23.15E3,1x))') XTEST(JN), YTEST(JN), IM_REF(JN), AIMAG(W1(JN)), &
            ABS((AIMAG(W1(JN))-IM_REF(JN))/IM_REF(JN))

          If (ERR_IM_W_MAX>ERR_W_MAX_OLD .Or. ERR_RE_W_MAX>ERR_W_MAX_OLD) Then
            JN_ERRMAX = JN
          End If
        End If
      End Do


      Write (*, '(A24,I7)') 'JN_ERRMAX           =', JN_ERRMAX
      If (RK/=QP) Then
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES23.15E3, 2X,ES23.15E3)') '  erfcxz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))
      Else
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') 'Z(JN_ERRMAX)        =', REAL(Z1(JN_ERRMAX)), &
          AIMAG(Z1(JN_ERRMAX))
        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfcxz(Z(JN_ERRMAX))     =', &
          REAL(W1(JN_ERRMAX)), AIMAG(W1(JN_ERRMAX))

        Write (*, '(A24,ES38.30E3, 2X,ES38.30E3)') '  erfcxz_ref(Z(JN_ERRMAX))     =', &
          RE_REF(JN_ERRMAX), IM_REF(JN_ERRMAX)
      End If

      Write (*, '(A24,ES11.3E3, 17X,A24,ES11.3E3)') 'MAX_ERR_RE_erfcxz PRESENT=', ERR_RE_W_MAX, &
        'MAX_ERR_IM_erfcxz PRESENT=', ERR_IM_W_MAX




      Close (Unit=12)
      Close (Unit=13)
      Close (Unit=14)
      Close (Unit=15)
      Close (Unit=16)
      Close (Unit=17)

      Close (Unit=20)
      Close (Unit=21)

    End Do

  End Program FADDEYEVA_V3_DRIVER_RK
