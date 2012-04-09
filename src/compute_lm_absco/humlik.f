**********************************************************************
      SUBROUTINE HUMLIK ( N, X, Y, K, L )
C*********************************************************************
C "HUMLIK": Complex Probability Function
C .........................................................
C         .       Subroutine to Compute the Complex       .
C         .        Probability Function W(z=X+iY)         .
C         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
C         .    Which Appears when Convoluting a Complex   .
C         .     Lorentzian Profile by a Gaussian Shape    .
C         .................................................
C
C		    N : Number of points 
C		    z = X+iY
C             K : Real Part of W(z)
C             L : Imaginary Part of W(z)
C
C This Routine was Taken from the Paper by R.J. Wells, page 29, Volume
C 62 of the 1999 Issue 1 of the "Journal of Quantitative Spectroscopy
C and Radiative Transfer"
C Please Refer to this Paper for More Information
C
C Accessed Files:  None
C --------------
C
C Called Routines: None                               
C ---------------                                 
C
C Called By: 'CompAbs' (COMPute ABSorpton)
C ---------
C
C Double Precision Version
C
C F. Niro, last change 15 Jan 2005
C*********************************************************************
C      

      real*8        R0,         R1 ! Region boundaries
      PARAMETER ( R0 = 146.7d0, R1 = 14.67d0 ) ! for R=4

* Arguments
      INTEGER N                 ! IN   Number of points
      real*8    X(0:N-1)        ! IN   Input x array
      real*8    Y               ! IN   Input y value >=0.0
      real*8    K(0:N-1)        ! OUT  real*8 (Voigt) array
      real*8    l(0:n-1)        ! OUT  Optional array

* Constants
      real*8        RRTPI       ! 1/SQRT(pi)
      PARAMETER ( RRTPI = 0.56418958d0 )
      real*8        Y0,       Y0PY0,         Y0Q ! for CPF12 algorithm
      PARAMETER ( Y0 = 1.5d0, Y0PY0 = Y0+Y0, Y0Q = Y0*Y0  )
      real*8  C(0:5), S(0:5), T(0:5)
      SAVE  C,      S,      T
*     SAVE preserves values of C, S and T (static) arrays between
c procedure calls
      DATA C / 1.0117281d0,     -0.75197147d0,        0.012557727d0,
     &     0.010022008d0,   -0.00024206814d0,    0.00000050084806d0/
      DATA S / 1.393237d0,       0.23115241d0,       -0.15535147d0,
     &     0.0062183662d0,   0.000091908299d0,  -0.00000062752596d0/
      DATA T / 0.31424038d0,     0.94778839d0,        1.5976826d0,
     &     2.2795071d0,      3.0206370d0,         3.8897249d0 /

* Local variables
      INTEGER I, J              ! Loop variables
      INTEGER RG1, RG2, RG3     ! y polynomial flags
      real*8 ABX, XQ, YQ, YRRTPI ! |x|, x^2, y^2, y/SQRT(pi)
      real*8 XLIM0, XLIM1, XLIM2, XLIM3, XLIM4 ! |x| on region boundaries
      real*8 A0, D0, D2, E0, E2, E4, H0, H2, H4, H6 ! W4 temporary variables
      real*8 P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
      real*8 B1, F1, F3, F5, Q1, Q3, Q5, Q7
      real*8 XP(0:5), XM(0:5), YP(0:5), YM(0:5) ! CPF12 temporary values
      real*8 MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
      real*8 D, YF, YPY0, YPY0Q  

***** Start of executable code *****************************************

      RG1 = 1                   ! Set flags
      RG2 = 1
      RG3 = 1
      YQ  = Y*Y                 ! y^2
      YRRTPI = Y*RRTPI          ! y/SQRT(pi)

*     Region boundaries when both K and L are required or when R<>4 
      XLIM0 = R0 - Y 
      XLIM1 = R1 - Y
      XLIM3 = 3.097d0*Y - 0.45d0
*     For speed the following 3 lines should replace the 3 above if R=4
c and L is not required   *
*     XLIM0 = 15100.0 + Y*(40.0 + Y*3.6)
c *
*     XLIM1 = 164.0 - Y*(4.3 + Y*1.8)
c *
*     XLIM3 = 5.76*YQ
c *

      XLIM2 = 6.8d0 - Y
      XLIM4 = 18.1d0*Y + 1.65d0
      IF ( Y .LE. 0.000001d0 ) THEN ! When y<10^-6
        XLIM1 = XLIM0           ! avoid W4 algorithm
        XLIM2 = XLIM0
      ENDIF
*.....
      DO I = 0, N-1             ! Loop over all points
        ABX = dABS ( X(I) )     ! |x|
        XQ  = ABX*ABX           ! x^2
        IF     ( ABX .GT. XLIM0 ) THEN ! Region 0 algorithm
          K(I) = YRRTPI / (XQ + YQ)
          L(I) = K(I)*X(I) / Y

        ELSEIF ( ABX .GT. XLIM1 ) THEN ! Humlicek W4 Region 1
          IF ( RG1 .NE. 0 ) THEN ! First point in Region 1
            RG1 = 0d0
            A0 = YQ + 0.5d0     ! Region 1 y-dependents
            D0 = A0*A0
            D2 = YQ + YQ - 1.0d0
            B1 = YQ - 0.5d0
          ENDIF
          D = RRTPI / (D0 + XQ*(D2 + XQ))
          K(I) = D*Y   *(A0 + XQ)
          L(I) = D*X(I)*(B1 + XQ)

        ELSEIF ( ABX .GT. XLIM2 ) THEN ! Humlicek W4 Region 2 
          IF ( RG2 .NE. 0 ) THEN ! First point in Region 2
            RG2 = 0
            H0 = 0.5625d0 + YQ*(4.5d0 + YQ*(10.5d0 + YQ*(6.0d0 + YQ))) ! Region 2 y-dependents
            H2 = -4.5d0    + YQ*(9.0d0 + YQ*( 6.0d0 + YQ* 4.0d0))
            H4 = 10.5d0    - YQ*(6.0d0 - YQ*  6.0d0)
            H6 = -6.0d0    + YQ* 4.0d0
            E0 =  1.875d0  + YQ*(8.25d0 + YQ*(5.5d0 + YQ))
            E2 =  5.25d0   + YQ*(1.0d0  + YQ* 3.0d0)
            E4 =  0.75d0*H6
            F1 = -1.875d0  + YQ*(5.25d0 + YQ*(4.5d0 + YQ))
            F3 =  8.25d0   - YQ*(1.0d0  - YQ* 3.0d0)
            F5 = -5.5d0    + YQ* 3.0d0
          ENDIF
          D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
          K(I) = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)))
          L(I) = D*X(I)*(F1 + XQ*(F3 + XQ*(F5 + XQ)))

        ELSEIF ( ABX .LT. XLIM3 ) THEN ! Humlicek W4 Region 3
          IF ( RG3 .NE. 0 ) THEN ! First point in Region 3
            RG3 = 0
            Z0 = 272.1014d0 + Y*(1280.829d0 + Y*(2802.870d0 + Y*(3764 ! Region 3 y-dependents
     &           .966d0+ Y*(3447.629d0+ Y*(2256.981d0 + Y*(1074.409d0 +
     &           Y*(369.1989d0+ Y*(88.26741d0+ Y*(13.39880d0 + Y))))))))
     &           )
            Z2 = 211.678d0  + Y*(902.3066d0 + Y*(1758.336d0 + Y*(2037
     &           .310d0+ Y*(1549.675d0 + Y*(793.4273d0 + Y*(266.2987d0+
     &           Y*(53.59518d0 + Y*5.0d0)))))))
            Z4 = 78.86585d0 + Y*(308.1852d0 + Y*(497.3014d0 + Y*(479
     &           .2576d0+ Y*(269.2916d0 + Y*(80.39278d0 + Y*10.0d0)))))
            Z6 = 22.03523d0 + Y*(55.02933d0 + Y*(92.75679d0 + Y*(53
     &           .59518d0+ Y*10.0d0)))
            Z8 = 1.496460d0   + Y*(13.39880d0 + Y*5.0d0)
            P0 = 153.5168d0 + Y*(549.3954d0 + Y*(919.4955d0 + Y*(946
     &           .8970d0+ Y*(662.8097d0+ Y*(328.2151d0 + Y*(115.3772d0 +
     &           Y*(27.93941d0+ Y*(4.264678d0 + Y*0.3183291d0))))))))
            P2 = -34.16955d0+ Y*(-1.322256d0+ Y*(124.5975d0 + Y*(189
     &           .7730d0+ Y*(139.4665d0 + Y*(56.81652d0 + Y*(12.79458d0+
     &           Y*1.2733163d0))))))
            P4 = 2.584042d0 + Y*(10.46332d0 + Y*(24.01655d0 + Y*(29
     &           .81482d0+ Y*(12.79568d0 + Y*1.9099744d0))))
            P6 = -0.07272979d0+Y*(0.9377051d0+Y*(4.266322d0 +Y*1
     &           .273316d0))
            P8 = 0.0005480304d0 + Y*0.3183291d0
            Q1 = 173.2355d0 + Y*(508.2585d0 + Y*(685.8378d0 + Y*(557
     &           .5178d0+ Y*(301.3208d0 + Y*(111.0528d0 + Y*(27.62940d0+
     &           Y*(4.264130d0 + Y*0.3183291d0)))))))
            Q3 = 18.97431d0 + Y*(100.7375d0 + Y*(160.4013d0 + Y*(130
     &           .8905d0+ Y*(55.88650d0 + Y*(12.79239d0+Y*1.273316d0))))
     &           )
            Q5 = 7.985877d0 + Y*(19.83766d0 + Y*(28.88480d0 + Y*(12
     &           .79239d0+ Y*1.909974d0)))
            Q7 = 0.6276985d0    + Y*(4.264130d0 + Y*1.273316d0)
          ENDIF
          D =1.7724538d0 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))
     &         ))
          K(I) =D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
          L(I) =D*X(I)*(Q1 + XQ*(Q3 + XQ*(Q5 + Xq*(Q7 + XQ*0.3183291d0))
     &         ))

        ELSE                    ! Humlicek CPF12 algorithm
          YPY0 = Y + Y0
          YPY0Q = YPY0*YPY0
          K(I) = 0.0d0
          L(I) = 0.0d0
          DO J = 0, 5
            D = X(I) - T(J)
            MQ(J) = D*D
            MF(J) = 1.0d0 / (MQ(J) + YPY0Q)
            XM(J) = MF(J)*D
            YM(J) = MF(J)*YPY0
            D = X(I) + T(J)
            PQ(J) = D*D
            PF(J) = 1.0d0 / (PQ(J) + YPY0Q)
            XP(J) = PF(J)*D
            YP(J) = PF(J)*YPY0
            L(I)  = L(I) + C(J)*(XM(J)+XP(J)) + S(J)*(YM(J)-YP(J))
          ENDDO

          IF ( ABX .LE. XLIM4 ) THEN ! Humlicek CPF12 Region I
            DO J = 0, 5
              K(I) = K(I) + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
            ENDDO

          ELSE                  ! Humlicek CPF12 Region II
            YF   = Y + Y0PY0
            DO J = 0, 5
              K(I) = K(I)
     &             + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)) /
     &             (MQ(J)+Y0Q)+ (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF
     &             *XP(J)) / (PQ(J)+Y0Q)
            ENDDO
            K(I) = Y*K(I) + dEXP ( -XQ )
          ENDIF
        ENDIF
      ENDDO
*.....
      END SUBROUTINE HUMLIK
**********************************************************************
