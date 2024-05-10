      SUBROUTINE GCT_HUMLIK ( N, X1, DX, Y, K )

c     To calculate the Voigt/Faddeeva function with relative error less than 10^(-4).
c  Based on the paper:
c        Rapid approximation to the Voigt/Faddeeva function and its derivatives 
c        Wells R.J., JQSRT, Volume 62, Number 1, May 1999 , pp. 29-48(20)
c        with corrections.

* Arguments
      implicit none
      INTEGER N                                       ! IN   Number of points
c      REAL    X(0:N-1)                               ! IN   Input x array
      REAL*8  X1, DX, XX
      REAL    Y                                       ! IN   Input y value >=0.0
      REAL    K(0:N-1)                                ! OUT  Real (Voigt) array

* Constants
      REAL        RRTPI                               ! 1/SQRT(pi)
      PARAMETER ( RRTPI = 0.56418958 )
      REAL        Y0,       Y0PY0,         Y0Q        ! for CPF12 algorithm
      PARAMETER ( Y0=1.5, Y0PY0=Y0+Y0, Y0Q=Y0*Y0 )
      REAL  C(0:5), S(0:5), T(0:5)
      SAVE  C, S, T    ! preserves values of C, S and T arrays between calls
* Local variables
      INTEGER I, J                                   ! Loop variables
      INTEGER RG1, RG2, RG3                          ! y polynomial flags
      REAL ABX, XQ, YQ, YRRTPI, xfudge               ! |x|, x^2, y^2, y/SQRT(pi)
      REAL XLIM0, XLIM1, XLIM2, XLIM3, XLIM4         ! |x| on region boundaries
      REAL A0, D0, D2, E0, E2, E4, H0, H2, H4, H6    ! W4 temporary variables
      REAL P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
      REAL XP(0:5), XM(0:5), YP(0:5), YM(0:5)        ! CPF12 temporary values
      REAL MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
      REAL D, YF, YPY0, YPY0Q  
c
      DATA C / 1.0117281,     -0.75197147,        0.012557727,
     &         0.010022008,   -0.00024206814,     0.00000050084806 /
      DATA S / 1.393237,       0.23115241,       -0.15535147,
     &         0.0062183662,   0.000091908299,   -0.00000062752596 /
      DATA T / 0.31424038,     0.94778839,        1.5976826,
     &         2.2795071,      3.0206370,         3.8897249 /

***** Start of executable code *****************************************

      XX=X1
      YQ  = Y*Y                                      ! y^2
      YRRTPI = Y*RRTPI                               ! y/SQRT(pi)

      IF ( Y .GE. 70.55 ) THEN                       ! All points
        DO I = 0, N-1                                ! Loop over all points
          XQ   = XX*XX
          xfudge=(xx/(10000*dx))**2
          K(I) = YRRTPI / (YQ + XQ*(1+xfudge))
          xx=xx+dx
        ENDDO
        RETURN
      ENDIF


      RG1 = 1                                        ! Set flags
      RG2 = 1
      RG3 = 1

      XLIM0 = SQRT ( 15100.0 + Y*(40.0 - Y*3.6) )
      IF ( Y .GE. 8.425 ) THEN
         XLIM1 = 0.0
      ELSE
         XLIM1 = SQRT ( 164.0 - Y*(4.3 + Y*1.8) )
      ENDIF

      XLIM2 = 6.8 - Y
      XLIM3 = 2.4*Y
      XLIM4 = 18.1*Y + 1.65
      IF ( Y .LE. 0.000001 ) THEN                    ! When y<10^-6
       XLIM1 = XLIM0                                 ! avoid W4 algorithm
       XLIM2 = XLIM0
      ENDIF
*.....
      DO I = 0, N-1                                  ! Loop over all points
c       ABX = ABS ( X(I) )                           ! |x|
       abx  =  abs(XX)
       XQ  =   XX*XX                                 ! x^2
       xfudge=(xx/(10000*dx))**2
       IF     ( ABX .GE. XLIM0 ) THEN                ! Region 0 algorithm
          K(I) = YRRTPI / (YQ + XQ*(1+xfudge))

       ELSEIF ( ABX .GE. XLIM1 ) THEN                ! Humlicek W4 Region 1
        IF ( RG1 .NE. 0 ) THEN                       ! First point in Region 1
         RG1 = 0
         A0 = YQ + 0.5                               ! Region 1 y-dependents
         D0 = A0*A0
         D2 = YQ + YQ - 1.0
        ENDIF
        D = RRTPI / (D0 + XQ*(D2 + XQ))
        K(I) = D*Y   *(A0 + XQ)

       ELSEIF ( ABX .GT. XLIM2 ) THEN                ! Humlicek W4 Region 2 
        IF ( RG2 .NE. 0 ) THEN                       ! First point in Region 2
         RG2 = 0
         H0 = 0.5625+YQ*(4.5+YQ*(10.5+YQ*(6.0+YQ)))   ! Region 2 y-dependents
         H2 = -4.5    + YQ*(9.0 + YQ*( 6.0 + YQ* 4.0))
         H4 = 10.5    - YQ*(6.0 - YQ*  6.0)
         H6 = -6.0    + YQ* 4.0
         E0 =  1.875  + YQ*(8.25 + YQ*(5.5 + YQ))
         E2 =  5.25   + YQ*(1.0  + YQ* 3.0)
         E4 =  0.75*H6
        ENDIF
        D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
        K(I) = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)))

       ELSEIF ( ABX .LT. XLIM3 ) THEN                 ! Humlicek W4 Region 3
        IF ( RG3 .NE. 0 ) THEN                        ! First point in Region 3
         RG3 = 0
         Z0 = 272.1014 + Y*(1280.829 + Y*(2802.870 + Y*(3764.966  ! Region 3 y-dependents
     &       + Y*(3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989
     &       + Y*(88.26741 + Y*(13.39880 + Y)))))))))
         Z2 = 211.678      + Y*(902.3066 + Y*(1758.336 + Y*(2037.310
     &                     + Y*(1549.675 + Y*(793.4273 + Y*(266.2987
     &                     + Y*(53.59518 + Y*5.0)))))))
         Z4 = 78.86585     + Y*(308.1852 + Y*(497.3014 + Y*(479.2576
     &                     + Y*(269.2916 + Y*(80.39278 + Y*10.0)))))
         Z6 = 22.03523     + Y*(55.02933 + Y*(92.75679 + Y*(53.59518
     &                     + Y*10.0)))
         Z8 = 1.496460     + Y*(13.39880 + Y*5.0)
         P0 = 153.5168     + Y*(549.3954 + Y*(919.4955 + Y*(946.8970
     &       + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941
     &                     + Y*(4.264678 + Y*0.3183291))))))))
         P2 = -34.16955    + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730
     &                     + Y*(139.4665 + Y*(56.81652 + Y*(12.79458
     &                     + Y*1.2733163))))))
         P4 = 2.584042     + Y*(10.46332 + Y*(24.01655 + Y*(29.81482
     &                     + Y*(12.79568 + Y*1.9099744))))
         P6 = -0.07272979  + Y*(0.9377051+ Y*(4.266322 + Y*1.273316))
         P8 = 0.0005480304 + Y*0.3183291
        ENDIF
        D = 1.7724538 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))))
        K(I) = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))

       ELSE                                          ! Humlicek CPF12 algorithm
        YPY0 = Y + Y0
        YPY0Q = YPY0*YPY0
        K(I) = 0.0
         DO J = 0, 5
c          D = X(I) - T(J)
          D = XX - T(J)
          MQ(J) = D*D
          MF(J) = 1.0 / (MQ(J) + YPY0Q)
          XM(J) = MF(J)*D
          YM(J) = MF(J)*YPY0
c          D = X(I) + T(J)
          D = XX + T(J)
          PQ(J) = D*D
          PF(J) = 1.0 / (PQ(J) + YPY0Q)
          XP(J) = PF(J)*D
          YP(J) = PF(J)*YPY0
         ENDDO

        IF ( ABX .LE. XLIM4 ) THEN                   ! Humlicek CPF12 Region I
         DO J = 0, 5
          K(I) = K(I) + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
         ENDDO

        ELSE                                         ! Humlicek CPF12 Region II
         YF   = Y + Y0PY0
         DO J = 0, 5
          K(I) = K(I)
     &     + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)) / (MQ(J)+Y0Q)
     &     + (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF*XP(J)) / (PQ(J)+Y0Q)
         ENDDO
         K(I) = Y*K(I) + EXP ( -XQ )
        ENDIF
       ENDIF
       XX=XX+DX
      ENDDO
*.....
      END
