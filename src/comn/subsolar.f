      SUBROUTINE subsolar(OBJECT,YEAR,MON,DAY,XJDFR,sslat,sslon,eorv)
C     Subroutine to get location (sslat,sslon) where object (Sun) is
c     overhead. Also computes Earth-Object radial velocity
C----------------------------------------------------------------------
C     Calculates the sub-solar latitude and longitude (in radians).
C     Solar ephemeris based on trignometric series developed by
C     VAN FLANDERN & PULKINNEN, Ap.J.Suppl. 41,391-411, 1979.
C     Claimed accuracy < 1 arcmin for +/- 300 years.
C----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,L,O-Z)
      INTEGER*4 YEAR,MON,DAY,zero,one,OBJECT,jd,n
C
      PARAMETER (H2R=0.2617993878D0,zero=0,one=1)
C
      CALL JULIAN(YEAR,MON,DAY,JD)
      XJD=DBLE(JD)-0.5D0
C
C     Get the offset from the standard epoch: Year 2000, Jan 1.5
      TAU = XJD+XJDFR - 2451545.0D0
C
C     and convert to Julian centuries from 1900.0
      T = 1.0D0 + TAU / 36525.0D0
C
c  Calculate Earth-Sun Radial velocity
      CALL COORS
     & (OBJECT,TAU-1.0d0,T-1.0d0/36525.d0,ALPHA,DELTA,RHOM,LON,BETA,RP)
      CALL COORS
     & (OBJECT,TAU+1.0d0,T+1.0d0/36525.d0,ALPHA,DELTA,RHOP,LON,BETA,RP)
c  RHOP & RHOM are in AU's
      eorv=865750.d0*(rhop-rhom)   ! Earth-Sun radial velocity in m/s
c
C  Get the solar coordinates
      CALL COORS (OBJECT,TAU,T,ALPHA,DELTA,RHO,LON,BETA,RP)
C
C     Get the Greenwich mean sidereal time for Jan 0.0 of year
      CALL JULIAN (YEAR,one,zero,JD)
      YJD = JD
      TP = (YJD - 2415020.5D0) / 36525.0D0
      YM = 2400.051262D0 * TP            
      N = int(YM/24.0D0)                
      YM = YM - 24.0D0 * N       
C
      GMST = 6.646065556D0 + YM + 2.58056D-05 * TP * TP
C      IF (GMST .GT. 24.0D0) GMST = GMST - 24.0D0
C
C     Correct for the number of days since the year began
      N = int(XJD-YJD+0.5D0) 
      GMST = (6.570982222D-02*N) + (24.06570982D0*XJDFR) + GMST
      GMST=DMOD(GMST,24.D0)
C     Units of GMST are fractional hours
C     Astronomical convention is that larger angles (times)
C     are EAST of smaller ones
      sslat=DELTA
      sslon=ALPHA-GMST*H2R
      RETURN
      END


C     COORDINATE SUBROUTINE
      SUBROUTINE COORS(OBJECT,TAU,T,ALPHA,DELTA,RHO,LON,BETA,RP)
      IMPLICIT REAL*8 (A-H,L,O-Z)
      INTEGER*4 OBJECT
      COMMON PLON,U,V,W,A1,A2,A3,A4,A5,A6,A7,A8,A9,
     &A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,
     &A24,A25,A26,A27,A28,A29,A30,A31,A32,A33
C
      DATA FF/6.28318530718D0/,GG/206264.8062D0/
C
C     GET THE FUNDAMENTAL ARGUMENTS
C
C     MOON
      LM=0.606434d0+0.03660110129D0*TAU
      GM=0.374897d0+0.03629164709D0*TAU
      FM=0.259091d0+0.03674819520D0*TAU
      D =0.827362d0+0.03386319198D0*TAU
      OMEGAM=0.347343d0-.00014709391D0*TAU
C
C     SUN
      LS=0.779072d0+0.00273790931D0*TAU
      GS=0.993126d0+0.00273777850D0*TAU
C
C     MERCURY
      L1=0.700695d0+0.01136771400D0*TAU
      G1=0.485541d0+0.01136759566D0*TAU
      F1=0.566441d0+0.01136762384D0*TAU
C
C     VENUS
      L2=0.505498d0+0.00445046867D0*TAU
      G2=0.140023d0+0.00445036173D0*TAU
      F2=0.292498d0+0.00445040017D0*TAU
C
C     MARS
      L4=0.987353d0+0.00145575328D0*TAU
      G4=0.053856d0+0.00145561327D0*TAU
      F4=0.849694d0+0.00145569465D0*TAU
C
C     JUPITER
      L5=0.089608d0+0.00023080893D0*TAU
      G5=0.056531d0+0.00023080893D0*TAU
      F5=0.814794d0+0.00023080893D0*TAU
C
C     SATURN
      L6=0.133295d0+0.00009294371D0*TAU
      G6=0.882987d0+0.00009294371D0*TAU
      F6=0.821218d0+0.00009294371D0*TAU
C
C     URANUS
      L7=0.870169d0+0.00003269438D0*TAU
      G7=0.400589d0+0.00003269438D0*TAU
      F7=0.664614d0+0.00003265562D0*TAU
C
C     NEPTUNE
      L8=0.846912d0+0.00001672092D0*TAU
      G8=0.725368d0+0.00001672092D0*TAU
      F8=0.480856d0+0.00001663715D0*TAU
C
C     PLUTO
      L9=0.663854d0+0.00001115482D0*TAU
      G9=0.041020d0+0.00001104864D0*TAU
      F9=0.357355d0+0.00001104864D0*TAU
C
C     REDUCE TO FRACTIONAL REVOLUTION AND CONVERT TO RADIANS
C
      A1=FF*(LM-DINT(LM))
      A2=FF*(GM-DINT(GM))
      A3=FF*(FM-DINT(FM))
      A4=FF*(D-DINT(D))
      A5=FF*(OMEGAM-DINT(OMEGAM))
      A7=FF*(LS-DINT(LS))
      A8=FF*(GS-DINT(GS))
      A9=FF*(L1-DINT(L1))
      A10=FF*(G1-DINT(G1))
      A11=FF*(F1-DINT(F1))
      A12=FF*(L2-DINT(L2))
      A13=FF*(G2-DINT(G2))
      A14=FF*(F2-DINT(F2))
      A15=FF*(L4-DINT(L4))
      A16=FF*(G4-DINT(G4))
      A17=FF*(F4-DINT(F4))
      A18=FF*(L5-DINT(L5))
      A19=FF*(G5-DINT(G5))
      A20=FF*(F5-DINT(F5))
      A21=FF*(L6-DINT(L6))
      A22=FF*(G6-DINT(G6))
      A23=FF*(F6-DINT(F6))
      A24=FF*(L7-DINT(L7))
      A25=FF*(G7-DINT(G7))
      A26=FF*(F7-DINT(F7))
      A27=FF*(L8-DINT(L8))
      A28=FF*(G8-DINT(G8))
      A29=FF*(F8-DINT(F8))
      A31=FF*(L9-DINT(F9))
      A32=FF*(G9-DINT(G9))
      A33=FF*(F9-DINT(F9))
c10    FORMAT(1X,I3,F10.5)
C
      GO TO(100,200,300,400,500,600,700,800,900,1000),OBJECT
C
C     MOON CODE
100   CONTINUE
C
      CALL MOON (BETA,RP,T)
      LON=PLON+A1*GG
      L=A1
      DBAR=60.40974D0
      GO TO 2000
C
C     SUN CODE
200   CONTINUE
      CALL SUN (BETA,RP,T)
      LON=PLON+A7*GG
      L=A7
      DBAR=1.00021D0
      GO TO 2000
C
C  MERCURY CODE
300   CONTINUE
cc      CALL MERCURY (BETA,RP,T)
      LON=PLON+A9*GG
      L=A7
      DBAR=1.07693D0
      GO TO 2000
C
C  VENUS CODE
400   CONTINUE
cc      CALL VENUS (BETA,RP,T)
      LON=PLON+A12*GG
      L=A7
      DBAR=1.23437D0
      GO TO 2000
C
C  MARS CODE
500   CONTINUE
cc      CALL MARS (BETA,RP,T)
      LON=PLON+A15*GG
      L=A15
      DBAR=1.83094D0
      GO TO 2000
C
C  JUPITER CODE
600   CONTINUE
cc      CALL JUPITER (BETA,RP,T)
      LON=PLON+A18*GG
      L=A18
      DBAR=5.30693D0
      GO TO 2000
C
C  SATURN CODE
700   CONTINUE
cc      CALL SATURN (BETA,RP,T)
      LON=PLON+A21*GG
      L=A21
      DBAR=9.61711D0
      GO TO 2000
C
C  URANUS CODE
800   CONTINUE
cc      CALL URANUS (BETA,RP,T)
      LON=PLON+A24*GG
      L=A24
      DBAR=19.24877D0
      GO TO 2000
C
C  NEPTUNE CODE
900   CONTINUE
cc      CALL NEPTUNE (BETA,RP,T)
      LON=PLON+A27*GG
      L=A27
      DBAR=30.08900D0
      GO TO 2000
C
C  PLUTO CODE
1000  CONTINUE
cc      CALL PLUTO (BETA,RP,T)
      LON=PLON+A31*GG
      L=A31
      DBAR=41.32680D0
2000  CONTINUE
C
C  GET GEOCENTRIC COORDINATES IN DEGREES
      IF(L.LT.0d0) L=FF+L
      ALPHA=(L+DASIN(W/DSQRT(U-V**2)))
      IF(ALPHA.LT.0d0)ALPHA=ALPHA+FF
      DELTA=DASIN(V/DSQRT(U))
      RHO=DBAR*DSQRT(U)
      LON=LON/3600.d0
      IF(LON.LT.0.)LON=360.d0+LON
      BETA=BETA/3600.d0
c      WRITE(6,*)PLON,RP,V,U,W,ALPHA,DELTA,RHO,LON
      RETURN
      END


C     SUN CODE
      SUBROUTINE SUN (BETA,RP,T)
C
      IMPLICIT REAL*8 (A-H,L,O-Z)
      COMMON PLON,U,V,W,A1,A2,A3,A4,A5,A6,A7,A8,A9,
     &A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,
     &A24,A25,A26,A27,A28,A29,A30,A31,A32,A33
C
C
      PLON=6910*DSIN(A8)+72*DSIN(2*A8)-17*T*DSIN(A8)-7*DCOS(A8-A19)
     A+6*DSIN(A1-A7)+5*DSIN(4*A8-8*A16+3*A19)-5*DCOS(2*A8-2*A13)
     B-4*DSIN(A8-A13)+4*DCOS(4*A8-8*A16+3*A19)+3*DSIN(2*A8-2*A13)
     C-3*DSIN(A19)-3*DSIN(2*A8-2*A19)
C      DONE
      BETA=0.0d0
C      DONE
      RP=1.00014D0-0.01675d0*DCOS(A8)-0.00014d0*DCOS(2*A8)
C      DONE
      V=0.39785D0*DSIN(A7)-0.01d0*DSIN(A7-A8)+0.00333d0*DSIN(A7+A8)
     A-.00021d0*T*DSIN(A7)+0.00004d0*DSIN(A7+2*A8)-0.00004d0*DCOS(A7)
     B-.00004d0*DSIN(A5-A7)+.00003d0*T*DSIN(A7-A8)
C      DONE
      U=1.d0-0.03349d0*DCOS(A8)-0.00014d0*DCOS(2*A8)
     A+0.00008d0*T*DCOS(A8)-0.00003d0*DSIN(A8-A19)
C      DONE
      W=-0.04129d0*DSIN(2*A7)
     &+0.03211d0*DSIN(A8)+0.00104d0*DSIN(2*A7-A8)
     A-0.00035d0*DSIN(2*A7+A8)-0.0001d0-0.00008d0*T*DSIN(A8)
     B-0.00008d0*DSIN(A5)+0.00007d0*DSIN(2*A8)+0.00005d0*T*DSIN(2*A7)
     C+0.00003d0*DSIN(A1-A7)-0.00002d0*DCOS(A8-A19)
     D+0.00002d0*DSIN(4*A8-8*A16+3*A19)-0.00002d0*DSIN(A8-A13)
     E-0.00002d0*DCOS(2*A8-2*A13)
C      DONE
      RETURN
      END


C     MOON CODE
      SUBROUTINE MOON (BETA,RP,T)
C
      IMPLICIT REAL*8 (A-H,L,O-Z)
      COMMON PLON,U,V,W,A1,A2,A3,A4,A5,A6,A7,A8,A9,
     &A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,
     &A24,A25,A26,A27,A28,A29,A30,A31,A32,A33
C
      PLON=22640.D0*DSIN(A2)-4586.D0*DSIN(A2-2*A4)+2370*DSIN(2*A4)
     A+769*DSIN(2*A2)-668*DSIN(A8)-412*DSIN(2*A3)-212*DSIN(2*A2-2*A4)
     B-206*DSIN(A2-2*A4+A8)+192*DSIN(A2+2*A4)+165*DSIN(2*A4-A8)
     C+148*DSIN(A2-A8)-125*DSIN(A4)-110*DSIN(A2+A8)-55*DSIN(2*A3-2*A4)
     D-45*DSIN(A2+2*A3)+40*DSIN(A2-2*A3)-38*DSIN(A2-4*A4)
     E+36*DSIN(3*A2)-31*DSIN(2*A2-4*A4)+28*DSIN(A2-2*A4-A8)
     F-24*DSIN(2*A4+A8)+19*DSIN(A2-A4)+18*DSIN(A4+A8)
     G+15*DSIN(A2+2*A4-A8)+14*DSIN(2*A2+2*A4)+14*DSIN(4*A4)
     H-13*DSIN(3*A2-2*A4)-11*DSIN(A2+16*A7-18*A12)+10*DSIN(2*A2-A8)
     I+9*DSIN(A2-2*A3-2*A4)+9*DCOS(A2+16*A7-18*A12)
     J-9*DSIN(2*A2-2*A4+A8)-8*DSIN(A2+A4)+8*DSIN(2*A4-2*A8)
     K-8*DSIN(2*A2+A8)-7*DSIN(2*A8)-7*DSIN(A2-2*A4+2*A8)
     L+7*DSIN(A5)-6*DSIN(A2-2*A3+2*A4)-6*DSIN(2*A3+2*A4)
     M-4*DSIN(A2-4*A4+A8)+4*T*DCOS(A2+16*A7-18*A12)-4*DSIN(2*A2+2*A3)
     N+4*T*DSIN(A2+16*A7-18*A12)+3*DSIN(A2-3*A4)-3*DSIN(A2+2*A4+A8)
     O-3*DSIN(2*A2-4*A4+A8)+3*DSIN(A2-2*A8)+3*DSIN(A2-2*A4-2*A8)
     P-2*DSIN(2*A2-2*A4-A8)-2*DSIN(2*A3-2*A4+A8)+2*DSIN(A2+4*A4)
     Q+2*DSIN(4*A2)+2*DSIN(4*A4-A8)+2*DSIN(2*A2-A4)
C      DONE
      BETA=18461.D0*DSIN(A3)+1010*DSIN(A2+A3)+1000*DSIN(A2-A3)
     A-624*DSIN(A3-2*A4)-199*DSIN(A2-A3-2*A4)-167*DSIN(A2+A3-2*A4)
     B+117*DSIN(A3+2*A4)+62*DSIN(2*A2+A3)+33*DSIN(A2-A3+2*A4)
     C+32*DSIN(2*A2-A3)-30*DSIN(A3-2*A4+A8)-16*DSIN(2*A2+A3-2*A4)
     D+15*DSIN(A2+A3+2*A4)+12*DSIN(A3-2*A4-A8)-9*DSIN(A2-A3-2*A4+A8)
     E-8*DSIN(A3+A5)+8*DSIN(A3+2*A4-A8)-7*DSIN(A2+A3-2*A4+A8)
     F+7*DSIN(A2+A3-A8)-7*DSIN(A2+A3-4*A4)-6*DSIN(A3+A8)
     G-6*DSIN(3*A3)+6*DSIN(A2-A3-A8)-5*DSIN(A3+A4)-5*DSIN(A2+A3+A8)
     H-5*DSIN(A2-A3+A8)+5*DSIN(A3-A8)+5*DSIN(A3-A4)+4*DSIN(3*A2+A3)
     I-4*DSIN(A3-4*A4)-3*DSIN(A2-A3-4*A4)+3*DSIN(A2-3*A3)
     J-2*DSIN(2*A2-A3-4*A4)-2*DSIN(3*A3-2*A4)+2*DSIN(2*A2-A3+2*A4)
     K+2*DSIN(A2-A3+2*A4-A8)+2*DSIN(2*A2-A3-2*A4)+2*DSIN(3*A2-A3)
C      DONE
      RP=60.36298D0-3.27746D0*DCOS(A2)-0.57994D0*DCOS(A2-2*A4)
     A-0.46357d0*DCOS(2*A4)-0.08904d0*DCOS(2*A2)
     B+0.03865d0*DCOS(2*A2-2*A4)-0.03237d0*DCOS(2*A4-A8)
     &-0.02688d0*DCOS(A2+2*A4)-0.02358d0*DCOS(A2-2*A4+A8)
     &-0.02030d0*DCOS(A2-A8)+0.01719d0*DCOS(A4)
     D+0.01671d0*DCOS(A2+A8)+0.01247d0*DCOS(A2-2*A3)
     &+0.00704d0*DCOS(A8)
     E+0.00529d0*DCOS(2*A4+A8)-0.00524d0*DCOS(A2-4*A4)
     F+0.00398d0*DCOS(A2-2*A4-A8)-0.00366d0*DCOS(3*A2)
     G-0.00295d0*DCOS(2*A2-4*A4)-0.00263d0*DCOS(A4+A8)
     H+0.00249d0*DCOS(3*A2-2*A4)-0.00221d0*DCOS(A2+2*A4-A8)
     I+0.00185d0*DCOS(2*A3-2*A4)-0.00161d0*DCOS(2*A4-2*A8)
     J+0.00147d0*DCOS(A2+2*A3-2*A4)-0.00142d0*DCOS(4*A4)
     K+0.00139d0*DCOS(2*A2-2*A4+A8)-0.00118d0*DCOS(A2-4*A4+A8)
     L-0.00116d0*DCOS(2*A2+2*A4)-0.00110d0*DCOS(2*A2-A8)
C      DONE
      V=0.39558D0*DSIN(A3+A5)
     &+0.082d0*DSIN(A3)+0.03257d0*DSIN(A2-A3-A5)
     A+0.01092d0*DSIN(A2+A3+A5)+0.00666d0*DSIN(A2-A3)
     B-0.00644d0*DSIN(A2+A3-2*A4+A5)-0.00331d0*DSIN(A3-2*A4+A5)
     C-0.00304d0*DSIN(A3-2*A4)-0.0024d0*DSIN(A2-A3-2*A4-A5)
     D+0.00226d0*DSIN(A2+A3)-0.00108d0*DSIN(A2+A3-2*A4)
     E-0.00079d0*DSIN(A3-A5)+0.00078d0*DSIN(A3+2*A4+A5)
     F+0.00066d0*DSIN(A3+A5-A8)-0.00062d0*DSIN(A3+A5+A8)
     G-0.00050d0*DSIN(A2-A3-2*A4)+0.00045d0*DSIN(2*A2+A3+A5)
     H-0.00031d0*DSIN(2*A2+A3-2*A4+A5)-0.00027d0*DSIN(A2+A3-2*A4+A5+A8)
     I-0.00024d0*DSIN(A3-2*A4+A5+A8)-0.00021d0*T*DSIN(A3+A5)
     J+0.00018d0*DSIN(A3-A4+A5)+0.00016d0*DSIN(A3+2*A4)
     K+0.00016d0*DSIN(A2-A3-A5-A8)-0.00016d0*DSIN(2*A2-A3-A5)
     L-0.00015d0*DSIN(A3-2*A4+A8)-0.00012d0*DSIN(A2-A3-2*A4-A5+A8)
     M-0.00011d0*DSIN(A2-A3-A5+A8)+0.00009d0*DSIN(A2+A3+A5-A8)
     N+0.00009d0*DSIN(2*A2+A3)+0.00008d0*DSIN(2*A2-A3)
     O+0.00008d0*DSIN(A2+A3+2*A4+A5)-0.00008d0*DSIN(3*A3-2*A4+A5)
     P+0.00007d0*DSIN(A2-A3+2*A4)-0.00007d0*DSIN(2*A2-A3-2*A4-A5)
     Q-0.00007d0*DSIN(A2+A3+A5+A8)-0.00006d0*DSIN(A3+A4+A5)
     R+0.00006d0*DSIN(A3-2*A4-A8)+0.00006d0*DSIN(A2-A3+A5)
     S+0.00006d0*DSIN(A3+2*A4+A5-A8)-0.00005d0*DSIN(A2+A3-2*A4+A8)
     T-0.00004d0*DSIN(2*A2+A3-2*A4)+0.00004d0*DSIN(A2-3*A3-A5)
     U+0.00004d0*DSIN(A2-A3-A8)-0.00003d0*DSIN(A2-A3+A8)
     V+0.00003d0*DSIN(A3-A4)+0.00003d0*DSIN(A3-2*A4+A5-A8)
     W-0.00003d0*DSIN(A3-2*A4-A5)+0.00003d0*DSIN(A2+A3-2*A4+A5-A8)
     X+0.00003d0*DSIN(A3-A8)-0.00003d0*DSIN(A3-A4+A5-A8)
     Y-0.00002d0*DSIN(A2-A3-2*A4+A8)-0.00002d0*DSIN(A3+A8)
     Z+0.00002d0*DSIN(A2+A3-A4+A5)-0.00002d0*DSIN(A2+A3-A5)
     A+0.00002d0*DSIN(3*A2+A3+A5)-0.00002d0*DSIN(2*A2-A3-4*A4-A5)
     B+0.00002d0*DSIN(A2-A3-2*A4-A5-A8)-0.00002d0*T*DSIN(A2-A3-A5)
     C-0.00002d0*DSIN(A2-A3-4*A4-A5)-0.00002d0*DSIN(A2+A3-4*A4)
     D-0.00002d0*DSIN(2*A2-A3-2*A4)+0.00002d0*DSIN(A2+A3+2*A4)
     E+0.00002d0*DSIN(A2+A3-A8)
C     DONE
      U=1.d0-0.10828D0*DCOS(A2)
     &-0.0188d0*DCOS(A2-2*A4)-0.01479d0*DCOS(2*A4)
     A+0.00181d0*DCOS(2*A2-2*A4)-0.00147d0*DCOS(2*A2)
     B-0.00105d0*DCOS(2*A4-A8)-0.00075d0*DCOS(A2-2*A4+A8)
     C-0.00067d0*DCOS(A2-A8)+0.00057d0*DCOS(A4)+0.00055d0*DCOS(A2+A8)
     D-0.00046d0*DCOS(A2+2*A4)+0.00041d0*DCOS(A2-2*A3)
     &+0.00024d0*DCOS(A8)
     E+0.00017d0*DCOS(2*A4+A8)+0.00013d0*DCOS(A2-2*A4-A8)
     F-0.00010d0*DCOS(A2-4*A4)-0.00009d0*DCOS(A4+A8)
     G+0.00007d0*DCOS(2*A2-2*A4+A8)+0.00006d0*DCOS(3*A2-2*A4)
     H+0.00006d0*DCOS(2*A3-2*A4)-0.00005d0*DCOS(2*A4-2*A8)
     I-0.00005d0*DCOS(2*A2-4*A4)+0.00005d0*DCOS(A2+2*A3-2*A4)
     J-0.00005d0*DCOS(A2-A4)-0.00004d0*DCOS(A2+2*A4-A8)
     K-0.00004d0*DCOS(3*A2)-0.00003d0*DCOS(A2-4*A4+A8)
     L-0.00003d0*DCOS(2*A2-2*A3)-0.00003d0*DCOS(2*A3)
C     DONE
      W=0.10478D0*DSIN(A2)
     &-0.04105d0*DSIN(2*A3+2*A5)-0.0213d0*DSIN(A2-2*A4)
     A-0.01779d0*DSIN(2*A3+A5)+0.01774d0*DSIN(A5)+0.00987d0*DSIN(2*A4)
     B-0.00338d0*DSIN(A2-2*A3-2*A5)-0.00309d0*DSIN(A8)
     C-0.00190d0*DSIN(2*A3)-0.00144d0*DSIN(A2+A5)
     D-.001440d0*DSIN(A2-2*A3-A5)-0.00113d0*DSIN(A2+2*A3+2*A5)
     E-0.00094d0*DSIN(A2-2*A4+A8)-0.00092d0*DSIN(2*A2-2*A4)
     F+0.00071d0*DSIN(2*A4-A8)+0.00070d0*DSIN(2*A2)
     G+0.00067d0*DSIN(A2+2*A3-2*A4+2*A5)+0.00066d0*DSIN(2*A3-2*A4+A5)
     H-0.00066d0*DSIN(2*A4+A5)+0.00061d0*DSIN(A2-A8)-0.00058d0*DSIN(A4)
     I-0.00049d0*DSIN(A2+2*A3+A5)-0.00049d0*DSIN(A2-A5)
     J-0.00042d0*DSIN(A2+A8)+0.00034d0*DSIN(2*A3-2*A4+2*A5)
     K-0.00026d0*DSIN(2*A3-2*A4)+0.00025d0*DSIN(A2-2*A3-2*A4-2*A5)
     L+0.00024d0*DSIN(A2-2*A3)+0.00023d0*DSIN(A2+2*A3-2*A4+A5)
     M+0.00023d0*DSIN(A2-2*A4-A5)+0.00019d0*DSIN(A2+2*A4)
     N+0.00012d0*DSIN(A2-2*A4-A8)+0.00011d0*DSIN(A2-2*A4+A5)
     O+0.00011d0*DSIN(A2-2*A3-2*A4-A5)-0.00010d0*DSIN(2*A4+A8)
     P+0.00009d0*DSIN(A2-A4)+0.00008d0*DSIN(A4+A8)
     Q-0.00008d0*DSIN(2*A3+2*A4+2*A5)-0.00008d0*DSIN(2*A5)
     R-0.00007d0*DSIN(2*A3+2*A5-A8)+0.00006d0*DSIN(2*A3+2*A5+A8)
     S-0.00005d0*DSIN(A2+2*A3)+0.00005d0*DSIN(3*A2)
     T-0.00005d0*DSIN(A2+16*A7-18*A12)-0.00005d0*DSIN(2*A2+2*A3+2*A5)
     U+0.00004d0*T*DSIN(2*A3+2*A5)+0.00004d0*DCOS(A2+16*A7-18*A12)
     V-0.00004d0*DSIN(A2-2*A3+2*A4)-0.00004d0*DSIN(A2-4*A4)
     X-0.00004d0*DSIN(3*A2-2*A4)-0.00004d0*DSIN(2*A3+2*A4+A5)
     Y-0.00004d0*DSIN(2*A4-A5)-0.00003d0*DSIN(2*A8)
     Z-0.00003d0*DSIN(A2-2*A4+2*A8)+0.00003d0*DSIN(2*A3-2*A4+A5+A8)
     A-0.00003d0*DSIN(2*A4+A5-A8)+0.00003d0*DSIN(2*A2+2*A3-2*A4+2*A5)
     B+0.00003d0*DSIN(2*A4-2*A8)-0.00003d0*DSIN(2*A2-2*A4+A8)
     C+0.00003d0*DSIN(A2+2*A3-2*A4+2*A5+A8)-0.00003d0*DSIN(2*A2-4*A4)
     D+0.00002d0*DSIN(2*A3-2*A4+2*A5+A8)-0.00002d0*DSIN(2*A2+2*A3+A5)
     E-0.00002d0*DSIN(2*A2-A5)+0.00002d0*T*DCOS(A2+16*A7-18*A12)
     F+0.00002d0*DSIN(4*A4)-0.00002d0*DSIN(2*A3-A4+2*A5)
     G-0.00002d0*DSIN(A2+2*A3-2*A4)-0.00002d0*DSIN(2*A2+A5)
     H-0.00002d0*DSIN(2*A2-2*A3-A5)+0.00002d0*DSIN(A2+2*A4-A8)
     I+0.00002d0*DSIN(2*A2-A8)-0.00002d0*DSIN(A2-4*A4+A8)
     J+0.00002d0*T*DSIN(A2+16*A7-18*A12)-0.00002d0*DSIN(A2-2*A3-2*A5-A8)
     K+0.00002d0*DSIN(2*A2-2*A3-2*A5)-0.00002d0*DSIN(A2+2*A4+A5)
     L-0.00002d0*DSIN(A2-2*A3+2*A4-A5)
C      DONE
      RETURN
      END