      subroutine julian(y,m,d,jd)
C  Subroutine converts Gregorian calendar date
C  to Julian Day number at Greenwich Mean Noon
      implicit none
      INTEGER*4 Y,M,D,JD
C******************************************************
C        Input:
C     Y            Integral calendar year
C     M            Integral calendar month
C     D            Integral calendar day
C
C        Output:
C     JD           Integral Julian Day number
C******************************************************
      jd=367*y-7*(y+(m+9)/12)/4-3*((y+(m-9)/7)/100+1)/4
     $+275*m/9+d+1721029
      return
      end

