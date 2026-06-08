    integer function julian(y,m,d)
!   Function returns Julian Day number at Greenwich Mean Noon
!   given Gregorian calendar date

    implicit none
    INTEGER(4):: Y,M,D
!******************************************************
!        Input:
!     Y            Integral calendar year
!     M            Integral calendar month
!     D            Integral calendar day
!
!        Output:
!     JD           Integral Julian Day number
!******************************************************
    julian=367*y-7*(y+(m+9)/12)/4-3*((y+(m-9)/7)/100+1)/4+275*m/9+d+1721029
    return
    end
