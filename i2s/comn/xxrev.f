      subroutine i2rev(i2var)
c
c  Performs an in-place byte reversal of I*2 variable i2var.
c
      implicit none
      integer*2 i2var,i2loc
      integer*1 b2loc(2),btemp
      equivalence (i2loc,b2loc)

      i2loc=i2var

      btemp=b2loc(1)
      b2loc(1)=b2loc(2)
      b2loc(2)=btemp

      i2var=i2loc
      return
      end
c ===================================================================
      subroutine i4rev(i4var)
c
c  Performs an in-place byte reversal of I*4 variable i4var.
c
      implicit none
      integer*4 i4var,i4loc
      integer*1 b4loc(4),btemp
      equivalence (i4loc,b4loc)

      i4loc=i4var

      btemp=b4loc(1)
      b4loc(1)=b4loc(4)
      b4loc(4)=btemp

      btemp=b4loc(2)
      b4loc(2)=b4loc(3)
      b4loc(3)=btemp

      i4var=i4loc
      return
      end
c ===================================================================
      subroutine r4rev(r4var)
c
c  Performs an in-place byte reversal of R*4 variable r4var.
c
      implicit none
      real*4 r4var,r4loc
      integer*1 b4loc(4),btemp
      equivalence (r4loc,b4loc)

      r4loc=r4var

      btemp=b4loc(1)
      b4loc(1)=b4loc(4)
      b4loc(4)=btemp

      btemp=b4loc(2)
      b4loc(2)=b4loc(3)
      b4loc(3)=btemp

      r4var=r4loc
      return
      end
c ===================================================================
      subroutine r8rev(r8var)
c
c  Performs an in-place byte reversal of R*8 variable r8var.
c
      implicit none
      real*8 r8var,r8loc
      integer*1 b8loc(8),btemp
      equivalence (r8loc,b8loc)

      r8loc=r8var

      btemp=b8loc(1)
      b8loc(1)=b8loc(8)
      b8loc(8)=btemp

      btemp=b8loc(2)
      b8loc(2)=b8loc(7)
      b8loc(7)=btemp

      btemp=b8loc(3)
      b8loc(3)=b8loc(6)
      b8loc(6)=btemp

      btemp=b8loc(4)
      b8loc(4)=b8loc(5)
      b8loc(5)=btemp

      r8var=r8loc
      return
      end
