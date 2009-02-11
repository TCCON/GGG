      subroutine endian(iend)
c  Determines the "endian" of the host computer
c     Inputs:
c              None
c     Outputs:
c              iend   i*4
c
c  iend=0 on big-endian CPU (e.g. Sun); iend=1 on little-endian CPU (e.g. PC)
c
      implicit none
      integer*2 i2
      integer*4 i4,iend
      equivalence (i2,i4)
c
      i4=1
      iend=i2  !   iend = 0 on SUN;    iend = 1 on PC
      return
      end
