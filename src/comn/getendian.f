      subroutine getendian(iend)
c  Determines the "endian" of the host computer
c     Inputs:
c              None
c
c     Outputs:
c              iend   i*4
c
c     iend=+1 on big-endian CPU (e.g. Sun)
c     iend=-1 on little-endian CPU (e.g. PC)
c
      implicit none
      integer*2 i2
      integer*4 i4,iend
      equivalence (i2,i4)
c
c 131071 is the number whose first 2 bytes = +1 and whose second 2 bytes = -1.
      i4=131071   ! =  2**17-1
      iend=i2
      return
      end
