      subroutine rbyte(b,lenw,nword)
c  Performs an in-place byte reversal of array/string B containing
c  NWORD consecutive data words each of length LENW bytes.
c  In the calling program C can be of any data type ( I*2, I*4, R*4,
c  R*8). It doesn't matter because only the starting
c  address of B is actually passed to this subroutine.
c
c  INPUTS:
c      B        Data array to be reversed
c      LENW     Word Length (bytes)
c      NWORD    Number of consecutive words to be reversed
c
c   OUTPUT:
c      B        Data array/string (byte reversed)
c
c  Note that although byte array b(lenw,nword) is 2-dimensional
c  here, this is just to simplify the indexing.   In the main
c  program, it might only be one-dimensional or a scalar. 
c
      implicit none
      integer*4 lenw,nword,iword,ibyte,kbyte
      byte b(lenw,nword),temporary

      do iword=1,nword
         do ibyte=1,lenw/2
            kbyte=lenw+1-ibyte
            temporary=b(kbyte,iword)
            b(kbyte,iword)=b(ibyte,iword)
            b(ibyte,iword)=temporary
         end do              !  ibyte=1,lenw/2
      end do               !  iword=1,nword
      return
      end
