      subroutine byterev(c,nbytes)
c  Performs an in-place byte reversal of variable C (length NBYTE)
c  In the calling program C can be I*2, I*4, R*4, R*8 or anything else
      integer*4 nbytes,k,j
      character c*(*),chtemp*1
      do k=1,nbytes/2
        j=nbytes+1-k
        chtemp=c(k:k)
        c(k:k)=c(j:j)
        c(j:j)=chtemp
      end do
      return
      end
