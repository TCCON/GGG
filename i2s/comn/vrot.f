      subroutine vrot(y,kshift,npts)
c  Performs an in-place, circular shift/rotation of the contents of Y(NPTS)
c  such that  Y(i) = Y(i+kshift)
c  If i+kshift > npts   then NPTS is subtracted from the index
c  If i+kshift < 1     then NPTS is added to the index
c
c  Inputs:
c      Y(NPTS)  R*4  Vector of data values
c      KSHIFT   I*4  Shift (right=-ve; left=+ve)
c      NPTS     I*4  Length of vector Y
c
c  Output:
c      Y(NPTS)  R*4  Vector of shifted data values
c
      implicit none
      integer*4 kshift,npts,gcd,kg,i,j,jj,jwas
      real*4 y(npts),ystart

      kg=iabs(gcd(npts,kshift))
c      write(*,*)npts,kshift,kg

      do i=1,kg
         jwas=i
         jj=i
         ystart=y(jwas)
         do j=2,npts/kg
            jj=mod(jj+kshift-1,npts)+1
            if(jj.le.0) jj=jj+npts
            y(jwas)=y(jj)
            jwas=jj
         end do
         y(jj)=ystart
      end do
      return
      end

      function gcd(a,b)
c  Greatest Common Divisor
      integer a1,b1,a,b,t,gcd
      a1=a
      b1=b
      do while (b1.ne.0)
         t=b1
         b1=mod(a1,b1)
         a1=t
      end do
      gcd=a1
      return
      end
