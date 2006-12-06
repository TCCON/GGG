c  VMOV:  moves the contents of a vector "vecin" to the vector "vecout".
c
c     calling sequence:
c     call vmov ( vecin,incrin,vecout,incrout,nele )
c
c     input parameters:
c     vecin : real*4 vector to move FROM
c     incrin : integer array index increment for "vecin"
c     incrout : integer array index increment for "vecout"
c     nele : integer number of elements on which to perform
c            the operation
c
c     output parameters:
c     vecout : real*4 vector to contain the values from "vecin"
c     vecout(1+m*incrout)=vecin(1+m*incrin)  for m=0 to nele-1
c
      subroutine vmov(vecin,incrin,vecout,incrout,nele)
      implicit none
      integer i,j,k,incrin,incrout,nele
      real*4 vecin(*),vecout(*)
      j=1
      k=1
      do 100 i=1,nele
        vecout(k)=vecin(j)
        j=j+incrin
        k=k+incrout
100   continue
      return
      end
c==========================================================================
c  VSWAP:  swaps the contents of two vectors.
c
c     calling sequence:
c     call vswap ( v1,inc1,v2,inc2,nele )
c
c     input parameters:
c     v1    : real*4 vector 
c     inc1  : integer array index increment for V1
c     v2    : real*4 vector 
c     inc2  : integer array index increment for V2
c     nele : integer number of elements to swap
c
c     output parameters:
c     v1    : real*4 vector 
c     v2    : real*4 vector 
c
c     v2(1+m*inc2)=v1(1+m*inc1)  for m=0 to nele-1
c     v1(1+m*inc1)=v2(1+m*inc2)  for m=0 to nele-1
c
c  Note:  VSWAP may be used to reverse a vector in place by calling
c     CALL VSWAP(V1,INC1,V1(1+INC1*(NELE-1)),-INC1,NELE/2)
c
      subroutine vswap(v1,inc1,v2,inc2,nele)
      implicit none
      integer i,j1,j2,inc1,inc2,nele
      real*4 v1(*),v2(*),dum
      j1=1
      j2=1
      do 100 i=1,nele
         dum=v1(j1)
         v1(j1)=v2(j2)
         v2(j2)=dum
         j1=j1+inc1
         j2=j2+inc2
100   continue
      return
      end
c==================================================================
c  VADD: subroutine to add two n x 1 vectors "v1" and "v2" together
c     and return the resulting n x 1 vector in "vout".
c
c     input parameters:
c     v1 : real*4(n) vector to add
c     v2 : real*4(n) vector to add
c     n : integer dimension of all vectors
c
c     output parameters:
c     vout : real*4(n) vector containing sum of v1 and v2
c
      subroutine vadd(v1,i1,v2,i2,vout,iout,nele)
      implicit none
      integer k1,k2,kout,i,i1,i2,iout,nele
      real*4 v1(1+i1*(nele-1)),v2(1+i2*(nele-1)),vout(1+iout*(nele-1))
      k1=1
      k2=1
      kout=1
      do 100 i=1,nele
        vout(kout) = v1(k1) + v2(k2)
        k1=k1+i1
        k2=k2+i2
        kout=kout+iout
100   continue
      return
      end
c====================================================================
c  VMUL:  multiplies the elements of a vector "v1" by the corresponding
c     elements of the vector "v2", and return the resulting vector in "v3".
c 
c     calling sequence:
c     call vmul (v1,i1,v2,i2,v3,i3,nele)
c
c     input parameters:
c     v1, v2 : real*4 vectors to multiply
c     i1, i2 : integer array index increments for "v1" and "v2" respectively
c     i3 : integer array index increment for "v3"
c     nele : integer number of elements on which to perform the operation
c
c     Output parameters:
c     v3(1+m*i3)=v1(1+m*i1)*v2(1+m*i2) for m=0 to nele-1
c
      subroutine vmul(v1,i1,v2,i2,vec3,i3,nele)
      implicit none
      integer i,j,k,l,i1,i2,i3,nele
      real*4 v1(1+i1*(nele-1)),v2(1+i2*(nele-1)),vec3(1+i3*(nele-1))
      j=1
      k=1
      l=1
      do 100 i=1,nele
        vec3(l)=v1(j)*v2(k)
        j=j+i1
        k=k+i2
        l=l+i3
100   continue
      return
      end
c====================================================================
c  VDIV:  divides the elements of a vector "v1" by the corresponding
c     elements of the vector "v2", and return the resulting vector in "v3".
c 
c     calling sequence:
c     call vmul (v1,i1,v2,i2,v3,i3,nele)
c
c     input parameters:
c     v1, v2 : real*4 vectors to divide
c     i1, i2 : integer array index increments for "v1" and "v2" respectively
c     i3 : integer array index increment for "v3"
c     nele : integer number of elements on which to perform the operation
c
c     Output parameters:
c     v3(1+m*i3)=v1(1+m*i1)/v2(1+m*i2) for m=0 to nele-1
c
      subroutine vdiv(v1,i1,v2,i2,vec3,i3,nele)
      implicit none
      integer i,j,k,l,i1,i2,i3,nele
      real*4 v1(1+i1*(nele-1)),v2(1+i2*(nele-1)),vec3(1+i3*(nele-1))
      j=1
      k=1
      l=1
      do 100 i=1,nele
        vec3(l)=v1(j)/v2(k)
        j=j+i1
        k=k+i2
        l=l+i3
100   continue
      return
      end
c=====================================================================
c  VRAMP:  Generates a ramp or arithmetic progression of floating point
c  values whose extrema differ by 1.0 and in which point (NMP+1)/2 is zero
c 
c     input parameters:
c     nele : number of elements on which to perform the operation
c
c     Output parameters:
c     vec,: real*4 array containing ramp values
c     inc,: integer array address increment for "vec" 
c     vec(i)=-0.5+float(i-1)/(nele-1) for i=1,nele
c
      subroutine vramp(vec,inc,nele)
      implicit none
      integer i,k,nele,nop,inc
      real*4 vec(1+inc*(nele-1))
      nop=(nele+1)/2
      k=1
      do 100 i=1,nele
        vec(k)=+float(i-nop)/((nele-1)+1.e-37) ! adding 1e+38 stops NaN when nele=1
        k=k+inc
100   continue
      return
      end
c===========================================================
c   VSUB:  subtracts two n x 1 vectors "v1" and "v2"
c     and return the resulting n x 1 vector in "vout".
c
c     input parameters:
c     v1 : real*4(nele) vector to add
c     v2 : real*4(nele) vector to add
c     nele : integer dimension of all vectors
c
c     output parameters:
c     vout : real*4(nele) vector containing V1-V2
c
      subroutine vsub(v1,i1,v2,i2,vout,iout,nele)
      implicit none
      integer k1,k2,kout,i,i1,i2,iout,nele
      real*4 v1(1+i1*(nele-1)),v2(1+i2*(nele-1)),vout(1+iout*(nele-1))
      k1=1
      k2=1
      kout=1
      do 100 i=1,nele
        vout(kout) = v1(k1) - v2(k2)
        k1=k1+i1
        k2=k2+i2
        kout=kout+iout
100   continue
      return
      end
c===================================================================
c  VDOT:  performs the dot-product of a vector "v1" with the vector "v2",
c     and returns the resulting scalar in "prod".
c
c     calling sequence:
c     call vdot (v1,incr1,v2,incr2,prod,nele)
c
c     input parameters:
c     v1, v2 : real*4 vectors to take dot product of
c     incr1, incr2 : integer array index increments for "v1" and "v2",
c                    respectively
c     nele : integer number of elements on which to perform
c            the operation
c
c     output parameters:
c     prod : real*4 dot product result
c     prod = sum_over_i {v1(1+i*incr1)*v2(1+i*incr2)} for i=0 to nele-1
c
      subroutine vdot(v1,incr1,v2,incr2,prod,nele)
      implicit none
      integer i,j,k,incr1,incr2,nele
      real*8 dp
      real*4 prod,v1(*),v2(*)
      j=1
      k=1
      dp=0.0d0
      do 100 i=1,nele
        dp = dp + dprod(v1(j),v2(k))
        j=j+incr1
        k=k+incr2
100   continue
      prod=sngl(dp)
      return
      end
c================================================================
c  VEXP: exponentiates the vector VECIN(nele)
c
c   input parameters:
c       vecin :     real*4     vector to exponentiate
c       incrin  :   integer    array index increment for "vecin"
c       incrout :   integer    array index increment for "vecout"
c       nele :      integer    number of times to perform the operation
c
c   output parameters:
c       vecout : real*4 vector to contain the exponentiated values
c
      subroutine vexp(vecin,incrin,vecout,incrout,nele)
      implicit none
      integer i,j,k,incrin,incrout,nele
      real*4 vecin(*),vecout(*),x
      j=1
      k=1
      do i=1,nele
        x=vecin(j)
        if(x.gt. 0.0) then
           vecout(k)=1.+x*(1.+0.5*x) ! quad expansion prevents overflow for large x
        elseif(x.lt.-80.0) then
           vecout(k)=0.0
        else
           vecout(k)=exp(x)
        endif
        j=j+incrin
        k=k+incrout
      end do
      return
      end
c==============================================================
c  VLOGe: takes the natural logarithm of the vector VECIN(nele)
c
c   input parameters:
c       vecin :     real*4     vector to exponentiate
c       incrin  :   integer    array index increment for "vecin"
c       incrout :   integer    array index increment for "vecout"
c       nele :      integer    number of times to perform the operation
c
c   output parameters:
c       vecout : real*4 vector to contain the exponentiated values
c
      subroutine vloge(vecin,incrin,vecout,incrout,nele)
      implicit none
      integer i,j,k,incrin,incrout,nele
      real*4 vecin(*),vecout(*),x
      j=1
      k=1
      do i=1,nele
        x=vecin(j)
        vecout(k)=log(x)
        j=j+incrin
        k=k+incrout
      end do
      return
      end
c==============================================================
c  VSIN:  calculates the sine of the vector VECIN(nele)
c
c   input parameters:
c       vecin :     real*4     vector to exponentiate
c       incrin  :   integer    array index increment for "vecin"
c       incrout :   integer    array index increment for "vecout"
c       nele :      integer    number of times to perform the operation
c
c   output parameter:
c       vecout : real*4 vector to contain the sined values
c
      subroutine vsin(vecin,incrin,vecout,incrout,nele)
      implicit none
      integer i,j,k,incrin,incrout,nele
      real*4 vecin(*),vecout(*)
      j=1
      k=1
      do i=1,nele
        vecout(k)=sin(vecin(j))
        j=j+incrin
        k=k+incrout
      end do
      return
      end
c==============================================================
c  VCOS:  calculates the cosine of the vector VECIN(nele)
c
c   input parameters:
c       vecin :     real*4     vector to exponentiate
c       incrin  :   integer    array index increment for "vecin"
c       incrout :   integer    array index increment for "vecout"
c       nele :      integer    number of times to perform the operation
c
c   output parameter:
c       vecout : real*4 vector to contain the cosined values
c
      subroutine vcos(vecin,incrin,vecout,incrout,nele)
      implicit none
      integer i,j,k,incrin,incrout,nele
      real*4 vecin(*),vecout(*)
      j=1
      k=1
      do i=1,nele
        vecout(k)=cos(vecin(j))
        j=j+incrin
        k=k+incrout
      end do
      return
      end
c==============================================================
c  VSQRT: takes the square root of the NELE elements of the vector VECIN
c     calling sequence:
c     call vsqrt ( vecin,incrin,vecout,incrout,nele )
c
c     input parameters:
c     vecin : real*4 vector to exponentiate
c     incrin : integer array index increment for "vecin"
c     incrout : integer array index increment for "vecout"
c     nele : integer number of elements on which to perform
c            the operation
c
c     output parameters:
c     vecout : real*4 vector to contain the exponentiated values
c     vecout(1+m*incrout)=sqrt(vecin(1+m*incrin))  for m=0 to nele-1
c
c  For +ve X (-ve absorber amounts) use a simple quadratic expansion.
c  This prevents overflow for very large X.
c
      subroutine vsqrt(vecin,incrin,vecout,incrout,nele)
      implicit none
      integer i,j,k,incrin,incrout,nele
      real*4 vecin(*),vecout(*)
      j=1
      k=1
      do 100 i=1,nele
        vecout(k)=sign(sqrt(abs(vecin(j))),vecin(j))
        j=j+incrin
        k=k+incrout
100   continue
      return
      end
c===============================================================
c  VSMA: Multiplies vector V1 by the scalar SCAL, then adds the
c  result to vector V2.
c     calling sequence:
c     call vsma (v1,i1,scal,v2,i2,vout,iout,nele)
c
c     input parameters:
c     v1 : real*4 vector
c     i1 : integer array index iement for "v1"
c     scal : real*4 scalar to multiply the input vector by
c     v2 : real*4 vector to add to the above product
c     i2 : real*4 array index iement for "v2"
c     iout : integer array index iement for "vecout"
c     nele : integer number of elements on which to perform
c            the operation
c
c     output parameters:
c     vout : real*4 vector containing V1*SCAL+V2
c     vout(1+m*iout)=(v1(1+m*i1)*scal)+v2(1+m*i2) for m=0 to nele-1
c
      subroutine vsma(v1,i1,scal,v2,i2,vout,iout,nele)
      implicit none
      integer i,j,k,l,i1,i2,iout,nele
      real*4 v1(1+i1*(nele-1)),v2(1+i2*(nele-1)),vout(1+iout*(nele-1)),
     &scal
      j=1
      k=1
      l=1
      do 100 i=1,nele
        vout(l)=(v1(j)*scal)+v2(k)
        j=j+i1
        k=k+i2
        l=l+iout
100   continue
      return
      end
c===============================================================
c  VMAX:  finds largest element of vector.
c     v1 : real*4 vector
c     i1 : integer array index increment for "v1"
c     nele : integer number of elements on which to perform search
c
c     output parameters:
c     vmax : real*4 scalar containing largest value
c
      real*4 function vmax(v1,i1,nele)
      integer j,k,i1,nele
      real*4 v1(*)
      j=1
      vmax=v1(j)
      do k=2,nele
         j=j+i1
         if(v1(j).gt.vmax) vmax=v1(j)
      end do
      return
      end
c===============================================================
c  VMIN:  finds largest element of vector.
c     v1 : real*4 vector
c     i1 : integer array index increment for "v1"
c     nele : integer number of elements on which to perform search
c
c     output parameters:
c     vmin : real*4 scalar containing largest value
c
      real*4 function vmin(v1,i1,nele)
      integer j,k,i1,nele
      real*4 v1(*)
      j=1
      vmin=v1(j)
      do k=2,nele
         j=j+i1
         if(v1(j).lt.vmin) vmin=v1(j)
      end do
      return
      end
