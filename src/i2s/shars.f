      subroutine shars(kmode,nip,nhw,shbar,yin,youtshars)
c  Shifts the odd/even interferogram points with respect
c  to each other and re-samples the interferogram.
c  This removes the effects of laser Sampling Error (LSE).
c
c  Inputs:
c       KMODE   I*4  (-1, 0, +1)  
c         NIP   I*4   Number of igram points
c         NHW   I*4   half-width of re-sampling operator
c       SHBAR   R*4   Shift to be applied (fraction of point spacing)
c    YIN(NIP)   R*4   Raw Interferogram
c
c  Outputs:
c    YOUT(NIP)  R*4   Re-sampled igram
c
c  Theoretical basis:
c
c  Lagrange interpolation with unequally-spaced abscissae.
c
c     Y(X)= SUM_over_j [Yj.Pj]
c     Pj = PRODUCT_over_k [ (X-Xk)/(Xj-Xk) ]     k.ne.j
c  where Yj are the igram values measured at Xj
c  and Y(X) is the interpolated igram at any new X value
c
c  Assumes that we want the new X values to be equally spaced
c    Xi = dx.i
c  but that the old Xj had a periodic sampling error
c    Xj = dx.(j+shbar.(abs(kmode)+(-1)^j))
c    Xk = dx.(k+shbar.(abs(kmode)+(-1)^k))
c  where shbar is the sampling error expressed as a fraction
c  of the sampling interval, dx=lambda/2.
c
c  Implementation Notes:
c
c  Although an odd-even shift of SHBAR could be corrected in an
c  infinite number of ways (e.g. shift the odd points by 0.2*SHBAR
c  and the evens by -0.8*SHBAR), there are three special cases:
c    KMODE=-1: Shift odd points by  -SHBAR  and even points by  0
c    KMODE= 0: Shift odd points by -SHBAR/2 and even points by +SHBAR/2
c    KMODE=+1: Shift odd points by    0    and even points by  +SHBAR
c  These cases are equivalent as far as the resulting phase-corrected
c  spectra are concerned, but the igrams will be different in each case.
c
c  These three cases are special in the sense that only one re-sampling
c  operator is needed. For the KMODE=0 case, the -SHBAR/2 resampling
c  operator is the mirror image of the +SHBAR/2 operator. And for the
c  KMODE= +/- 1 cases, the unshifted points are simply left unchanged.
c  Any other case would require two operators (one for shifting 
c  the odds by +0.2*SHBAR and another for the evens by -0.8*SHBAR).
c
c  Since the numerical error scales as SHBAR^2, then KMODE=0 method
c  is slightly more accurate than the KMOD=+/-1 methods (for a given
c  operator size). But the KMODE=+/-1 methods are faster, requiring
c  only half the points be re-sampled.
c
c  Only the igram points from 1+nhw to nip-nhw are actually re-sampled.
c  The end points are left alone. Since NIP >> NHW, this is a negligible
c  source of error. In fact, the end-points on the short side of ZPD
c  are not even used.
c
      implicit none
      integer i,m,nhw,nip,kmode,kflip,iflip,imf
      real*4 yin(nip),youtshars(nip),shbar
      real*4 oper(-nhw:nhw)
      real*8 tyy

c  Pre-Compute shifting/re-sampling operator (OPER),
c  in fact, its difference from a delta-function.
      call glpio(nhw,kmode,shbar,oper)
c
c  The loop below performs convolution of YIN with OPER.
c  Small values at ends of OPER are done first,
c  large values last, for slightly improved accuracy.
c  For KMODE=1, odd points of output vector are skipped.
      iflip=1
      if(kmode.eq.-1) iflip=-1
      kflip=(-1)**(nhw-1)  ! = (-1)**(i-nhw)
      do i=1+nhw,nip-nhw
         if(kmode.eq.0) iflip=-iflip
         if(kmode*kflip.ge.0) then
            tyy=0.0d0
            do m=nhw,1,-1
               imf=m*iflip
               tyy = tyy + yin(i-m)*oper(-imf) + yin(i+m)*oper(imf)
            end do
            tyy = tyy + yin(i)*(1.0d0+oper(0)) ! Add back 1.0 to OPER(0)
            youtshars(i)=sngl(tyy)
         else
            youtshars(i)=yin(i)  ! skip point
         endif
         kflip=-kflip
      end do
      return
      end

      subroutine glpio(nhw,kmode,shbar,oper)
c  Generates an operator OPER which, when convolved
c  with the interferogram, shifts the points by +/- SHBAR
c  (or 0 & 2*SHBAR) and re-samples at the new point locations.
c
c  Inputs:
c      NHW    I*4  Half-width of operator
c    KMODE    I*4  -1,0,+1
c    SHBAR    R*4  The shift to be imparted (fraction of point spacing)
c
c  Output:
c    OPER(-NHW:NHW) R*4  Convolution operator
c
c  Theoretical Basis:
c
c  Lagrange interpolation with unequally-spaced abscissae.
c
c     Y(X)= SUM_over_j [Yj.Pj]
c     Pj = PRODUCT_over_k [ (X-Xk)/(Xj-Xk) ]     k.ne.j
c  where Yj are the igram values measured at Xj
c  and Y(X) is the interpolated igram at any new X value
c
c  Assumes that we want the new X values to be equally spaced
c     Xi = dx.i
c  but that the old Xj had a periodic sampling error
c     Xj = dx.(j+shbar.(abs(kmode)+(-1)^j))
c     Xk = dx.(k+shbar.(abs(kmode)+(-1)^k))
c  where shbar is the sampling error expressed as a fraction
c  of the sampling interval, dx=lambda/2.
c
c  Implementation Notes:
c
c    KMODE=-1: Shift odd points by  -SHBAR  and even points by  0
c    KMODE= 0: Shift odd points by -SHBAR/2 and even points by +SHBAR/2
c    KMODE=+1: Shift odd points by      0   and even points by +SHBAR
c
c  For small SHBAR values, the operator is close to a delta-function:
c  with a value of 1.0 at the center and very small everywhere else.
c  So to improve the numerical accuracy, 1.0 is subtracted from the
c  middle point, and must be added back when the operator is used.
c  This trick allows precision of better than 10^-6 to be achieved
c  with a single-precision operator.

      implicit none
      integer j,k,nhw,kflip0,kflip,jflip,kmode,iabkm
      real*4 oper(-nhw:nhw),shbar
      real*8 ww,xk,xj

      iabkm=iabs(kmode)
      kflip0=(-1)**nhw
      kflip=kflip0
      do k=-nhw,+nhw
         ww=1.0d0
         xk=k+0.5*shbar*(kflip+iabkm)
         jflip=kflip0
         do j=-nhw,+nhw
            if(j.ne.k) then
               xj=j+0.5*shbar*(jflip+iabkm)
               ww=ww*xj/(xj-xk)
c         write(34,'(i4,f9.5,2i8,2f9.5,e13.5)')kmode,shbar,k,j,xk,xj,ww
            endif
            jflip=-jflip
         end do
c         write(33,*)kmode,shbar,k,ww
         if(k.eq.0) ww=ww-1.0d0  ! subtract 1 from middle point of OPER
         oper(k)=sngl(ww)
         kflip=-kflip
      end do
c      write(*,*)'oper=',oper(-1),oper(0),oper(1)
      return
      end
