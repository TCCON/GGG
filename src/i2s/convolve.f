      subroutine convolve(method,nfftpc,ac_igram,workspace,lpco,pco)
c  Convolves of the raw interferogram, AC_IGRAM, of length NFFTPC
c  with the operator, PCO, of length LPCO.
c
c  The resulting (phase-corrected) interferogram is written in-place 
c  in AC_IGRAM, with the fringes lined up with those in the raw igram.
c
c  The CONVOLVE subroutine doesn't know or care where ZPD is located.
c
c  INPUTS:
c
c     METHOD            1=slow, direct, convolution
c                       2=fast, FFT, convolution
c
c     NFFTPC            Size of FFT
c
c     AC_IGRAM(NFFTPC)  Raw interferogram
c
c     WORKSPACE(NFFTPC) Used to FFT/interpolate PCO
c
c     LPCO              Length of Phase Correction Operator
c
c     PCO(LPCO)         FFT-packed Phase Correction Operator
c
c
c
c  OUTPUTS:
c
c     AC_IGRAM          Phase-corrected interferogram containing:
c        (1:LPCO/2-1)             Garbage (spectral contamination)
c        (LPCO/2:NFFTPC-LPCO/2)   Corrected igram (NFFT-LPCO+1 points)
c        (NFFT-LPCO/2+1,NFFTPC)   Garbage (spectral contamination)
c
c  The subroutine doesn't care about the size of the interferogram
c  that partially or fully occupies AC_IGRAM. The full NFFTPC points
c  get convolved, whether they are real igram points or just padding.
c  So the calling program must ensure, if the true size of AC_IGRAM
c  is less than NFFTPC, that the extra elements are all zeroed.
c  Otherwise they will get sucked into the convolution and cause damage.
c
c  Convolution is defined as     f*g(k) = SUM_i [f(i)*g(k-i])
c  Note that the operator, g, must be reversed in the sense that it is
c  a function of -i rather than +i.
c
c  PCO is in a FFT-packed format, meaning that the largest/center
c  point of the operator, originally at PCO(1+LPCO/2), is at PCO(1).
c  The point before this, PCO(LPCO/2), is now at PCO(LPCO)

c  Note that there are two possible ways of doing the convolution:
c   1) By convolution in the time/opd domain, which is slow but simple
c   2) By complex multiplication in the spectral domain, which is fast
c      but messy because it requires forward and reverse FFT's, and
c      requires extra memory to accomodate the interpolated PCO.
c
c  Note that on exit, AC_IGRAM can contain garbage outside the range
c  LPCO/2 to NFFTPC-LPCO/2. In the case of slow convolution (method=1)
c  this is because array elements outside this range are never set.
c  In the case of fast convolution (method=2) these elements are set
c  but can contain garbage from so-called wrap-around effects.
c
c  Upon return, the vector ac_igram has valid values in:
c     ac_igram(l2:nfftpc-l2),
c  a total of nfftpc-l2-l2+1 = nfftpc-lpco+1 good points.
c  Possibly invalid values are left at the beginning and end:
c     ac_igram(1:l2-1)              (l2-1 points)
c     ac_igram(nfftpc-l2+1:nfftpc)  (l2 points)
c  So for lpco=4, l2=2, the first point and last two points are invalid.
c  So for lpco=8, l2=4, the first three and last four points are invalid.
c
c  Why are more points lost at the end of the vector than at the beginning?
c  Because PCO is an even length operator with its center point at LPCO/2+1.
c  PCO is reversed for convolution, which moves the center point to LPCO/2.
c  The first point of AC_IGRAM that can be sampled is therefore LPCO/2.
c  The last point of AC_IGRAM that can be sampled is NFFTPC- LPCO/2.
c  Points outside this range are invalid: possibly garbage or never computed.
c  Just because points are invalid, doesn't mean that they are definitely
c  wrong.  If LPCO = 8 but only 3 consecutive points of PCO are actually
c  non-zero, then LPCO is effectively only 3 points long and all output
c  points except the first and last will be good (fast convolution only).
c  
c  In the case that PCO is the vector [1,0,0,0,....] then the
c  convolution is a do-nothing operation.

      implicit none
      integer*4 lpco,l2,nfftpc,i,j,method
      real*4  ac_igram(nfftpc),workspace(nfftpc),pco(lpco),xr
      real*8 tot

      l2=lpco/2
c      write(*,*)' CONVOL: lpco,nfftpc,method = ',lpco,nfftpc,method

c================================================================
c  Time/OPD-domain convolution. Slow but does not requires WORKSPACE.
      if(method.eq.1) then  
c         write(*,*)' CONVOLVE: slow convolution'
         do i=1,nfftpc-lpco+1
            tot=0.0d0
            do j=0,l2-1
               tot=tot+
     &         ac_igram(i+j)*pco(l2-j)+      ! Upper half of PCO
     &         ac_igram(i+j+l2)*pco(lpco-j)  ! Lower half of PCO
            end do
            ac_igram(i)=sngl(tot)
         enddo

c  Right-shift corrected igram by L2-1 to line up with raw igram.
c  (This cannot be done as part of main loop)
         do i=nfftpc-lpco+1,1,-1
            ac_igram(i+l2-1)=ac_igram(i)
         end do

c================================================================
c  Perform fast convolution (complex multiplication in spectral domain).
      elseif(method.eq.2) then    ! Perform fast convolution

c  Move PCO to the ends of WORKSPACE. Zero uninitialized middle points.
         do j=1,l2
            workspace(j)=pco(j)
            workspace(nfftpc+1-j)=pco(lpco+1-j)
         end do
         do j=l2+1,nfftpc-l2
            workspace(j)=0.0
         end do
c         do j=1,l2
c            ac_igram(nfftpc+j)=pco(j)
c            ac_igram(2*nfftpc+1-j)=pco(lpco+1-j)
c         end do
c         do j=nfftpc+l2+1,2*nfftpc-l2
c            ac_igram(j)=0.0
c         end do

c         do j=1,nfftpc
c            write(*,*) j,ac_igram(j),ac_igram(j+nfftpc)
c         end do
 
c  Perform a full-size, real-to-complex FFT on PCO to produce an
c  interpolated complex phase spectrum.
         call ffak(workspace,nfftpc)
 
c  Perform real-to-complex FFT on raw interferogram to produce complex spectrum
         call ffak(ac_igram,nfftpc)
 
c  Multiply complex spectra, storing product in AC_IGRAM.
         ac_igram(1)=ac_igram(1)*workspace(1)  ! Real DC term (no imaginary)
         ac_igram(2)=ac_igram(2)*workspace(2)  ! Real Nyquist term (no imaginary)
         do i=4,nfftpc,2
            xr=ac_igram(i-1)*workspace(i-1)-ac_igram(i)*workspace(i) ! Real
            ac_igram(i)=
     &         ac_igram(i-1)*workspace(i)+workspace(i-1)*ac_igram(i) ! Imag
            ac_igram(i-1)=xr
c            cbuf(i/2)=cbuf(i/2)*cbuf(nfftpc/2+i)  ! equivalent to 4 previous lines
         enddo
c
c         ac_igram(1)=ac_igram(1)*ac_igram(nfftpc+1)  ! Real DC term (no imaginary)
c         ac_igram(2)=ac_igram(2)*ac_igram(nfftpc+2)  ! Real Nyquist term (no imaginary)
c         do i=4,nfftpc,2
c            j=nfftpc+i
c            xr=ac_igram(i-1)*ac_igram(j-1)-ac_igram(i)*ac_igram(j) ! Real
c            ac_igram(i)=
c     &         ac_igram(i-1)*ac_igram(j)+ac_igram(j-1)*ac_igram(i) ! Imag
c            ac_igram(i-1)=xr
cc            cbuf(i/2)=cbuf(i/2)*cbuf(nfftpc/2+i)  ! equivalent to 5 previous lines
c         enddo
 
c  Inverse FFT complex product to yield phase corrected igram in lower half of AC_IGRAM.
         call ffsk(ac_igram,nfftpc)
      else
         write(*,*)' convolve.f: error: method = ',method
         stop 'convolve.f: Unknown method value'
      endif    !  if(method.eq.1)
c================================================================

      return
      end
