      subroutine compute_ils(apo,nhw,ldec,resnog,rectog,off,a)
c  Convolves a SINC function (SIN(X)/X) of half-width RESNOG with a
c  rectangle (box-car) of full-width RECTOG in order to represent the
c  Instrumental Line Shape (ILS) of a perfect FTIR spectrometer.
c
c  The resulting ILS is over-sampled by a factor of LDEC with respect
c  to the spectral grid of the function that it will eventually be
c  convolved with, e.g. the infinite resolution computed spectrum.
c  This oversampling allows accurate interpolation to abscissas
c  that are not on the primative spectral grid.
c
c  This is necessary because, in general, the input vector (infinite
c  -resolution transmittance on a sub-doppler grid) that is being
c  convolved with the ILS is on a different wavenumber grid to the
c  output vector (matching the grid of the measured spectrum).
c
c  So the returned vector (A) actually contains LDEC interleaved ILS's
c  each representing the ILS with a slightly different wavenumber offset.
c  Only one of which is to convolved with a given primitive spectrum
c  (the one closest to the desired output wavenumber).
c  Each of the LDEC ILS subsets has been normalized to unit area.
c  If you don't need this multiple-ILS capability, set LDEC=1.
c
c  Since the SINC function is infinite in extent, the result must be
c  truncated to some finite half-width (NHW points). To minimize the
c  resulting discontinuities, Connes-like windowing is applied to the ILS
c      (1.d0-(x/xmax)^2)^2  
c
c  Various apodization functions are supported including the Norton-Beer
c  ones. 
c
c  In normal use, compute_ils will be called twice: Once for the synthetic
c  spectrum with the actual values of RESNOG and RECTOG, and once for the
c  measured spectrum with RECTOG=0. Convolving the measured spectrum with
c  its own (infinite) SINC function would be a do-nothing operation.
c  However, convolving the measured spectrum with a truncated and
c  windowed version of its own SINC function will improve the agreement
c  with the synthetic spectrum.  In other words, you should apply the same
c  apodization to the measured spectrum as is appled to the calculated
c  spectrum, to maintain an apples-apples comparison.
c
c  INPUTS:
c      APO   I*4  Desired apodization function (0, 1, 2, 3, 4)
c      NHW   I*4  Half-width of operator in # of points (should be > 36*RESNOG)
c   INTERP   I*4  Integer factor by wich ILS is to be oversampled.
c   RESNOG   R*8  = 0.5/OPD/GRID for an FTIR Spectrometer
c   RECTOG   R*8  = FREQ*FOVD**2/8/GRID for an FTIR Spectrometer
c      OFF   R*8  = Frequency Offset / GRID (usually zero)
c
C  OUTPUT:
c     A(1+2*NHW*INTERP)  R*4  Vector containing resulting ILS operator
c
c  OPD  is the Optical Path Difference (cm, usually from runlog)
c  GRID is the desired point spacing (cm-1) of the resulting slit function
c  FREQ is the frequency (cm-1) of interest (usually center of window)
c  FOVD is the diameter of the field-of-view (radians, usually from runlog)
c
      implicit none

      integer*4 i,j,k,apo,np,jp,nhw,ldec
      real*4 a(1+2*nhw*ldec)
      real*8 resnog,rectog,off,c(0:3,0:3),del,can,xx,hwid,tot
      real*8 p,t,t2,t4,t6,t8,q0,q1,q2,q4,tr,dpi
      parameter(dpi=3.14159265359D0)
c
c  Array C contains the coefficients for the Norton-Beer apodization
      save c
      data c/1.0,0.5480,0.2600,0.0900,
     &       0.0,-.0833,-.154838,0.00,
     &       0.0,0.5353,.894838,.5875,
     &       0.0,0.0000,0.0000,0.3225/

c      write(*,*) 'Compute_ILS:',apo,nhw,ldec,resnog,rectog,off

      hwid=dfloat(nhw*ldec)
      if(abs(off)+resnog.gt.hwid) then
         write(*,*)'Warning from compute_ils: offset exceeds NK'
         off=dsign(hwid-resnog,off)
      endif
      if(apo.gt.4) stop 'maximum apo is 4'
c
c  NP = number of (equally weighted) points used to represent the RECT
c  contribution of the ILS. Their point spacing DEL is chosen to match the
c  first 3 moments of the continuous distribution (area, position, & width).
      np=2+int(4*rectog/resnog)
      del=rectog/dsqrt(np*dble(np)-1.d0)
      if(rectog.lt.0.0001d0) np=1   ! approximate RECT by a delta function
c      write(*,*)'np=',np,resnog,rectog
c
c  Calculate truncated instrumental function (sinx/x for apo=0)
      can=dpi/resnog
      do i=1,ldec
         tot=0.0d0
         do j=1,2*nhw  ! Loop over SINC samples
            k=i+ldec*(j-1)
            a(k)=0.0
            xx=dble(k)-1.d0-hwid
            do jp=-np+1,np-1,2  ! Loop over RECT samples
               T=can*(xx-off+jp*del/2)
               T2=T*T
               t4=t2*t2
               if (T2.GE.1.2D0) then
                  Q0=DSIN(T)/T
                  P=DCOS(T)
                  Q1=3*(Q0-P)/T2
                  tr=2*(1.d0-p)/t2     ! = sinc(T/2)**2
                  Q2=-(15*((1-3/T2)*Q0+3*P/T2)/T2)
                  Q4=945*((1-45/T2+105/t4)*Q0+5*(2-21/T2)*P/T2)/T4
               else
                  t6=t2*t4
                  t8=t2*t6
                  q0=1-t2/6 +t4/120 -t6/5040   +t8/362880
                  q1=1-t2/10+t4/280 -t6/15120  +t8/1330560
                  tr=1-t2/12+t4/360 -t6/20160  +t8/1814400
                  q2=1-t2/14+t4/504 -t6/33264  +t8/3459456
                  q4=1-t2/22+t4/1144-t6/102960 +t8/14002560
               endif
               if(apo.eq.4) then
                  a(k)=a(k)+sngl(tr)  ! 'TR' apodization  = sinc(T/2)**2
               elseif(apo.le.0) then
                  a(k)=a(k)+sngl(Q0)
               else
                  a(k)=a(k)+
     &            sngl(C(apo,0)*Q0+C(apo,1)*Q1+C(apo,2)*Q2+C(apo,3)*Q4)
               endif
            enddo   !  jp=-np+1,np-1,2  ! Loop over RECT samples
            a(k)=a(k)*sngl((1.d0-(xx/(hwid+0.0d0))**2)**2)  ! apodize weakly
            tot=tot+a(k)
         end do  !  j=1,2*nhw   ! Loop over SINC samples

c  Normalize each ILS subset to unity.
         k=i
         do j=1,2*nhw
            a(k)=a(k)/sngl(tot)
            k=k+ldec
         end do

      end do  !  i=1,ldec  ! Loop over SINC samples
      a(1+2*nhw*ldec)=0.0

      return
      end
