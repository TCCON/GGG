      subroutine profzl(apo,ns,resnog,rectog,off,a)
c  Convolves a SINC function (SIN(X)/X) of half-width RESNOG with a rectangle
c  (box-car) of full-width RECTOG in order to represent the Instrumental Line
c  Shape (ILS) of a perfect FTIR spectrometer.
cc
c  Since the SINC function is infinite in extent, the result must be truncated
c  to some finite number (NS) of points, and must therefore also be apodized
c  to avoid discontinuities at the ends of the operator.  Various apodization
c  functions can be selected including the Norton-Beer ones. However, note that
c  even if APO=0 is chosen, slight apodization will still be applied.
c
c  In normal use, PROFZL will be called twice: Once for the synthetic spectrum
c  with the actual values of RESNOG and RECTOG, and once for the measured
c  spectrum with RECTOG=0. Convolving the measured spectrum with its own
c  (infinite) SINC function would be a do-nothing operation. However, convolving
c  the measured spectrum with a finite and weakly apodized version of its own
c  SINC function will improve the agreement with the synthetic spectrum.
c
c  INPUTS:
c      APO  I*4  Desired apodization function (0, 1, 2, 3, 4)
c       NS  I*4  Number of points in the operator (should be > 36*RESNOG)
c   RESNOG  R*8  = 0.5/OPD/GRID for an FTIR Spectrometer
c   RECTOG  R*8  = FREQ*FOVD**2/8/GRID for an FTIR Spectrometer
c      OFF  R*8  = Frequency Offset / GRID (usually zero)
c
C  OUTPUT:
c     A(NS) R*4  Array containing resulting operator
c
c  OPD is the Optical Path Difference (cm)
c  GRID is the desired point spacing (cm-1) of the resulting slit function
c  FREQ is the frequency (cm-1) of interest 
c  FOVD is the diameter of the field-of-view (radians)
c
      implicit none
      INTEGER*4 k,apo,ns,np,jp
      REAL*4 a(ns),zero
      REAL*8 resnog,rectog,off,c(0:3,0:3),del,can,pi,xx,hwid
      real*8 p,t,t2,t4,t6,t8,q0,q1,q2,q4,tr
      SAVE C
      parameter (pi=3.14159265d0,zero=0.0)
      data c/1.0,0.5480,0.2600,0.0900,
     &       0.0,-.0833,-.154838,0.00,
     &       0.0,0.5353,.894838,.5875,
     &       0.0,0.0000,0.0000,0.3225/
C
      hwid=0.5d0*(ns-1)
      if(abs(off)+resnog.gt.hwid) then
         write(*,*)'Warning from PROFZL: offset exceeds NK'
         off=dsign(hwid-resnog,off)
      endif
      if(apo.gt.4) stop 'maximum apo is 4'
c
c      if(apo.lt.0) then
c         call vmov(zero,0,a,1,ns)
c         a((ns+1)/2)=1.0
c         return
c      endif
      
c  NP= number of (equally weighted) points used to represent the rectangular
c  contribution of the ILS. Their point spacing DEL is chosen to match the first
c  three moments of the continuous distribution (area, position, and hwidth).
      np=2+int(4*rectog/resnog)
      del=rectog/dsqrt(np*dble(np)-1.d0)
      if(rectog.lt.0.0001d0) np=1   ! approximate RECT by a delta function
c      write(*,*)'np=',np,resnog,rectog
c
c  Calculate truncated instrumental function (sinx/x for apo=0)
      can=pi/resnog
      DO 1040 k=1,ns
      a(k)=0.0
      xx=dble(k)-1.d0-hwid
      do 1041 jp=-np+1,np-1,2
      T=can*(xx-off+jp*del/2)
      T2=T*T
      t4=t2*t2
      IF (T2.GE.1.2D0) THEN
        Q0=DSIN(T)/T
        P=DCOS(T)
        Q1=3*(Q0-P)/T2
        tr=2*(1.d0-p)/t2     ! = sinc(T/2)**2
        Q2=-(15*((1-3/T2)*Q0+3*P/T2)/T2)
        Q4=945*((1-45/T2+105/t4)*Q0+5*(2-21/T2)*P/T2)/T4
      ELSE
        t6=t2*t4
        t8=t2*t6
        q0=1-t2/6 +t4/120 -t6/5040   +t8/362880
        q1=1-t2/10+t4/280 -t6/15120  +t8/1330560
        tr=1-t2/12+t4/360 -t6/20160  +t8/1814400
        q2=1-t2/14+t4/504 -t6/33264  +t8/3459456
        q4=1-t2/22+t4/1144-t6/102960 +t8/14002560
      ENDIF
      if(apo.eq.4) then
         a(k)=a(k)+sngl(tr)  ! 'TR' apodization  = sinc(T/2)**2
      else
         a(k)=a(k)+sngl(C(apo,0)*Q0+C(apo,1)*Q1+C(apo,2)*Q2+C(apo,3)*Q4)
      endif
1041  CONTINUE
      a(k)=a(k)*sngl((1.d0-(xx/(hwid+0.0d0))**2)**2)  ! apodize weakly
c      a(k)=a(k)*(cos(3.14159265d0*xx/hwid/2))**4  ! apodize weakly
1040  CONTINUE
      return
      end
