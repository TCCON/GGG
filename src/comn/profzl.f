      subroutine profzl(apo,ns,resnog,rectog,off,a)
c  Convolves a SINC function (SIN(X)/X) of half-width RESNOG with a
c  rectangle (box-car) of full-width RECTOG in order to represent
c  the Instrumental Line Shape (ILS) of a perfect FTIR spectrometer.
c
c  In a measured spectrum the ILS is infinite in extent. So the ringing
c  from a sharp feature can sometimes be seen to extend for thousands
c  of points, if the spectrum is otherwise sufficiently featureless,
c  and the self-apodization (due to the finite FOV) weak. In a
c  low-res spectrum this ringing can represent hundreds of cm-1.
c  To avoid having to use such a wide computed ILS, we truncate it at
c  a distance where the ringing is only a few % of the central peak.
c  To minimize the resulting discontinuities, a Connes-like apodization
c  is then applied
c        ILS(v-vo) = ILS(v-vo) * (1-((v-vo)/a)^2)^2  
c  where a is the wavenumber from center at which truncation occurs.
c  The ILS is then normalized to unit area, following this apodization.
c
c  For an apples/apples comparison between measured and computed spectra
c  the same truncation/apodizing function must be applied to both.
c  So since the computed spectrum will have the Connes-like apodization,
c  so must the measured spectrum.  This means that the "measured"
c  spectrum that appears in the GFIT spectral fits is never the original
c  measured spectrum. It is always slightly smoothed in order to kill
c  the distant ringing. This is always the case, independent of what
c  additional apodizations (e.g. Norton-Beer) may have been selected
c  in the .ggg file.
c
c  In normal use, PROFZL will be called twice: Once for the synthetic
c  spectrum with the actual values of RESNOG and RECTOG, and once for
c  the measured spectrum with RECTOG=0. 
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
c  OPD  is the Optical Path Difference (cm)
c  GRID is the desired point spacing (cm-1) of the resulting slit function
c  FREQ is the frequency (cm-1) of interest 
c  FOVD is the diameter of the field-of-view (radians)
c
      implicit none

      INTEGER*4 k,apo,ns,np,jp
      REAL*4 a(ns)
      REAL*8 resnog,rectog,off,c(0:3,0:3),del,can,xx,hwid
      real*8 p,t,t2,t4,t6,t8,q0,q1,q2,q4,tr
c
c  Array C contains the coefficients for the Norton-Beer apodization
      SAVE C
      data c/1.0,0.5480,0.2600,0.0900,
     &       0.0,-.0833,-.154838,0.00,
     &       0.0,0.5353,.894838,.5875,
     &       0.0,0.0000,0.0000,0.3225/
c
c      write(*,*)'profzl:',apo,ns,resnog,rectog
      hwid=0.5d0*(ns-1)
      if(abs(off)+resnog.gt.hwid) then
         write(*,*)'Warning from PROFZL: offset exceeds NK'
         off=dsign(hwid-resnog,off)
      endif
      if(apo.gt.4) stop 'maximum apo is 4'
c
cc  ILS is a delta-function
c      if(apo.lt.0) then
c         call vmov(zero,0,a,1,ns)
c         a((ns+1)/2)=1.0
c         return
c      endif
      
c  NP = number of (equally weighted) points used to represent the RECT
c  contribution of the ILS. Their point spacing DEL is chosen to match the first
c  three moments of the continuous distribution (area, position, and hwidth).
      np=2+int(4*rectog/resnog)
      del=rectog/dsqrt(np*dble(np)-1.d0)
      if(rectog.lt.0.0001d0) np=1   ! approximate RECT by a delta function
c      write(*,*)'np=',np,resnog,rectog
c
c  Calculate truncated instrumental function (sinx/x for apo=0)
      can=3.14159265359D0/resnog
      DO 1040 k=1,ns  ! Loop over SINC samples
      a(k)=0.0
      xx=dble(k)-1.d0-hwid
      do 1041 jp=-np+1,np-1,2  ! Loop over RECT samples
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
      elseif(apo.le.0) then
         a(k)=a(k)+Q0
      else
         a(k)=a(k)+sngl(C(apo,0)*Q0+C(apo,1)*Q1+C(apo,2)*Q2+C(apo,3)*Q4)
      endif
1041  CONTINUE
      a(k)=a(k)*sngl((1.d0-(xx/(hwid+0.0d0))**2)**2)  ! apodize weakly
1040  CONTINUE
      return
      end
