      subroutine profzl(apo,ns,resnog,rectog,off,a)
      implicit none

      INTEGER*4 k,apo,ns
      REAL*4 a(ns)
      REAL*8 resnog,rectog,off,can,xx,hwid
      real*8 ss,t,sinc,alpha
c
c      write(*,*)'prof:',apo,ns,resnog,rectog
      hwid=0.5d0*(ns-1)
      if(abs(off)+resnog.gt.hwid) then
         write(*,*)'Warning from PROFZL: offset exceeds NK'
         off=dsign(hwid-resnog,off)
      endif
      if(apo.gt.4) stop 'maximum apo is 4'
      
c  Calculate truncated instrumental function (sinx/x for apo=0)
      alpha=0.0
      if(rectog.gt.0.0) alpha=0.12
      can=dpi/resnog
      do k=1,ns
          xx=dble(k)-1.d0-hwid
          t=can*(xx-off)
          ss=sinc(t)
c          a(k)=ss*(1+alpha*(1-ss))
          a(k)=ss*(1+alpha) - alpha/(1+(0.50*t)**2)
          a(k)=a(k)*sngl((1.d0-(xx/(hwid+0.0d0))**2)**2)  ! apodize weakly
      end do
      return
      end

      real*8 function sinc(t)
      real*8 t,tt
      if (dabs(t).lt.0.2d0) then
         tt=t*t
         sinc=1.0d0-tt*(1.0d0-tt/20)/6  !  Approx good to 10^-8 at t=0.2
      else
         sinc=dsin(t)/t
      endif
      return
      end

