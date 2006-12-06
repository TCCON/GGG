      function height(pobs,p,t,z,nl)
c  Finds the altitude (HEIGHT) at which the pressure POBS occurs.
c  Interpolates linearly in log(p) to determine z corresponding to pobs
c  FTA is a correction for the temperature non-uniformity.
      implicit none
      real*4 fta,ftawas,tobs,fr,pobs,height,p(*),t(*),z(*)
      integer*4 k,nl,it
      do 11 k=2,nl
        if(p(k).le.pobs) goto 12
 11   continue
      write(6,*)'nl,Z(top),P(top), Pobs =',k-1,z(k-1),p(k-1),pobs
      stop 'Error in function HEIGHT: insufficiently small p'
 12   fr=log(p(k-1)/pobs)/log(p(k-1)/p(k))
c      write(6,*) p(k),z(k),pobs,p(k-1),z(k-1)
      fta=1.
      do it=1,10
        ftawas=fta
        tobs=t(k-1)+fta*fr*(t(k)-t(k-1))
        fta=(t(k-1)+tobs)/(t(k-1)+t(k))
        if(abs(fta-ftawas).lt.2.e-7) go to 13
      end do
      write(*,*) pobs,tobs,fr,fta,ftawas,z(k-1)+fta*fr*(z(k)-z(k-1))
      stop 'HEIGHT failed to converge'
 13   height=z(k-1)+fta*fr*(z(k)-z(k-1))
c      write(6,*) pobs,tobs,height,p(k-1),p(k)
      return
      end
