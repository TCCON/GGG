      function height(pobs,p,t,z,nl)
c  Finds the altitude (HEIGHT) at which the pressure POBS occurs
c  in the atmospheric model defined by z/p/t.
c  Interpolates linearly in log(p) to determine z corresponding
c  to pobs.  FTA is a correction for the temperature non-uniformity.
      implicit none
      integer*4 k,nl,it
      real*4 fta,ftawas,tobs,fr,pobs,height,p(nl),t(nl),z(nl)
      do k=2,nl
c        write(*,*)k,p(k),pobs
        if(p(k).le.pobs) goto 12
      end do
      write(6,*)'nl,Z(top),P(top), Pobs =',k-1,z(k-1),p(k-1),pobs
      write(6,*) 'Warning from function HEIGHT: insufficiently small p'
      height=z(nl)
      return
 12   fr=log(p(k-1)/pobs)/log(p(k-1)/p(k))
c      write(*,*) k,p(k),z(k),pobs,p(k-1),z(k-1)
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
