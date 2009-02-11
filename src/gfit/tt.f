      integer*4 ifirst,ilast
      real*8 nus,nue,graw,fcen,wid,sw

      graw=-0.01d0
      ifirst=-290000
      ilast=-300000
      sw=1.0d0
      write(*,*)' NUS, NUE'
      read(*,*) nus,nue

      vbar=graw*(ilast+ifirst)/2
      vwid=graw*(ilast-ifirst)/2
      write(*,*)nus,nue,vbar-vwid,vbar+vwid
      if(nus-sw/2 .lt. vbar-vwid) 
     & write(*,*)'Error1:',nus-sw/2,vbar-vwid
      if(nue+sw/2 .gt. vbar+vwid) 
     & write(*,*)'Error2:',nue+sw/2,vbar+vwid
      stop
      end

