      function plint(ax,ay,xx,np)
c  Performs piece-wise linear interpolation of function AY.
c
c  Inputs:
c    AX(N)   R*8  Array of X-values
c    AY(N)   R*8  Array of Y-values
c    XX      R*8  x-value at which to evaluate function
c    NP      I*4  Number of points on AX and AY vectors.
c
c  Outputs:
c    PLINT   R*8  value of Y interpolated to XX (i.e. Y(XX)).
c
c  AX must be monotonically increasing or decreasing function.
c  If AX is equally spaced, you should not use PLINT. There are
c  much more efficient methods available in this case.
c
       implicit none
       integer*4 ip,np
       real*8 ax(np), ay(np), xx, fr, dx, dwas, plint
c
       if(np.lt.2) stop 'PLINT: NP < 2'
       dx=ax(2)-ax(1)
       do ip=2,np
         dwas=dx
         dx=ax(ip)-ax(ip-1)
         if(dx/dwas.le.0) stop 'PLINT: Non-monotonic x-array. '
         fr=(ax(ip)-xx)/dx
         if (fr.ge.0.0) go to 99
       end do
       ip=np
99     plint=fr*ay(ip-1)+(1.-fr)*ay(ip)
       if(fr.lt.0 .or. fr.gt.1) then
          write(*,*)'PLINT: Out of range. '
          write(*,*)xx,ax(1),ax(np)
       endif
       return
       end
