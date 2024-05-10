      function nlo_bs(nlev,z,zz)

c  Finds the model levels which bracket the altitude zz.
c      z(nlo) < zz < z(nhi)
c  or
c      z(nlo) > zz > z(nhi)
c
c  z(i) must be ordered, monotonically increasing or decreasing.
c  Binary search is used.
c  After convergence, nhi_bs=nlo_bs+1, so no need to return both

      integer*4 nlev
      real*4 z(nlev),zz

      if(nlev.le.1) stop 'nlo_bs:  nlev <= 1'
      nlo_bs=1
      nhi_bs=nlev
      dz=z(nhi_bs)-z(nlo_bs)         ! Does z(i) increase with i or decrease ?
      do while(nhi_bs-nlo_bs.gt.1)
         new=nlo_bs+(nhi_bs-nlo_bs)/2   ! Bisect remaining range, avoiding overflow
         if((z(new)-zz)/dz.lt.0.0) then  ! Move in the redundant end point
            nlo_bs=new
         else
            nhi_bs=new
         endif
      end do

      return
      end
