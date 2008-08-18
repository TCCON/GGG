      subroutine read_mav(lunm,mlev,nlev,nspeci,z,t,p,d,vmr)
c  Reads the contents of the .mav file MAVNAME into the appropriate arrays.
c
c  Inputs:
c       LUNM           Logical unit number
c       MLEV           Declared dimension of arrays Z, T, P, D
c       NLEV           Number of levels of MAVFILE to be read
c       NLOW           The lowest model level that needs to be read
c       NSPECI         Actual number of different gases in MAVFILE
c
c Outputs:
c       Z              Array of altitudes
c       T              Array of temperatures
c       P              Array of pressures
c       D              Array of densities
c       VMR(nspeci,mlev)   Array of vmrs
c
      implicit none
      integer*4 lunm,mlev,nlev,nlevmav,nspeci,nlhead,j,ilev,ncol
      real*4 z(mlev),t(mlev),p(mlev),d(mlev),vmr(nspeci,mlev)
c
      read(lunm,*)nlhead,ncol,nlevmav
      if (nlevmav.ne.nlev) then
           write(*,*)'read_mav: nlevmav .ne. nlev',
     &     nlevmav,nlev
           stop ' READ_MAV: mismatch in nlev'
      endif
      if(ncol-4 .ne. nspeci) then
           write(*,*)'ncol-4 nspeci =',ncol-4,nspeci
           stop ' READMAV: mismatch in NSPECI'
      endif
      call skiprec(lunm,nlhead-1)  !   Skip header

c      write(*,*)'READMAV: nlev,nspeci=',nlev,nspeci
      do ilev=1,nlev
         read(lunm,'(2f7.2,2(1pe11.3),150(e11.3))')
     &   z(ilev),t(ilev),p(ilev),d(ilev),(vmr(j,ilev),j=1,nspeci)
         if(t(ilev).eq.0.0) stop 'zero temperature'
      end do
      return
      end
