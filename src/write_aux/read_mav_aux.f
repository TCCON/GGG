      subroutine read_mav_aux(lunm,nlhead,nlev,nspeci,z,t,p,d,vmr,
     & vmrfile,modfile,head)
c  Reads the contents of the .mav file MAVNAME into the appropriate arrays.
c
c  Inputs:
c       LUNM           Logical unit number
c       NLHEAD         Number of header lines
c       NLEV           Number of levels of MAVFILE to be read
c       NSPECI         Actual number of different gases in MAVFILE
c
c Outputs:
c       Z              Array of altitudes
c       T              Array of temperatures
c       P              Array of pressures
c       D              Array of densities
c       VMR(nspeci,mlev)   Array of vmrs
c       HEAD           Header line describing contents of file
c
      implicit none
      integer*4 lunm,nlev,nspeci,nlhead,j,ilev
      real*4 z(nlev),t(nlev),p(nlev),d(nlev),vmr(nspeci,nlev)
      character head*1800,vmrfile*1800,modfile*1800

c     call skiprec(lunm,nlhead-1)  !   Skip header
      call skiprec(lunm,nlhead-4)  !   Skip header DW 20100322
      read(lunm,'(a)')vmrfile
      read(lunm,'(a)')modfile
      read(lunm,'(a)')head ! DW 20100322
c     write(*,*)head   ! DW 20100322

c      write(*,*)'READ_MAV: nlev,nspeci=',nlev,nspeci
      do ilev=1,nlev
         read(lunm,'(2f7.2,2(1pe11.3),200(e11.3))')
     &   z(ilev),t(ilev),p(ilev),d(ilev),(vmr(j,ilev),j=1,nspeci)
         if(t(ilev).eq.0.0) stop 'zero temperature'
      end do
      return
      end
