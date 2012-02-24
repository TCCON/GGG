      subroutine read_mavfile_body(lun_mav,nlev,nspeci,z,t,p,d,vmr)
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
c       VMR(nspeci,nlev)   Array of vmrs
c
      implicit none
      integer*4 lun_mav,nlev,nlevmav,nspeci,nlhead,j,ilev,ncol
      real*4 z(nlev),t(nlev),p(nlev),d(nlev),vmr(nspeci,nlev)
c
      read(lun_mav,*)nlhead,ncol,nlevmav
      if (nlevmav.ne.nlev) then
           write(*,*)'read_mavfile_body: nlevmav .ne. nlev',
     &     nlevmav,nlev
           stop ' READ_MAV: mismatch in nlev'
      endif
      if(ncol-4 .ne. nspeci) then
           write(*,*)'ERROR: ncol-4 nspeci =',ncol-4,nspeci
           write(*,*)'The number of species in isotopolog.dat differs'
           write(*,*)'from the number in your .mav file. Re-run GSETUP.'
           stop ' READMAV: mismatch in NSPECI'
      endif
      call skiprec(lun_mav,nlhead-1)  !   Skip header

c      write(*,*)'READ_MAVFILE_BODY: nlev,nspeci=',nlev,nspeci
      do ilev=1,nlev
c         read(lun_mav,'(2f7.2,2(1pe11.3),150(e11.3))')
         read(lun_mav,*)
     &   z(ilev),t(ilev),p(ilev),d(ilev),(vmr(j,ilev),j=1,nspeci)
         if(t(ilev).eq.0.0) stop 'zero temperature'
      end do
      return
      end
