      subroutine read_mavfile_body(lunr_mav,nlev_ray,nspeci_iso,
     & z,t,p,d,vmr,parfile)
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
c       VMR(nspeci_iso,nlev_ray)   Array of vmrs
c
      implicit none
      integer*4 lunr_mav,nlev_ray,nlev_mav,nspeci_iso,jspeci,
     & nlhead,j,ilev,lunr_iso,
     & ncol,kcol,nspeci_mav,kgas,kiso

      include '../ggg_int_params.f'
      parameter (lunr_iso=19)
      integer*4 sindex(mspeci)
      real*4 z(nlev_ray),t(nlev_ray),p(nlev_ray),d(nlev_ray),
     & vmr(nspeci_iso,nlev_ray),
     & zero

      character header_string_mav*2000,header_vec_mav(mspeci+4)*10,
     & sss*8,
     & shortname(mspeci)*9,parfile*(*)
c
      zero=0.0
      read(lunr_mav,*)nlhead,ncol,nlev_mav
      if (nlev_mav.gt.nlev_ray) then
           write(*,*)'read_mavfile_body: nlev_mav .gt. nlev_ray',
     &     nlev_mav,nlev_ray
           stop ' READ_MAV: mismatch in nlev'
      endif

c      call skiprec(lunr_mav,nlhead-1)  !   Skip header
      call skiprec(lunr_mav,nlhead-2)  !   Skip header
      read(lunr_mav,'(a)') header_string_mav
      call substr(header_string_mav,header_vec_mav,mspeci,kcol)
c      write(*,*) kcol,ncol
      if(kcol.ne.ncol) stop 'Format Error in .mav file'
      nspeci_mav=kcol-4

c  Read isotopologs.dat to get gas names
c  Initialize the VMR array with zeros.
      open(lunr_iso,file=parfile,status='old')
      do jspeci=1,mspeci
         sindex(jspeci)=0
         read(lunr_iso,'(i3,i2,x,a8)',end=97) kgas,kiso,sss
         call lowercase(sss)
         shortname(jspeci)=char(kiso+48)//sss
         call vmov(zero,0,vmr(jspeci,1),nlev_ray,nlev_ray)
      end do
97    close(lunr_iso)
      if(jspeci-1.ne.nspeci_iso) stop 'jspeci-1 .ne. nspeci_iso'
 
c  Index specie names in .mav file header to those in isotopologs.dat
c  First 4 elements of header_vec_mav are Z,T,P,D
      call clistindex(nspeci_mav,header_vec_mav(5),nspeci_iso,shortname,
     & sindex)

c  Run checks on .mav file header content (specie names)
      do jspeci=1,nspeci_mav

c  Warn against different column order between .mav & isotopologs files
         if(jspeci.ne.sindex(jspeci)) then
         write(*,*)'Column order in isotopologs.dat does not match .mav'
            write(*,*)jspeci,shortname(jspeci),header_vec_mav(jspeci+4),
     &      sindex(jspeci)
            write(*,*)'Re-run GSETUP to get rid of this warning'
         endif

c  Check that species in .mav file header are present in isotopologs.dat
         if(sindex(jspeci).le.0) then
            write(*,*) j,header_vec_mav(jspeci+4)
            stop 'Unidentified specie in .mav file'
         endif

c  Check for duplicated specie entrys in .mav file
         do j=1,jspeci-1
            if(header_vec_mav(jspeci).eq.header_vec_mav(j)) then
            write(*,*)jspeci,header_vec_mav(jspeci),j,header_vec_mav(j)
            stop 'Duplicated specie entry in .mav file'
            endif
         end do

      end do  ! jspeci=1,nspeci_mav

      do ilev=1,nlev_mav
c         read(lunr_mav,'(2f7.2,2(1pe11.3),150(e11.3))')
         read(lunr_mav,*) z(ilev),t(ilev),p(ilev),d(ilev),
     &   (vmr(sindex(j),ilev),j=1,nspeci_mav)
         if(t(ilev).eq.0.0) stop 'zero temperature'
      end do
      return
      end
