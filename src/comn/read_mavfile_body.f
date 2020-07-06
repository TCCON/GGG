      subroutine read_mavfile_body(lunr_mav,nlev_ray,nspeci_iso,
     & nlev_mav,z,t,p,d,wmf,parfile)
c  Reads a block of data from the already open .mav file
c  into arrays Z, T, P, D, WMF.
c
c  Inputs:
c       LUNR_MAV     I*4    Logical unit number for reading .mav file
c       NLEV_RAY     I*4    Dimension of Z, T, P, D vectors
c       NSPECI_ISO   I*4    Number of species in isotopologs.dat
c       PARFILE      C*(*)  Isotopologs.dat path (to get gas names)  
c
c Outputs:
c       NLEV_MAV     I*4    Number of levels in .mav file
c       Z(NLEV_RAY)  R*4    Array of altitudes
c       T(NLEV_RAY)  R*4    Array of temperatures
c       P(NLEV_RAY)  R*4    Array of pressures
c       D(NLEV_RAY)  R*4    Array of densities
c       WMF(nspeci_iso,nlev_ray)   Wet Mole Fractions
c
      implicit none
      integer*4 lunr_mav,nlev_ray,nlev_mav,nspeci_iso,jspeci,
     & nlhead,ncol_iso,j,ilev,lunr_iso,
     & ncol_mav,kcol,nspeci_mav,kgas,kiso,idum

      include '../gfit/ggg_int_params.f'
      parameter (lunr_iso=19)
      integer*4 sindex(mspeci)
      real*4 z(nlev_ray),t(nlev_ray),p(nlev_ray),d(nlev_ray),
c     & dmf(nspeci_iso),
     & wmf(nspeci_iso,nlev_ray),
     & zero

      character header_string_mav*2248,header_vec_mav(mspeci+4)*10,
     & sss*8,col1*1,
     & shortname(mspeci)*9,parfile*(*)

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mgas      ! Avoid compiler warning (unused parameter)
      idum=mlev      ! Avoid compiler warning (unused parameter)
      idum=mvmode    ! Avoid compiler warning (unused parameter)
      idum=mauxcol   ! Avoid compiler warning (unused parameter)
      idum=mcolvav   ! Avoid compiler warning (unused parameter)
      idum=mrow_qc   ! Avoid compiler warning (unused parameter)
      idum=mspeci    ! Avoid compiler warning (unused parameter)
      idum=ncell     ! Avoid compiler warning (unused parameter)
      idum=nchar     ! Avoid compiler warning (unused parameter)

      zero=0.0
      read(lunr_mav,*)nlhead,ncol_mav,nlev_mav
c      write(*,*) 'nlhead,ncol_mav,nlev_mav=',nlhead,ncol_mav,nlev_mav
c      if (nlev_mav.gt.nlev_ray) then
c           write(*,*)'read_mavfile_body: nlev_mav .gt. nlev_ray',
c     &     nlev_mav,nlev_ray
c           stop ' READ_MAV: mismatch in nlev'
c      endif

c      call skiprec(lunr_mav,nlhead-1)  !   Skip header
      call skiprec(lunr_mav,nlhead-2)  !   Skip header
      read(lunr_mav,'(a)') header_string_mav
c      write(*,*) 'header_string_mav = '//header_string_mav
      call substr(header_string_mav,header_vec_mav,mspeci,kcol)
      if(kcol.ne.ncol_mav) then
         write(*,*)'kcol,ncol_mav=', kcol,ncol_mav
         stop 'Format Error in .mav file: kcol .ne. ncol_mav'
      endif
      nspeci_mav=kcol-4

c  Read isotopologs.dat to get gas names
c  Initialize the WMF array with zeros.
      open(lunr_iso,file=parfile,status='old')
      read(lunr_iso,*) nlhead,ncol_iso
      do j=2,nlhead
         read(lunr_iso,*)
      end do
      do jspeci=1,mspeci
         sindex(jspeci)=0
1        read(lunr_iso,'(a1,i2,i2,x,a8)',end=97) col1,kgas,kiso,sss
         if(col1.eq.':' .or. col1.eq.';') go to 1
         call lowercase(sss)
         if(kiso.le.9) then
            write(shortname(jspeci),'(i1,a8)')kiso,sss
         else
            write(shortname(jspeci),'(i2,a7)')kiso,sss
         endif
c         call vmov(zero,0,wmf(jspeci,1),nlev_ray,nlev_ray)
         do ilev=1,nlev_ray
            wmf(jspeci,ilev)=0.0
         end do

      end do  !  do jspeci=1,mspeci
97    close(lunr_iso)
      if(jspeci-1.ne.nspeci_iso) then
         write(*,*) jspeci-1, nspeci_iso
         stop 'jspeci-1 .ne. nspeci_iso'
      endif
 
c  Index specie names in .mav file header to those in isotopologs.dat
c  First 4 elements of header_vec_mav are Z,T,P,D
      call clistindex(nspeci_mav,header_vec_mav(5),nspeci_iso,shortname,
     & sindex)

c  Run checks on .mav file header content (specie names)
      do jspeci=1,nspeci_mav
c         write(*,*)'jspeci,sindex(jspeci)',jspeci,sindex(jspeci)
c  Warn against different column order between .mav & isotopologs files
         if(jspeci.ne.sindex(jspeci)) then
            write(*,*)'Warning: read_mavfile_body.f: Column order in'//
     &      ' isotopologs.dat does not match that in .mav'
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
               write(*,*)jspeci,header_vec_mav(jspeci),header_vec_mav(j)
               stop 'Duplicated specie entry in .mav file'
            endif
         end do

      end do  ! jspeci=1,nspeci_mav

      do ilev=1,nlev_mav
         read(lunr_mav,*) z(ilev),t(ilev),p(ilev),d(ilev),
     &   (wmf(sindex(j),ilev),j=1,nspeci_mav)
c     &   (dmf(sindex(j)),j=1,nspeci_mav)
         if(t(ilev).le.0.0) stop 'Bad temperature'
c         do j=1,nspeci_mav   ! Convert from Dry to Wet Mole Fractions
c            wmf(sindex(j),ilev)=dmf(sindex(j))/(1+dmf(1)) ! dmf(1) is H2O
c         end do
      end do
      return
      end
