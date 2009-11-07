C  Program GFIT
c  See ggg.history for description of latest changes

      implicit none
c
      integer*4
     & lun_col,lun_iso,lun_cs,nss,i,getpid,
     & istat, lnbc,lp,lc,
     & mmode,nmode,         ! Number of vibrational modes
     & nlhead_ggg,          ! number of header lines in .ggg file
     & kgas,kiso,
     & mtg, ntg, jtg,       ! number of target molecules
     & mspeci,nspeci,jspeci ! number of different species listed in PARFILE
c
      parameter (lun_col=14,lun_cs=16,lun_iso=18,
     & mspeci=141,mtg=16,mmode=30)

      integer*4
     & targmol(mspeci), ! Group assignment for each specie.
     & speci(mtg),      ! speci # of the parent isotopolog of each target gas
     & dgen(mmode),     ! degeneracy of vibrational modes
     & molewt,          ! Molar Mass
     & system           ! For calling system operations
c
      character
     & colfile*80,       ! output file containing column amounts
     & akfile*80,        ! output file of averaging kernels
     & sptfile*80,       ! output file of ascii spectral fits
     & mavfile*80,       ! file containing T/P & VMR at user-chosen levels
     & rayfile*80,       ! file of slant paths at user-chosen levels
     & parfile*80,       ! path to "molparam.dat"
     & apvalerr*80,      ! path to a priori variable values and uncertainties
c     & llsize*80,        ! path to llsize.dat file
     & linefiles*400,    ! paths to linelists
     & solarll*80,       ! solar linelist
     & gfit_version*64,  ! gfit version number
     & gsversion*64,     ! GSETUP version number
     & dplist*80,        ! Data Partition List (e.g. m4part.lst)
     & runlog*80,        ! name of occultation file
     & gasname*8,        ! names of species in PARFILE
c     & tname*8,  ! names of species in PARFILE prefixed by "t"
c     & fullname*10,! full names of species in PARFILE
     & winfo*160,         !  window information (command line)
     & pars(mtg)*9,     ! parameters to fit
     & linelists(9)*80,     ! linelists
     & checksum*32,      ! md5sum checksum value
     & csfilename*32,    ! file containing checksums
     & levfile*80,       ! name of file containing fitting levels
     & amodel*80,        ! name of atmospheric model
     & vmrset*80,        ! name of vmr set
     & window*80         ! name of window list
c
      real*4
     &   tdrpf,        ! T-Dependence of Rotational Partition Function
     &   fia,delta,epsilon,  ! Fractional Isotopic Abundance
     &   vibfrq(mmode) ! Array of vibrational frequencies
c
      character
     &   speci_id*24,
     &   root*128, dl*1

      data speci/mtg*0/
      
      gfit_version=
     & ' GFIT                     Version 4.4.10   07-Nov-2009    GCT '
      write(6,*)
      write(6,'(a)')gfit_version

c     Platform specification:      DG090519
      call get_ggg_environment(root, dl)
c
c  Read runlog, model, vmrset & window information from input file (.ggg)

      read(5,*)nlhead_ggg
      read(5,'(a)')gsversion
      read(5,'(a)')dplist
      read(5,'(a)')apvalerr
      read(5,'(a)')runlog
      read(5,'(a)')levfile
      read(5,'(a)')amodel
      read(5,'(a)')vmrset
      read(5,'(a)')mavfile
      read(5,'(a)')rayfile
      read(5,'(a)')parfile
      read(5,'(a)')window
c      read(5,'(a)')llsize
      read(5,'(a)')linefiles
      read(5,'(a)')solarll
      read(5,'(a)')akfile
      read(5,'(a)')sptfile
      read(5,'(a)')colfile
      read(5,'(a)')winfo

      lc=index(winfo,'#')
      if(lc.gt.0) winfo=winfo(:lc-1)
      lp=lnbc(winfo)
      lc=index(winfo,':')
      write(6,'(a)')winfo(:lp)
      call lowercase(winfo)
      call substr(winfo(lc+1:),pars,mtg,ntg)
      if(ntg.gt.mtg) then
          write(*,*)' gfit: Error: NTG > MTG ',ntg,mtg
          stop 'Increase parameter MTG in gfit.f'
      endif
      if( index(winfo,' cf ').gt.0)write(*,*)'Fitting channel fringes'
c
c  Read in names of isotopomers
      open(lun_iso,file=parfile,status='old')
      do jspeci=1,mspeci
         call read_isotop(lun_iso,kgas,kiso,gasname,speci_id,
     &   fia,delta,epsilon,molewt,tdrpf,vibfrq,dgen,nmode,mmode,istat)
         if(istat.ne.0) go to 77
         call lowercase(gasname)
         targmol(jspeci)=0
         do jtg=1,ntg
            if(char(kiso+48)//gasname.eq.pars(jtg)) targmol(jspeci)=jtg
            if( gasname//' '.eq.pars(jtg) .or. 
     &          'a'//gasname.eq.pars(jtg) .or.
     &          'b'//gasname.eq.pars(jtg) .or.
     &          't'//gasname.eq.pars(jtg)) then
                if(targmol(jspeci).eq.0) targmol(jspeci)=jtg
            endif
         end do
      end do
      read(lun_iso,*,end=77)
      stop ' GFIT: Number of species listed in PARFILE exceeds MSPECI'
77    close(lun_iso)
      nspeci=jspeci-1
c
      do jspeci=nspeci,1,-1
           if(targmol(jspeci).gt.0) speci(targmol(jspeci))=jspeci
      end do
c
c========================================================================
c The first time that spectrum_loop is called it skips the
c CPU intensive operations (e.g. computing VACs, spectral fitting)
c and simply checks that the spectra, .mav, and .ray files are okay
c before investing alot of time on computing absorption coeffs.
      write(*,*) ' Pre-screening input files...'
      call spectrum_loop(apvalerr,winfo,lun_col,runlog,
     & akfile,rayfile,mavfile,targmol,linefiles,parfile,
     & dplist,ntg,speci,nspeci,solarll,pars,sptfile)
c======================================================================
c Compute md5sum checksums of input files and write to "check_md5sums_012345.tmp"
      if(dl.eq.'/')then
          write(csfilename,'(a14,i6.6,a4)') 'check_md5sums_',getpid(),
     &    '.tmp'
          istat=system('md5sum '//dplist//' > '//csfilename)
          istat=system('md5sum '//apvalerr//' >> '//csfilename)
          istat=system('md5sum '//runlog//' >> '//csfilename)
          istat=system('md5sum '//levfile//' >> '//csfilename)
          istat=system('md5sum '//mavfile//' >> '//csfilename)
          istat=system('md5sum '//rayfile//' >> '//csfilename)
          istat=system('md5sum '//parfile//' >> '//csfilename)
          istat=system('md5sum '//window//' >> '//csfilename)
          call substr(linefiles,linelists,9,nss)
          do i=1,nss
              istat=system('md5sum '//linelists(i)//' >> '//csfilename)
          end do
          istat=system('md5sum '//solarll//' >> '//csfilename)

          open(lun_cs, file=csfilename,status='old')
          open(lun_col,file=colfile,status='unknown')
          write(lun_col,'(2i3)')  nlhead_ggg+1+nss,9+4*ntg
          write(lun_col,'(a)') gfit_version
          write(lun_col,'(a)') gsversion(:lnbc(gsversion))
          read(lun_cs,'(a32,2x,a)') checksum,dplist
          write(lun_col,'(a32,2x,a)') checksum,dplist(:lnbc(dplist))
          read(lun_cs,'(a32,2x,a)') checksum,apvalerr
          write(lun_col,'(a32,2x,a)') checksum,apvalerr(:lnbc(apvalerr))
          read(lun_cs,'(a32,2x,a)') checksum,runlog
          write(lun_col,'(a32,2x,a)') checksum,runlog(:lnbc(runlog))
          read(lun_cs,'(a32,2x,a)') checksum,levfile
          write(lun_col,'(a32,2x,a)') checksum,levfile(:lnbc(levfile))
          write(lun_col,'(32x,2x,a)') amodel(:lnbc(amodel))
          write(lun_col,'(32x,2x,a)') vmrset(:lnbc(vmrset))
          read(lun_cs,'(a32,2x,a)') checksum,mavfile
          write(lun_col,'(a32,2x,a)') checksum,mavfile(:lnbc(mavfile))
          read(lun_cs,'(a32,2x,a)') checksum,rayfile
          write(lun_col,'(a32,2x,a)') checksum,rayfile(:lnbc(rayfile))
          read(lun_cs,'(a32,2x,a)') checksum,parfile
          write(lun_col,'(a32,2x,a)') checksum,parfile(:lnbc(parfile))
          read(lun_cs,'(a32,2x,a)') checksum,window
          write(lun_col,'(a32,2x,a)') checksum,window(:lnbc(window))
          do i=1,nss
              read(lun_cs,'(a32,2x,a)') checksum,linelists(i)
              write(lun_col,'(a32,2x,a)')checksum,
     &        linelists(i)(:lnbc(linelists(i)))
          end do
          read(lun_cs,'(a32,2x,a)') checksum,solarll
          write(lun_col,'(a32,2x,a)') checksum,solarll(:lnbc(solarll))
          write(lun_col,'(34x,a)') akfile(:lnbc(akfile))
          write(lun_col,'(34x,a)') sptfile(:lnbc(sptfile))
          write(lun_col,'(34x,a)') colfile(:lnbc(colfile))
          write(lun_col,'(a)') winfo(:lp)
          close(lun_cs,status='delete')
	else
          call substr(linefiles,linelists,9,nss)
          open(lun_col,file=colfile,status='unknown')
          write(lun_col,'(2i3)')  nlhead_ggg+1+nss,9+4*ntg
          write(lun_col,'(a)') gfit_version
          write(lun_col,'(a)') gsversion(:lnbc(gsversion))
          write(lun_col,'(a)') dplist(:lnbc(dplist))
          write(lun_col,'(a)') apvalerr(:lnbc(apvalerr))
          write(lun_col,'(a)') runlog(:lnbc(runlog))
          write(lun_col,'(a)') levfile(:lnbc(levfile))
          write(lun_col,'(a)') amodel(:lnbc(amodel))
          write(lun_col,'(a)') vmrset(:lnbc(vmrset))
          write(lun_col,'(a)') mavfile(:lnbc(mavfile))
          write(lun_col,'(a)') rayfile(:lnbc(rayfile))
          write(lun_col,'(a)') parfile(:lnbc(parfile))
          write(lun_col,'(a)') window(:lnbc(window))
          do i=1,nss
              write(lun_col,'(a)')linelists(i)(:lnbc(linelists(i)))
          end do
          write(lun_col,'(a)') solarll(:lnbc(solarll))
          write(lun_col,'(a)') akfile(:lnbc(akfile))
          write(lun_col,'(a)') sptfile(:lnbc(sptfile))
          write(lun_col,'(a)') colfile(:lnbc(colfile))
          write(lun_col,'(a)') winfo(:lp)

	endif
c
c Do the real spectral fitting.
      call spectrum_loop(apvalerr,winfo,lun_col,runlog,
     & akfile,rayfile,mavfile,targmol,linefiles,parfile,
     & dplist,ntg,speci,nspeci,solarll,pars,sptfile)

      close(lun_col)
      stop
      end
