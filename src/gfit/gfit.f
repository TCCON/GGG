c  Program GFIT
c  See ggg.history for description of latest changes

      implicit none
c
      integer*4
     & lun_col,lun_iso,
     & istat, lnbc,lp,lc,
     & mmode,nmode,         ! Number of vibrational modes
     & nlhead_ggg,          ! number of header lines in .ggg file
     & kgas,kiso,
     & mtg, ntg, jtg,       ! number of target molecules
     & mspeci,nspeci,jspeci ! number of different species listed in PARFILE
c
      parameter (lun_col=14,lun_iso=18,mspeci=140,mtg=15,mmode=30)

      integer*4
     & targmol(mspeci), ! Group assignment for each specie.
     & speci(mtg),      ! speci # of the parent isotopolog of each target gas
     & dgen(mmode),     ! degeneracy of vibrational modes
     & molewt           ! Molar Mass
c
      character
     & colfile*80,       ! output file containing column amounts
     & akfile*80,        ! output file of averaging kernels
     & sptfile*80,       ! output file of ascii spectral fits
     & mavfile*80,       ! file containing T/P & VMR at user-chosen levels
     & rayfile*80,       ! file of slant paths at user-chosen levels
     & parfile*80,       ! path to "molparam.dat"
     & apvalerr*80,      ! path to a priori variable values and uncertainties
     & llsize*80,        ! path to llsize.dat file
     & linefiles*400,    ! paths to linelists
     & solarll*80,       ! solar linelist
     & version*64,       ! version number
     & gsversion*64,     ! GSETUP version number
     & dplist*80,        ! Data Partition List (e.g. m4part.lst)
     & runlog*80,        ! name of occultation file
     & gasname*8,        ! names of species in PARFILE
c     & tname*8,  ! names of species in PARFILE prefixed by "t"
c     & fullname*10,! full names of species in PARFILE
     & winfo*128,         !  window information (command line)
     & pars(mtg)*9,     ! parameters to fit
c     & fdate*26,         ! current date
c     & getlog*8,         ! investigator
c     & hostnam*12,       ! computer hostname
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
     &   speci_id*24

      data speci/mtg*0/
      data version/' GFIT      Version 4.1.2    6-Aug-2008    GCT  '/

      write(6,*)
      write(6,88)version
 88   format(a)
c
c  Read runlog, model, vmrset & window information from input file (.ggg)

c      go to 10
      read(5,*)nlhead_ggg
      read(5,88)gsversion
      read(5,88)dplist
      read(5,88)apvalerr
      read(5,88)runlog
      read(5,88)levfile
      read(5,88)amodel
      read(5,88)vmrset
      read(5,88)mavfile
      read(5,88)rayfile
      read(5,88)parfile
      read(5,88)window
      read(5,88)llsize
      read(5,88)linefiles
      read(5,88)solarll
      read(5,88)akfile
      read(5,88)sptfile
      read(5,88)colfile
      read(5,88)winfo

      lc=index(winfo,'#')
      if(lc.gt.0) winfo=winfo(:lc-1)
      lp=lnbc(winfo)
      lc=index(winfo,':')
      write(6,88)winfo(:lp)
      call lowercase(winfo)
      call substr(winfo(lc+1:),pars,mtg,ntg)
      if(ntg.gt.mtg) write(*,*)' Warning: Increase MTG to',ntg
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
     &          't'//gasname.eq.pars(jtg)) then
                if(targmol(jspeci).eq.0) targmol(jspeci)=jtg
            endif
         end do
      end do
      read(lun_iso,*,end=77)
      stop ' The number of species listed in PARFILE exceeds MSPECI'
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
      open(lun_col,file=colfile,status='unknown')
      write(lun_col,*)  nlhead_ggg+2,9+4*ntg
      write(lun_col,88) version
c     &//fdate()//getlog()//hostnam(:lnbc(hostnam))  ! cpp
      write(lun_col,88) gsversion(:lnbc(gsversion))
      write(lun_col,88) dplist(:lnbc(dplist))
      write(lun_col,88) apvalerr(:lnbc(apvalerr))
      write(lun_col,88) runlog(:lnbc(runlog))
      write(lun_col,88) levfile(:lnbc(levfile))
      write(lun_col,88) amodel(:lnbc(amodel))
      write(lun_col,88) vmrset(:lnbc(vmrset))
      write(lun_col,88) mavfile
      write(lun_col,88) rayfile
      write(lun_col,88) parfile
      write(lun_col,88) window(:lnbc(window))
      write(lun_col,88) llsize(:lnbc(llsize))
      write(lun_col,88) linefiles(:lnbc(linefiles))
      write(lun_col,88) solarll(:lnbc(solarll))
      write(lun_col,88) akfile(:lnbc(akfile))
      write(lun_col,88) sptfile(:lnbc(sptfile))
      write(lun_col,88) colfile(:lnbc(colfile))
      write(lun_col,88) winfo(:lp)
c
c Do the real spectral fitting.
      call spectrum_loop(apvalerr,winfo,lun_col,runlog,
     & akfile,rayfile,mavfile,targmol,linefiles,parfile,
     & dplist,ntg,speci,nspeci,solarll,pars,sptfile)

      close(lun_col)
      stop
      end
