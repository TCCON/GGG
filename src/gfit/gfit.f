C  Program GFIT
c  See ggg.history for description of latest changes

      implicit none
      include "../ggg_int_params.f"
      include "int_params.f"
c
      integer*4
     & lunw_col,lunr_iso,lunr_cs,nss,i,getpid,nch,lcl,
     & istat, lnbc,lp,lc,
     & nmode,               ! Number of vibrational modes
     & nlhead_ggg,          ! number of header lines in .ggg file
     & kgas,kiso,
     & ntg, jtg,            ! number of target molecules
     & nspeci,jspeci        ! number of different species listed in ISOTOPOLOG.DAT
c
      parameter (lunw_col=14,lunr_cs=15,lunr_iso=17)

      integer*4
     & targmol(mspeci), ! Group assignment for each specie.
     & speci(mtg),      ! speci # of the parent isotopolog of each target gas
     & dgen(mvmode),    ! degeneracy of vibrational modes
     & molewt,          ! Molar Mass
     & system           ! For calling system operations
c
      character
     & version*64,       ! gfit version number
     & rlgfile*(mfilepath),       ! name of occultation file
     & gggdir*(mpath),
     & dl*1,
     & gasname*8,        ! names of species in ISOTOPOLOG.DAT
     & winfile*(mfilepath),       ! name of window list
     & colfile*80,       ! output file containing column amounts
     & colabel*600,      ! Header labels for .col file
     & akfile*(mfilepath),        ! output file of averaging kernels
     & sptfile*(mfilepath),       ! output file of ascii spectral fits
     & mavfile*80,       ! file containing T/P & VMR at user-chosen levels
     & rayfile*80,       ! file of slant paths at user-chosen levels
     & parfile*80,       ! path to "molparam.dat"
     & apvalerr*(mfilepath),      ! path to a priori variable values and uncertainties
     & levfile*(mfilepath),       ! name of file containing fitting levels
     & modelpath*(mfilepath),     ! path to atmospheric model directory 
     & vmrsetpath*(mfilepath),    ! path to vmr set directory
     & speci_id*24,
     &
c     & llsize*80,       ! path to llsize.dat file
     & linefiles*(mfilepath*4),    ! path to linelists
     & solarll*(mfilepath),       ! solar linelist
     & gsversion*64,     ! GSETUP version number
     & dplist*(mfilepath),        ! Data Partition List (e.g. m4part.lst)
c     & tname*8,         ! names of species in ISOTOPOLOG.DAT prefixed by "t"
c     & fullname*10,     ! full names of species in ISOTOPOLOG.DAT
     & winfo*160,        ! window information (command line)
     & pars(mtg)*9,      ! parameters to fit
     & linelists(9)*(mfilepath),  ! linelists
     & checksum*32,      ! md5sum checksum value
     & csfilename*32     ! file containing checksums
c
      real*4
     &   tdrpf,        ! T-Dependence of Rotational Partition Function
     &   fia,delta,epsilon,  ! Fractional Isotopic Abundance
     &   vibfrq(mvmode)! Array of vibrational frequencies
c
      character
     &   colfile_format*86,
     &   cdum*12

      character*24 md5sum
      parameter(md5sum = "$GGGPATH/bin/gfit_md5sum")

      data speci/mtg*0/
      
      version=
     & ' GFIT                     Version 4.8.3    09-Jan-2012    GCT '
      write(6,*)
      write(6,'(a)')version

      colabel='Nit  CL   CT   CC   FS   S-G  ZO   RMS/CL   Zpres'
      lcl=lnbc(colabel)
      call substr(colabel,cdum,1,nch)
c     write(*,*)'nch=',nch

      colfile_format='(1x,a,1x,i2,1x,f5.3,4(1x,f4.1),1x,f5.3,1x,f6.4,'//
     &'f8.3,15(0pf7.3,1pe11.4,0pf9.4,1pe8.1))'

c     Platform specification:      DG090519
      call get_ggg_environment(gggdir, dl)
c
c  Read runlog, model, vmrset & window information from input file (.ggg)

      read(5,*)nlhead_ggg
      read(5,'(a)')gsversion
      read(5,'(a)')dplist
      read(5,'(a)')apvalerr
      read(5,'(a)')rlgfile
      read(5,'(a)')levfile
      read(5,'(a)')modelpath
      read(5,'(a)')vmrsetpath
      read(5,'(a)')mavfile
      read(5,'(a)')rayfile
      read(5,'(a)')parfile
      read(5,'(a)')winfile
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
      open(lunr_iso,file=parfile,status='old')
      do jspeci=1,mspeci
         call read_isotop(lunr_iso,kgas,kiso,gasname,speci_id,
     &   fia,delta,epsilon,molewt,tdrpf,vibfrq,dgen,nmode,mvmode,istat)
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
      read(lunr_iso,*,end=77)
      stop ' GFIT: Number of species in ISOTOPOLOG.DAT exceeds MSPECI'
77    close(lunr_iso)
      nspeci=jspeci-1
c
      do jspeci=nspeci,1,-1
         if(targmol(jspeci).gt.0) speci(targmol(jspeci))=jspeci
      end do
c
c========================================================================
c The first time that spectrum_loop is called it skips the
c CPU intensive operations (e.g. computing VACs, spectral fitting)
c and simply checks that the spectra, .mav, and .ray files are okay.
c This prevents GFIT spending 20 minutes processing and then crashing.
      write(*,*) ' Pre-screening input files...'
      call spectrum_loop(apvalerr,winfo,
     & lunw_col,lcl,colabel,colfile_format,
     & rlgfile,
     & akfile,rayfile,mavfile,targmol,linefiles,parfile,
     & dplist,ntg,speci,nspeci,solarll,pars,sptfile)
c======================================================================
c Compute md5sum checksums of input files and write to "check_md5sums_012345.tmp"
      if(dl.eq.'/')then
          write(csfilename,'(a14,i6.6,a4)') 'check_md5sums_',getpid(),
     &    '.tmp'
c          istat=system('md5sum '//dplist//' > '//csfilename)
c          istat=system('md5sum '//apvalerr//' >> '//csfilename)
c          istat=system('md5sum '//rlgfile//' >> '//csfilename)
c          istat=system('md5sum '//levfile//' >> '//csfilename)
c          istat=system('md5sum '//mavfile//' >> '//csfilename)
c          istat=system('md5sum '//rayfile//' >> '//csfilename)
c          istat=system('md5sum '//parfile//' >> '//csfilename)
c          istat=system('md5sum '//winfile//' >> '//csfilename)
           istat=system(md5sum//' '//dplist//' > '//csfilename)
           istat=system(md5sum//' '//apvalerr//' >> '//csfilename)
           istat=system(md5sum//' '//rlgfile//' >> '//csfilename)
           istat=system(md5sum//' '//levfile//' >> '//csfilename)
           istat=system(md5sum//' '//mavfile//' >> '//csfilename)
           istat=system(md5sum//' '//rayfile//' >> '//csfilename)
           istat=system(md5sum//' '//parfile//' >> '//csfilename)
           istat=system(md5sum//' '//winfile//' >> '//csfilename)

          call substr(linefiles,linelists,9,nss)
          do i=1,nss
c              istat=system('md5sum '//linelists(i)//' >> '//csfilename)
              istat=
     &        system(md5sum//' '//linelists(i)//' >> '//csfilename)
          end do
c          istat=system('md5sum '//solarll//' >> '//csfilename)
           istat=system(md5sum//' '//solarll//' >> '//csfilename)
        
          open(lunr_cs, file=csfilename,status='old')
          open(lunw_col,file=colfile,status='unknown')
          write(lunw_col,'(2i3)')  nlhead_ggg+2+nss,nch+1+4*ntg
          write(lunw_col,'(a)') version
          write(lunw_col,'(a)') gsversion(:lnbc(gsversion))
          read(lunr_cs,'(a32,2x,a)') checksum,dplist
          write(lunw_col,'(a32,2x,a)') checksum,dplist(:lnbc(dplist))
          read(lunr_cs,'(a32,2x,a)') checksum,apvalerr
          write(lunw_col,'(a32,2x,a)')checksum,apvalerr(:lnbc(apvalerr))
          read(lunr_cs,'(a32,2x,a)') checksum,rlgfile
          write(lunw_col,'(a32,2x,a)') checksum,rlgfile(:lnbc(rlgfile))
          read(lunr_cs,'(a32,2x,a)') checksum,levfile
          write(lunw_col,'(a32,2x,a)') checksum,levfile(:lnbc(levfile))
          write(lunw_col,'(32x,2x,a)') modelpath(:lnbc(modelpath))
          write(lunw_col,'(32x,2x,a)') vmrsetpath(:lnbc(vmrsetpath))
          read(lunr_cs,'(a32,2x,a)') checksum,mavfile
          write(lunw_col,'(a32,2x,a)') checksum,mavfile(:lnbc(mavfile))
          read(lunr_cs,'(a32,2x,a)') checksum,rayfile
          write(lunw_col,'(a32,2x,a)') checksum,rayfile(:lnbc(rayfile))
          read(lunr_cs,'(a32,2x,a)') checksum,parfile
          write(lunw_col,'(a32,2x,a)') checksum,parfile(:lnbc(parfile))
          read(lunr_cs,'(a32,2x,a)') checksum,winfile
          write(lunw_col,'(a32,2x,a)') checksum,winfile(:lnbc(winfile))
          do i=1,nss
              read(lunr_cs,'(a32,2x,a)') checksum,linelists(i)
              write(lunw_col,'(a32,2x,a)')checksum,
     &        linelists(i)(:lnbc(linelists(i)))
          end do
          read(lunr_cs,'(a32,2x,a)')checksum,solarll
          write(lunw_col,'(a32,2x,a)')checksum,solarll(:lnbc(solarll))
          write(lunw_col,'(34x,a)')akfile(:lnbc(akfile))
          write(lunw_col,'(34x,a)')sptfile(:lnbc(sptfile))
          write(lunw_col,'(34x,a)')colfile(:lnbc(colfile))
          write(lunw_col,'(34x,a)')colfile_format(:lnbc(colfile_format))
          write(lunw_col,'(a)')winfo(:lp)
          close(lunr_cs,status='delete')
        else
          call substr(linefiles,linelists,9,nss)
          open(lunw_col,file=colfile,status='unknown')
          write(lunw_col,'(2i3)')  nlhead_ggg+1+nss,9+4*ntg
          write(lunw_col,'(a)') version
          write(lunw_col,'(a)') gsversion(:lnbc(gsversion))
          write(lunw_col,'(a)') dplist(:lnbc(dplist))
          write(lunw_col,'(a)') apvalerr(:lnbc(apvalerr))
          write(lunw_col,'(a)') rlgfile(:lnbc(rlgfile))
          write(lunw_col,'(a)') levfile(:lnbc(levfile))
          write(lunw_col,'(a)') modelpath(:lnbc(modelpath))
          write(lunw_col,'(a)') vmrsetpath(:lnbc(vmrsetpath))
          write(lunw_col,'(a)') mavfile(:lnbc(mavfile))
          write(lunw_col,'(a)') rayfile(:lnbc(rayfile))
          write(lunw_col,'(a)') parfile(:lnbc(parfile))
          write(lunw_col,'(a)') winfile(:lnbc(winfile))
          do i=1,nss
              write(lunw_col,'(a)')linelists(i)(:lnbc(linelists(i)))
          end do
          write(lunw_col,'(a)') solarll(:lnbc(solarll))
          write(lunw_col,'(a)') akfile(:lnbc(akfile))
          write(lunw_col,'(a)') sptfile(:lnbc(sptfile))
          write(lunw_col,'(a)') colfile(:lnbc(colfile))
          write(lunw_col,'(a)') colfile_format(:lnbc(colfile_format))
          write(lunw_col,'(a)') winfo(:lp)
        endif
c
c Do the real spectral fitting.
      call spectrum_loop(apvalerr,winfo,
     & lunw_col,lcl,colabel,colfile_format,
     & rlgfile,akfile,rayfile,mavfile,targmol,linefiles,parfile,
     & dplist,ntg,speci,nspeci,solarll,pars,sptfile)

      close(lunw_col)
      stop
      end
