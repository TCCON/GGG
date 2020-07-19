c  Program GFIT
c  See ggg.history for description of latest changes

      implicit none
      include "ggg_int_params.f"
      include "int_params.f"
c
      integer*4
     & nlhead,ncol,
     & lunw_col,lunw_cbf,lunr_iso,lunr_cs,lunr_apx,
     & j,getpid,nch,lcl,isv,lbf,lnblnk,
     & iptg,ipcl,ipfs,ipsg,ipzo,ipcf,lspmax,
     & istat,lnbc,lp,lw,lc,idum,
     & nmode,           ! Number of vibrational modes
     & nlhead_ggg,      ! number of header lines in .ggg file
     & kgas,kiso,
     & ntg, jtg,        ! number of target molecules
     & ncbf,            ! number of continuum terms (basis functions)
     & nfp,             ! number of fited parameters = ntg+ncbf+n
     & nspeci_iso,jspeci    ! number of different species listed in ISOTOPOLOG.DAT
c
      parameter (lunw_col=14,lunr_cs=15,lunr_iso=17,lunr_apx=18,
     & lunw_cbf=20)

      integer*4
     & targmol(mspeci),      ! Group assignment for each specie.
     & speci(mtg),           ! speci # of the parent isotopolog of each target gas
     & dgen(mvmode),         ! degeneracy of vibrational modes
     & icode,                ! Isotopolog code
     & molewt,               ! Molar Mass
     & system                ! For calling system operations
c
      character
     & version*64,           ! gfit version number
     & rlgfile*(mfilepath),  ! name of occultation file
     & gggdir*(mpath),
     & dl*1,
     & gasname*8,            ! names of species in ISOTOPOLOG.DAT
     & winfile*(mfilepath),  ! name of window list
     & colfile*80,           ! output file containing column amounts
     & colabel*1024,          ! Header labels for .col file
     & akpath*(mfilepath),   ! output file of averaging kernels
     & sptfile*(mfilepath),  ! output file of ascii spectral fits
     & mavfile*80,           ! file containing T/P & VMR at user-chosen levels
     & rayfile*80,           ! file of slant paths at user-chosen levels
     & parfile*(mfilepath),       ! path to "molparam.dat"
     & apvalerr*(mfilepath),      ! path to a priori variable values and uncertainties
     & levfile*(mfilepath),       ! name of file containing fitting levels
     & modelpath*(mfilepath),     ! path to atmospheric model directory 
     & vmrsetpath*(mfilepath),    ! path to vmr set directory
     & speci_id*20,
     & isoformat*80,              ! Format of 'isotopologs.dat' file
     & tll_file*(mfilepath), !
c     & llsize*80,           ! path to llsize.dat file
c     & linefiles*(mfilepath*4),    ! path to linelists
     & solarll*(mfilepath),       ! solar linelist
     & gsversion*64,         ! GSETUP version number
     & dplist*(mfilepath),        ! Data Partition List (e.g. m4part.lst)
c     & tname*8,             ! names of species in ISOTOPOLOG.DAT prefixed by "t"
c     & fullname*10,         ! full names of species in ISOTOPOLOG.DAT
     & winfo*160,            ! window information (command line)
     & pars(mtg)*8,          ! parameters to fit
c     & linelists(9)*(mfilepath),  ! linelists
     & checksum*32,          ! md5sum checksum value
     & csfilename*32,        ! file containing checksums
     & gggfile*80
c
      real*4
     &   atc,                ! Additional T-correction
     &   tdrpf,              ! T-Dependence of Rotational Partition Function
     &   ewvb,       !
     &   apx(mfp),apu(mfp),
     &   fia,delta,lnfrd,  ! Fractional Isotopic Abundance
     &   vibfrq(mvmode)! Array of vibrational frequencies
c
      character
     &   colfile_format*109,
     &   ddd*12,
     &   ss(mfp)*4,
     &   cdum*12

      logical debug,solarll_exist
      character*24 md5sum
      parameter(md5sum = "$GGGPATH/bin/gfit_md5sum")

      data speci/mtg*0/
      
      idum=mauxcol  ! avoid compiler warning (unused variable)
      idum=mcolvav  ! avoid compiler warning (unused variable)
      idum=mgas     ! avoid compiler warning (unused variable)
      idum=mlev     ! avoid compiler warning (unused variable)
      idum=mmp      ! avoid compiler warning (unused variable)
      idum=mrow_qc  ! avoid compiler warning (unused variable)
      idum=ncell    ! avoid compiler warning (unused variable)
      idum=nchar    ! avoid compiler warning (unused variable)

      lspmax=12
      version=
     & ' GFIT                     Version 5.28        2020-04-24   GCT '
      write(6,*)
      write(6,'(a)')version

      winfo=':'
      colabel='Nit  CL    CT    CC   FS    SG    ZO   RMS/CL   Zpres'
      lcl=lnbc(colabel)
      call substr(colabel,cdum,1,nch)
c     write(*,*)'nch=',nch

      colfile_format='(1x,a99,1x,i2,1x,f5.3,1x,f5.1,1x,f4.1,'//
     &'1x,f5.2,'//  ! FS
     &'1x,f5.2,'//  ! SG
     &'1x,f6.4,'//  ! ZO
     &'1x,f6.4,'//  ! RMS/CL
     &'f8.3,15(0pf7.3,1pe11.4,0pf10.5,1pe8.1))'

c     Platform specification:      DG090519
      call get_ggg_environment(gggdir, dl)
c
c  Read runlog, model, vmrset & window information from input file (.ggg)
      call getarg(1, gggfile)
      open(10, file=gggfile, status='old')
      read(10,*)nlhead_ggg
      read(10,'(a)')gsversion
      read(10,'(a)')dplist
      read(10,'(a)')apvalerr
      read(10,'(a)')rlgfile
      read(10,'(a)')levfile
      read(10,'(a)')modelpath
      read(10,'(a)')vmrsetpath
      read(10,'(a)')mavfile
      read(10,'(a)')rayfile
      read(10,'(a)')parfile
      read(10,'(a)')winfile
      read(10,'(a)')tll_file
      read(10,'(a)')solarll
      read(10,'(a)')akpath
      read(10,'(a)')sptfile
      read(10,'(a)')colfile
      do while(winfo(1:1).eq.':')
         read(10,'(a)')winfo
      end do
      close(10)
      akpath=akpath(:lnbc(akpath))//'_'//colfile(:lnbc(colfile)-4)

      if( index(winfo,' sg ').gt.0 .or. index(winfo,' so ').gt.0 .or.
     & index(winfo,' so/').gt.0 ) then
         inquire(file=solarll,exist=solarll_exist)
         if(solarll_exist .eqv. .false.) then
            write(*,*) solarll(:lnblnk(solarll))
            stop 'solar linelist file doesnt exist'
         endif
      endif

      if( index(winfo,' debug ') .gt. 0 ) then
         debug=.true.
      else
         debug=.false.
      endif

c  Look at the WINFO string to index the state vector.
c  The structure of the state vector is as follows. Elements:
c      1 to NTG            are the Target Gases
c  NTG+1 to NTG+NCBF       are the Continuum Basis Functions Coefficients
c  NTG+NCBF+1              is the FS 
c  NTG+NCBF+2              is the ZO
c  
      iptg=0
      ipcl=0
      ipfs=0
      ipsg=0
      ipzo=0
      ipcf=0
      lc=index(winfo,':')
      write(*,*) winfo(:lnbc(winfo))
      call lowercase(winfo)
      call substr(winfo(lc+1:),pars,mtg,ntg)
      if(ntg.gt.mtg) then
         write(*,*)' gfit: Error: NTG > MTG ',ntg,mtg
         stop 'Increase parameter MTG in gfit.f'
      endif
      if(ntg.gt.0) iptg=1

      lw=lnbc(winfo)
      isv=ntg     !   ISV is index into the State Vector (SV)
      if(index(winfo(:lc),' ncbf=').gt.0) then
         lbf = index(winfo(:lc),' ncbf=')
         read(winfo(lbf+6:),*)ncbf
      else
         ncbf=0
         if(index(winfo(:lc),' cl ').gt.0) ncbf=ncbf+1
         if(index(winfo(:lc),' ct ').gt.0) ncbf=ncbf+1
         if(index(winfo(:lc),' cc ').gt.0) ncbf=ncbf+1
      endif
      if(ncbf.gt.0) ipcl=isv+1
      isv=isv+ncbf
      if(index(winfo(:lc),' fs ').gt.0) then
         isv=isv+1
         ipfs=isv
      endif
      if(index(winfo(:lc),' sg ').gt.0) then
         isv=isv+1
         ipsg=isv
      endif
      if(index(winfo(:lc),' zo ').gt.0) then
         isv=isv+1
         ipzo=isv
      endif
      nfp=isv

      if(nfp.gt.mfp) then
         write(*,*) 'mfp,nfp=', mfp,nfp
         stop 'GFIT: NFP > MFP'
      endif

      if(debug) then
         write(*,*)'ip = iptg,   ipcl,   ipfs,    ipsg,   ipzo,     nfp'
         write(*,*)iptg,ipcl,ipfs,ipsg,ipzo,nfp
      endif

      if( index(winfo,' cf ').gt.0) then
         write(*,*)'Fitting channel fringes'
         ipcf=nfp+1
      endif

c----------------------------------------------------------------------
c  Read A Priori Values and Uncertainties. Assign to appropriate SV elements.
c  Read the a priori values and uncertainties of the parameters to be fitted.
      open(lunr_apx,file=apvalerr,status='old')
      read(lunr_apx,*)     ! Skip header line
      isv=ntg
      do j=1,3        ! Read CL, CT, CC a prioris
         if(j.le.ncbf) then
            isv=isv+1
            read(lunr_apx,'(f5.0,f9.0,1x,a4)') apx(isv),apu(isv),ss(isv)
         else
            read(lunr_apx,*)  !  Skip if j>ncbf
         endif
      end do

c  For higher order continuum basis functions, use the CC a priori.
      do j=4,ncbf
         isv=isv+1
         apx(isv)=apx(isv-1)
         apu(isv)=apu(isv-1)
         ss(isv) =ss(isv-1)
      end do
c
c  Read FS priors
      if(ipfs.gt.0) then
         read(lunr_apx,'(f5.0,f9.0,1x,a4)') apx(ipfs),apu(ipfs),ss(ipfs)
      else
         read(lunr_apx,*)                   ! Skip if FS not fitted
      endif
c
c  Read SG priors
      if(ipsg.gt.0) then
         read(lunr_apx,'(f5.0,f9.0,1x,a4)') apx(ipsg),apu(ipsg),ss(ipsg)
      else
         read(lunr_apx,*)                   ! Skip if FS not fitted
      endif
c
c  Read ZO priors
      if(ipzo.gt.0) then
         read(lunr_apx,'(f5.0,f9.0,1x,a4)') apx(ipzo),apu(ipzo),ss(ipzo)
      else
         read(lunr_apx,*)                   ! Skip if ZO not fitted
      endif

      read(lunr_apx,*)                      ! Solar Scaling

      if(ntg.ge.1) read(lunr_apx,*) apx(1),apu(1),ss(1)   ! First Target Gas
      if(ntg.ge.2) read(lunr_apx,*) apx(2),apu(2),ss(2)   ! Other Target Gases
      close(lunr_apx)

c      call vmov(apx(2),0,apx(3),1,ntg-2)   ! absorber amount
c      call vmov(apu(2),0,apu(3),1,ntg-2)   ! absorber uncertainty
      do j=3,ntg
         apx(j)=apx(j-1)
         apu(j)=apu(j-1)
         ss(j)=ss(j-1)
      end do

      if(debug) then
         do isv=1,nfp
            write(*,*) isv,apx(isv),apu(isv),ss(isv)
         end do
      endif
c----------------------------------------------------------
c
c  Read in names of isotopomers
      open(lunr_iso,file=parfile,status='old')
      read(lunr_iso,*) nlhead,ncol
      read(lunr_iso,'(7x,a)') isoformat
      do j=3,nlhead
         read(lunr_iso,*)               ! Column headers
      end do
      do jspeci=1,mspeci
c         write(*,*) 'gfit: calling read_isotopolog: jspeci =',jspeci
         call read_isotopolog(lunr_iso,isoformat,kgas,kiso,gasname,
     &   speci_id,icode,fia,delta,lnfrd,molewt,ewvb,atc,tdrpf,vibfrq,
     &   dgen,nmode,mvmode,istat)
c         write(*,*) 'jspeci,istat = ',jspeci,istat
         if(istat.ne.0) exit
         call lowercase(gasname)
         targmol(jspeci)=0
         do jtg=1,ntg
            if(kiso.le.9) then
               write(ddd,'(i1,a)')kiso,gasname
            else
               write(ddd,'(i2,a)')kiso,gasname
            endif
c            write(*,*)'jtg,pars(jtg),ddd = ',jtg,pars(jtg),ddd
            if(ddd.eq.pars(jtg)) targmol(jspeci)=jtg
            if( gasname//' '.eq.pars(jtg) .or. 
     &         'a'//gasname.eq.pars(jtg) .or.
     &         'b'//gasname.eq.pars(jtg) .or.
     &         'f'//gasname.eq.pars(jtg) .or. ! fringes from fco2
     &         'z'//gasname.eq.pars(jtg) .or. ! zero level offset from co2
     &         'w'//gasname.eq.pars(jtg) .or. ! weak co2
     &         'l'//gasname.eq.pars(jtg) .or. ! strong co2 (using 's' would confuse so2)
     &         't'//gasname.eq.pars(jtg)) then
               if(targmol(jspeci).eq.0) targmol(jspeci)=jtg
            endif
c            write(*,*)'gfit: jspeci,jtg,targmol(jspeci)=',jspeci,jtg,
c     &      targmol(jspeci)
         end do    !  do jtg=1,ntg
      end do    !  do jspeci=1,mspeci
      close(lunr_iso)
      nspeci_iso=jspeci-1
      write(*,*)'nspeci_iso=',nspeci_iso
c
      do jspeci=nspeci_iso,1,-1
         if(targmol(jspeci).gt.0) speci(targmol(jspeci))=jspeci
c         write(*,*)'gfit: jspeci,targmol,speci(targmol)=',
c     &   jspeci,targmol(jspeci),speci(targmol(jspeci))
      end do
c
      do jtg=1,ntg
         if(speci(jtg).eq.0) then
            write(*,*) jtg,pars(jtg),' unrecognized/duplicated'
            stop
         endif
         lp=lnbc(pars(jtg))
         colabel=colabel(:lnbc(colabel))//
     &   ' AM_'//pars(jtg)(:lp)//
     &   ' OVC_'//pars(jtg)//' VSF_'//pars(jtg)(:lp)//
     &   ' VSF_'//pars(jtg)(:lp)//'_error'
      end do
c========================================================================
c The first time that spectrum_loop is called it skips the
c CPU intensive operations (e.g. computing VACs, spectral fitting)
c and simply checks that the spectra, .mav, and .ray files are okay.
c This prevents GFIT spending 20 minutes processing and then crashing.
      write(*,*) ' Pre-screening input files...'
      write(*,*) ' Calling spectrum loop: nspeci_iso=',nspeci_iso

      call spectrum_loop(winfo,debug,
     & lunw_col,lunw_cbf,lcl,colabel,colfile_format,lspmax,
     & rlgfile,akpath,rayfile,mavfile,targmol,tll_file,parfile,
     & apx,apu,dplist,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & ntg,ncbf,nfp,speci,nspeci_iso,solarll,pars,sptfile)
c      write(*,*) ' Exited spectrum loop'
      write(colfile_format(6:7),'(i2.2)') lspmax
c======================================================================
c Compute md5sum checksums of input files and write to "check_md5sums_012345.tmp"
      if(dl.eq.'/')then
         write(csfilename,'(a14,i6.6,a4)') 'check_md5sums_',getpid(),
     &   '.tmp'
c         istat=system('md5sum '//dplist//' > '//csfilename)
c         istat=system('md5sum '//apvalerr//' >> '//csfilename)
c         istat=system('md5sum '//rlgfile//' >> '//csfilename)
c         istat=system('md5sum '//levfile//' >> '//csfilename)
c         istat=system('md5sum '//mavfile//' >> '//csfilename)
c         istat=system('md5sum '//rayfile//' >> '//csfilename)
c         istat=system('md5sum '//parfile//' >> '//csfilename)
c         istat=system('md5sum '//winfile//' >> '//csfilename)
         istat=system(md5sum//' '//dplist//' > '//csfilename)
         istat=system(md5sum//' '//apvalerr//' >> '//csfilename)
         istat=system(md5sum//' '//rlgfile//' >> '//csfilename)
         istat=system(md5sum//' '//levfile//' >> '//csfilename)
         istat=system(md5sum//' '//mavfile//' >> '//csfilename)
         istat=system(md5sum//' '//rayfile//' >> '//csfilename)
         istat=system(md5sum//' '//parfile//' >> '//csfilename)
         istat=system(md5sum//' '//winfile//' >> '//csfilename)
         istat=system(md5sum//' '//tll_file//' >> '//csfilename)
c         call substr(linefiles,linelists,9,nss)
c         do i=1,nss
c             istat=
c     &       system(md5sum//' '//linelists(i)//' >> '//csfilename)
c         end do
         istat=system(md5sum//' '//solarll//' >> '//csfilename)
        
         open(lunw_cbf,file=colfile(:lnbc(colfile)-4)//'.cbf',
     &   status='unknown')
         write(lunw_cbf,*) 9,ncbf+4
         write(lunw_cbf,'(a)') 'Channel Fringe info and  Continuum Basis
     & Function Coefficients'
         write(lunw_cbf,'(a)')'CF_Amp is divided by continuum level'
         write(lunw_cbf,'(a)')'CF_Period is in units of cm-1'
         write(lunw_cbf,'(a)')'CF_Phase is in units of rads'
         write(lunw_cbf,'(a)')'CBF_1 is the continuum level   (CL)'
         write(lunw_cbf,'(a)')'CBF_2 is the continuum tilt     (CT)'
         write(lunw_cbf,'(a)')'CBF_3 is the continuum curvature (CC)'
         write(lunw_cbf,'(a,30(a5,i2.2,2x))') 
     &   '     Spectrum_Name    CF_Amp/CL CF_Period CF_Phase ',
     &   (' CBF_',j,j=1,ncbf)

         open(lunr_cs, file=csfilename,status='old')
         open(lunw_col,file=colfile,status='unknown')
         write(lunw_col,'(2i3)')  nlhead_ggg+3,nch+1+4*ntg
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
         read(lunr_cs,'(a32,2x,a)') checksum,tll_file
         write(lunw_col,'(a32,2x,a)')checksum,tll_file(:lnbc(tll_file))
         read(lunr_cs,'(a32,2x,a)')checksum,solarll
         write(lunw_col,'(a32,2x,a)')checksum,solarll(:lnbc(solarll))
         write(lunw_col,'(34x,a)')akpath(:lnbc(akpath))
         write(lunw_col,'(34x,a)')sptfile(:lnbc(sptfile))
         write(lunw_col,'(34x,a)')colfile(:lnbc(colfile))
         write(lunw_col,'(34x,a)')colfile_format(:lnbc(colfile_format))
         write(lunw_col,'(a)')winfo(:lw)
         close(lunr_cs,status='delete')
      else
c         call substr(linefiles,linelists,9,nss)
         open(lunw_col,file=colfile,status='unknown')
         write(lunw_col,'(2i3)')  nlhead_ggg+3,9+4*ntg
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
         write(lunw_col,'(a)') tll_file(:lnbc(tll_file))
         write(lunw_col,'(a)') solarll(:lnbc(solarll))
         write(lunw_col,'(a)') akpath(:lnbc(akpath))
         write(lunw_col,'(a)') sptfile(:lnbc(sptfile))
         write(lunw_col,'(a)') colfile(:lnbc(colfile))
         write(lunw_col,'(a)') colfile_format(:lnbc(colfile_format))
         write(lunw_col,'(a)') winfo(:lw)
      endif
c
c Do the real spectral fitting.
      call spectrum_loop(winfo,debug,
     & lunw_col,lunw_cbf,lcl,colabel,colfile_format,lspmax,
     & rlgfile,akpath,rayfile,mavfile,targmol,tll_file,parfile,
     & apx,apu,dplist,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & ntg,ncbf,nfp,speci,nspeci_iso,solarll,pars,sptfile)

      close(lunw_col)
      close(lunw_cbf)
c     stop
      end
