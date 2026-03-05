!  Program to create a GGG-compatible sunrun from an ascii list
!  of Wollongong solar spectra, either DA8 or Bruker.
!  Adapted from Bremen version by DG, MAy/June 2007
!  Contains corrections specific to Bruker and DA8 @ Wollongong.
!
!  The .gnd file is intended to reside in the local directory.
!  which is usually the directory containing the spectra.
!  The sunrun program can then be immediately run from there.
!  It will search the local directory for the .gnd file, and
!  write the .gop file to the sunruns/ directory.
!
!  Revisions
!  =========
!   02 Jul 2008 - added date-dependent pressure calibration to Read_Oscar_Log
!   April 2009 - added ipp-generated spectrum files
!   19 May 2009 - PC must use environment variable GGGPATH
!       Determine if PC or UNIX from that:
!       If GGGPATh contains \, assume PC, otherwise (default) it is assumed UNIX
!   March 2011
!   - tidy/unify opus and ipp cases in one version
!   - add wind speed and direction from Oscar log
!       - note this required search for corresponding OPUS file by time when crearing sunrun for IPP spectra
!       - use true Davis outside temperature not room temperature for Tout
!   August 2020
!   - Major re-working for compatability with GGG2020 and reading _sunrun input file
!   - Tried to follow as much as possible the standard GGG create_sunrun format
!   - NMD 20200811

      implicit none

      integer*4 iy,im,id,hh,mm,ss,ms,pkl,prl,gfw,gbw, &
       fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nip,dfr,bytepw, &
       mdet,ndet,idet,    &  ! Number of detectors
       ifirst,ilast,  &
       mpath,mfilepath,nchar,  &
       lr,lrt,ls,lunr,luns,lunt,doy,  &
       possp,istat,object,mcol,ncol,instr,lst
      parameter (lunr=14,luns=15,lunt=16,mcol=40,mdet=4,mpath=128,  &
       mfilepath=mpath+40,nchar=57)   

      integer filter
!    integer*4 i_dl,ierror, ispec, julian
!    integer*4 idate,jd, &
!    nss, &
!    one, &
!    
      integer specmissing, daylogmissing
    
!    character*256   gggpath
!
      real*8 tins,pins,hins,tout,pout,pout_corr,hout,  &
       wspd,wdir,gmt,fovi,opd,snr,asza,delwav,   &
       wavtkr,oblat,oblon,obalt,lfl,hfl,foc,   &
       fsf,tcorr,vdc,lse,lsu,lsf,dip,mvd,  &
       nus(mdet),nue(mdet),   &
       lasf,sia,sis,fvsi,aipl,tel_mag,   &
       fxv,lxv,apt,dur,vel,phr,res,ptrue

      integer k1(mdet), k2(mdet),irec, nhead, lsss

!      real*8 
!       lwn, &
!       pout_cor, &
!     & sigflo(8), sigfhi(8)!,  mode, FileTimeUT, 
!
      character  &
       sss(mdet)*2,   &
       xx_sunrun_file*(mfilepath),   &
       header*512,outarr(mcol)*20,   &
       col1*1,    &                !first column of runlog record
       apf*2,     &                !apodization function (e.g. BX N2, etc)
       dl*1,      &                !forward or backward slash
       ext*3,     &                !geometry ['air','bal','gnd','lab',orb','syn']
       logfile*64,   &
       string*256,   &
       outfile*(mpath+40),   &
       dplist*(mpath+25),   &
       path*(mpath),       &       !ggg directory path (GGGPATH?)
       gggdir*(mpath),     &       !ggg directory path (GGGPATH?)
       specname*(nchar),    &      !spectrum name
       version*64                 !current program version
      
      character  &
       spfmt*2                    !format (OPUS, DA8, etc)
     
!     character apf*2,spfmt*2,listfile*40, &
!    specpath*256,root*94,logname*256, &
!     logpath*256, lowercase*3, logline*1024

!    logical*1   foundinlog

!    DATA sigflo /3950.0, 2800.0, 2350.0, 2050.0, 1850.0,  700.0,  980.0,  700.0/
!    DATA sigfhi /4450.0, 3500.0, 3120.0, 2600.0, 2250.0, 1350.0, 1350.0, 1025.0/

      version=' create_sunrun_from_Wgong     Version 4.0     2020-08-11     NMD'
      write(6,'(a)') version     

      col1=' '
      call getendian(iend)  ! iend=+1 on Sun; iend=-1 on PC
      vdc = 0.d0  ! initialize VDC
      
!    Prompt user for name of input list
      lr=0
      do while(lr.eq.0)
        if (iargc() == 0) then
          write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd):'
          read(*,'(a)') logfile
        elseif (iargc() == 1) then
          call getarg(1, logfile)
        else
          stop 'Usage: $gggpath/bin/create_sunrun pa2004.gnd'
        endif
        lr=lnbc(logfile)
      end do
      ext=logfile(lr-2:lr)
      spfmt=logfile(lr-1:lr)
      
      call get_ggg_environment(gggdir,dl)
      lrt=lnbc(gggdir)        !Length of root
      
!  Read file containing the invariant parameters.
      xx_sunrun_file=gggdir(:lrt)//'tccon/'//logfile(1:2)//'_sunrun.dat'
      write(*,*)'Opening xx_sunrun_file = '//xx_sunrun_file
      open(luns,file=xx_sunrun_file,status='old')

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!c     Root specification:        DG000909
!      call getenv('GGGPATH',root)
!      dl='/'
!      root=root(:lnbc(root))//dl
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      header= &
      ' Spectrum_File_Name                                        Obj'// &
      '  tcorr   oblat    oblon   obalt   tins   pins   hins'//  &
      '  tout  pout   hout    sia    fvsi   wspd   wdir'//  &
      '   Nus    Nue      FSF      lasf    wavtkr   AIPL   TM'
      call substr(header, outarr, mcol,ncol)

      read(luns,*) nhead
      do irec=2,nhead
         read(luns,*)
      end do
      read(luns,*) instr
      read(luns,*) object
      read(luns,*) tcorr
      read(luns,*) oblat
      read(luns,*) oblon
      read(luns,*) obalt
      read(luns,*) tins
      read(luns,*) pins
      read(luns,*) hins
      read(luns,*) tout
      read(luns,*) pout
      read(luns,*) pout_corr
      read(luns,*) hout
      read(luns,*) sia
      read(luns,*) fvsi
      read(luns,*) wspd
      read(luns,*) wdir
      read(luns,*) ndet
      do idet=1,ndet
         read(luns,*) k1(idet),k2(idet),sss(idet),nus(idet),nue(idet)
!         write(*,*) k1(idet),k2(idet),sss(idet),nus(idet),nue(idet)
      end do
      read(luns,*) fsf
      read(luns,*) lasf
      read(luns,*) wavtkr
      read(luns,*) aipl
      read(luns,*) tel_mag 
      close(luns)
      write(*,*) 'Closed xx_sunrun.dat file: '//xx_sunrun_file
      lsss=lnbc(sss(1))

!      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'m4'
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
!      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//spfmt
!  we don't yet implement the different formats

      open(lunt,file=outfile,status='unknown')
      write(lunt,*)3,ncol
      write(lunt,'(a)') version
      write(lunt,'(a)')header(:lnbc(header))

      open(lunr,file=logfile,status='old')
      
      specmissing = 0
!      daylogmissing = 0

      do ispe=1,999999  !--------Main loop over spectra----------
        read(lunr,'(a)',end=99) string
        
!  Remove any leading spaces/blanks.
        string=string(fnbc(string):)
        ls=fbc(string)
        specname=string(:ls-1)
        do idet=1,ndet
          if(specname(k1(idet):k2(idet)).eq.sss(idet)(1:lsss)) exit
        end do

!  find the spectral file, return the PATH to the spectrum
        dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
        call gindfile(dplist,specname,path)
        if(lnbc(path).le.0) then
          write(6,*) ' Not Found : ',specname
          stop
        endif
        
!   Do some instrument specific header reads.
!   Generic version caters for MkIV and OPUS, adding DA8
!   instr==1 is MkIV, 2 is TCCON, setting 8 to Bomem DA8
        if(instr.eq.1) then
          call read_mkiv_header(specname,path,iend,ifirst,ilast,possp, &
         bytepw,apf,delwav,opd,fovi,snr,oblat,oblon,obalt,pout,  &
         tout,hout,asza,iy,im,id,gmt,wavtkr,tins,pins,hins,lasf)

!        elseif(instr.eq.8) then ! Bomem DA8, most header values set manually as follows
!          lasf=15798.012D0     ! Defined by PCAT
!          fsf=1.0D0

!!       Use truncated start and end frequencies to exclude zero signal regions
!          read(specname(1:1),'(i1)')filter
!!          nus=sigflo(filter)
!!          nue=sigfhi(filter)

!!       Get P from header, hardwire T and H
!!          call read_DA8_Press(specname,specpath,pout,pins)
!          tins=22.0
!          tout=22.0
!          hins=5.0
!          hout=50.0
!          fvsi=0.0
!          sia=0.0
!          wspd=0.0
!          wdir=0.0        

        else ! TCCON/Bruker/OPUS
!        write(*,*) 'Calling read_opus_header: path = '//path
          call read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,  &
         id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lasf,foc,nip,dfr,  &
         pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,  &
         tins,pins,hins,tout,pout,hout,wspd,wdir,sia,sis,vdc,  &
         lst,lse,lsu,lsf,dip,mvd,snr)
     
        endif

!  Apply correction to measured surface pressure 
        ptrue = pout + pout_corr  

!  Calculate fvsi as sis/sia if sia is not zero or missing, and if 
!  it's not a lamp run:
          if(sis.le.0.0 .or. sia.le.0.0 .or. ext(:1).eq.'l') then
!  DW 20170809: I've changed sia to sis above, because the EM27/SUN 
!  community uses a computed SIA value from MXY and MNY to track mirror 
!  degredation. They do not compute SIS, so they depend on I2S to compute
!  FVSI from the VDC parameter.
            if(vdc.gt.0) then
              fvsi=vdc
            else
!             Use the missing value in xx_sunrun.dat
            endif
          else
            fvsi=sis/sia
          endif
          if(fvsi.ge.1.0)fvsi=1.
!  accounts for default value being set to something <1

!  If the TCCON spectrum has no red filter, limit the upper
!  frequency (nue) to 13500.0d0
          if(instr.eq.2) then ! TCCON spectrum
            if(specname(15:15).eq.'0') then ! no Si red filter
              nue=min(nue,13500.0d0)
            endif
          endif

!  Site-dependent corrections should go in here. For example, to correct 
!  a timing error in Lamont:
         
!  Apply correction for timing errors from July 15 - October 24, 2010
          if(instr.eq.2 .and. logfile(1:2).eq.'oc') then ! Lamont
            if(iy.eq.2010 .and. ((im.eq.7 .and. id.ge.15) .or. &
             im.eq.8 .or. im.eq.9 .or.  &
             (im.eq.10 .and. id.le.24))) then
!              write(*,*)'Applying Lamont timing error correction.'
!              doy 181 is June 30, 196 is July 15, 297 is October 24
!              31 is approximate length of month in days
               doy = 181 + (im-7)*31 + id
!              linearly decreasing tcorr from 0 to -37 seconds
               tcorr = (doy-196)*37/(196-297)
            else
               tcorr = 0.0d0
            endif
          elseif(instr.eq.2 .and. logfile(1:2).eq.'wg') then ! Wollongong Bruker
            call pressure_calibration(Pout, specname(1:),pout_corr)
            ptrue = pout
          endif

          call write_sunrun(lunt,col1,specname,object,tcorr,oblat,  &
            oblon,obalt,tins,pins,hins,tout,ptrue,hout,sia,fvsi,wspd,wdir,  &
            nus(idet),nue(idet),fsf,lasf,wavtkr,aipl,tel_mag,istat)

        if(mod(ispe,1000).eq.0) then
          write(*,*)ispe
        endif
      end do ! -------------Main loop over spectra----------------

 99   close(lunr)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
