!  Program to create a GGG-compatible sunrun from an ascii list
!  of Wollongong solar spectra, either DA8 or Bruker.
!  Adapted from Bremen version by DG, MAy/June 2007
!  Contains corrections specific to Bruker and DA8 @ Wollongong.
!
!  The .gnd file is intended to reside in the local directory,
!  which is usually the directory containing the spectra.
!  The sunrun program can then be immediately run from there.
!  It will search the local directory for the .gnd file, and
!  write the .gop file to the sunruns/gnd/ directory.
!
!  Revisions
!  =========
!  02 Jul 2008 - added date-dependent pressure calibration to Read_Oscar_Log
!

    implicit none

    integer*4 i,ierror, ispec
    integer*4 idate,iy,im,id,hh,mm,ss,ms,pkl,prl,gbw,gfw, &
    fnbc,lnbc,fbc,iend,dtype,nsp,nip,dfr, &
    lr,lrt,ls,luns,lunt,possp,istat,object, &
    filter,specmissing, logmissing,platform
    
    parameter (luns=15,lunt=16)
!
    real*8 tins,pins,hins,tout,pout,hout, &
    obalt,wavtkr,oblat,oblon,lfl,hfl,foc, &
    fsf,tcorr,nus,nue,lwn,sia,sis,aipl,tel_mag, &
    fxv,lxv,apt,dur,vel,phr,res, &
    sigflo(8), sigfhi(8),fvsi,wspd,wdir
!
    character apf*2,dl*1,ext*3,spfmt*2,listfile*40, &
    outfile*256,specpath*256,root*94,dplist*80,specname*256,logname*256, &
    string*163,user*8,col1*1, logpath*256, lowercase*3

    logical*1   foundinlog

    DATA sigflo /3950.0, 2800.0, 2350.0, 2050.0, 1850.0,  700.0,  980.0,  700.0/
    DATA sigfhi /4450.0, 3500.0, 3120.0, 2600.0, 2250.0, 1350.0, 1350.0, 1025.0/

!   Begin
!   =====
    write(6,*)'SUNRUN     Wollongong Version 3.0    12-Feb-2009  ND'
    col1=' '
    call getendian(iend)  ! iend=+1 on Sun; iend=-1 on PC
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!   Platform specification:      DG000909
!   Also need to select path for m4head.inc in gethead.f
    call getenv('LOGNAME',user)
    if(user.ne.'      ')then
       platform=0           !0=Sun, 1=PC-Linux, 2=PC-Win32
       call getenv('GGGPATH',root)
       dl='/'
       root=root(:lnbc(root))//dl
    else
       platform=2           !0=Sun, 1=PC-Linux, 2=PC-Win32
       dl=char(92)  ! backslash ('\')
       root='g:'//dl
       user='PC-Win'
    endif
    lrt=lnbc(root)     !Length of root
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    lr=0
    do while(lr.eq.0)
         if (iargc() == 0) then
            write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd):'
            read(*,'(a)') listfile
         elseif (iargc() == 1) then
            call getarg(1, listfile)
         else
            stop 'Usage: $gggpath/bin/create_sunrun pa2004.gnd'
         endif

       lr=lnbc(listfile)
    end do
    ext=listfile(lr-2:lr)
    spfmt=listfile(lr-1:lr)
    open(luns,file=listfile,status='old')
!
    object=2 ! Sun
    wavtkr=9900.0d0
    oblat = -34.4061
    oblon = 150.8793
    obalt = 0.030
    aipl=0.005      ! Airmass-Independent Path Length (km)
    tel_mag=1.0
    sia= 0.0
    sis=0.0

!-----Main loop over spectra----------
    specmissing=0
    logmissing=0
    
!    do ispe=1,999999    
!    read(luns,'(a)',end=99) string 
     ispec=0
    do     
      read(luns,'(a)',iostat=ierror) string 
      if(ierror<0)exit      !end of file
      ispec=ispec+1
!      print *, ispec

!     Remove any leading spaces/blanks.
      string=string(fnbc(string):)
      ls=fbc(string)
      specname=string(:ls-1)

      i=index(specname,'\')  ! This allows for directory names to be prefixed to filenames

      if(ispec.eq.1)then
        if(specname(i+1:i+2).eq.'wg')then
          outfile=root(:lrt)//'sunruns'//dl//ext//dl//listfile(:lr-2)//'op'
        elseif(ispec.eq.1.and.lowercase(specname(i+10:i+12)).eq.'spc')then
          outfile=root(:lrt)//'sunruns'//dl//ext//dl//listfile(:lr-2)//'D8'
        else
          stop 'Unrecognised spectrum type'
        endif
        open(lunt,file=outfile,status='unknown')
        write(lunt,*) '3 23'
        write(lunt,'(a)') ' Create_sunrun_from_Wgong Version 3.0 '
        write(lunt,'(a)') ' Spectrum_File_Name             '// &
            ' Obj tcorr   oblat   oblon   obalt   tins   pins'// &
            ' hins    tout    pout     hout   sia    fvsi   wspd   wdir'// &
            '  Nus   Nue    FSF  lasf       wavtkr   AIPL     TM'
      endif

!     Find the spectral file, return the PATH to the spectrum
      dplist=root(:lrt)//'config'//dl//'data_part.lst'
      call gindfile(dplist,specname,specpath)
      if(lnbc(specpath).le.0) then
        write(6,*) ' Not Found : ',specname
        specmissing=specmissing+1
        cycle
      endif

!     Choose if Bruker or Bomem spectrum
      if(specname(i+1:i+2).eq.'wg')then                             ! Bruker spectrum
        call read_opus_header(specpath,iend,dtype,nsp,fxv,lxv,iy, &
          im,id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lwn,foc,nip,dfr, &
          pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt, &
          tins,pins,hins,tout,pout,hout,wspd,wdir,sia,sis)

!       to be included
!       NMD 20090211
        wspd=0.0
        wdir=0.0

        nus=fxv
        nue=lxv
!       lwn is read from OPUS header
        fsf=1.0d0       !0.99999850d0


!       If an OPUS/OSCAR spectrum, get Pout, Tout, SIA & SIS from Oscar log file
!       For ipp spectra, Pout, Tout, SIA & SIS are in the OPUS header already

        if(specname(i+11:i+11)=='_') then
            tcorr=-36000.0d0         ! OPUS times are local, Wollongong is UT+10 hrs
            read(specname(i+5:i+10), '(I6)') idate
            if (idate < 080312)then
                logname=specname(i+5:i+10)//'BrukerHR.csv'
            else
                logname=specname(i+5:i+10)//'Solar.csv'
            endif
            call gindfile(dplist,logname,logpath)
            if(lnbc(logpath).gt.0) then
                call Read_Oscar_log(specname,specpath,logname,logpath,Pout,Tout,SIA,SIS,foundinlog) 
                if(.not.foundinlog)then
!                    print *, 'Create_Sunrun: Spectrum ', trim(specpath), ' not found in log ', trim(logpath)
                    logmissing=logmissing+1
                    cycle
                endif
            else
                print *, 'Create_Sunrun: ', 'Oscar log file ', trim(logname), ' not found'
                cycle
            endif
        else
            tcorr=0.0d0         ! ipp times are UT
        endif
        
        hout = 50.
        if(sia.gt.0.0.and.sis.gt.0.0)then
            fvsi=sis/sia
            if(fvsi.ge.1.0)fvsi=0.999
        else
            fvsi=0.0
        endif

      elseif(lowercase(specname(i+10:i+12)).eq.'spc')then ! DA8 spectrum

        lwn=15798.012D0     ! Defined by PCAT
       fsf=1.0D0

!       Use truncated start and end frequencies to exclude zero signal regions
        read(specname(1:1),'(i1)')filter
        nus=sigflo(filter)
        nue=sigfhi(filter)

!       Get P from header, hardwire T and H
        call read_DA8_Press(specname,specpath,pout,pins)
        tins=22.0
        tout=22.0
        hins=5.0
        hout=50.0
        sia=0.0
        fvsi=0.0

      else
        stop 'Unrecognised spectrum type'
      endif

      call write_sunrun(lunt,col1,specname,object,tcorr,oblat, &
        oblon,obalt,tins,pins,hins,tout,pout,hout,sia,fvsi, &
        wspd,wdir,nus,nue,fsf,lwn,wavtkr,aipl,tel_mag,istat)

        if(mod(ispec,1000).eq.0) write(*,*)ispec

    enddo  ! Main loop over spectra
!
    close(luns)
    close(lunt)
    write(*,*) trim(outfile)
    write(*,*) ispec, ' spectra found'
    write(*,*) specmissing, ' spectra missing'
    write(*,*) logmissing, ' log entries missing'
    stop
    end

!****************************************************************
!   Perform an in-place conversion of string SS to lower-case.
!****************************************************************
    function lowercase(ss)
    implicit none

    integer(4)::    i,ic
    character(*)::  lowercase, ss

    do i=1,len(ss)
      ic=ichar(ss(i:i))
      if(ic.ge.65 .and. ic.le.90) ss(i:i)=char(ic+32)
    end do

    lowercase = ss
    end function lowercase

