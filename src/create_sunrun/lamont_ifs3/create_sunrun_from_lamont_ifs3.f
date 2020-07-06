c  Program to create a GGG-compatible sunrun from an ascii list
c  of OPUS spectra.
c
c  To create the input file for this program run list_maker.
c
c  The sunrun program can then be immediately run from there.
c  It will search the local directory for the .gnd file, and
c  write the .gop file to the sunruns/ directory.
c
      implicit none
      include "../../gfit/ggg_int_params.f"
      include "../params.f"

      integer*4 iy,im,id,hh,mm,ss,ms,pkl,prl,gfw,gbw,
     & fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nip,dfr,
     & lr,lrt,ls,lom,doy,warning_flag,
     & possp,istat,object,mcol,ncol,lst,idum
      parameter (mcol=40)
c
      real*8 tins,pins,hins,tout,pout,hout,wspd,wdir,
     & wavtkr,oblat,oblon,obalt,lfl,hfl,foc,
     & fsf,tcorr,nus,nue,lwn,sia,sis,fvsi,vdc,aipl,tel_mag,
     & fxv,lxv,apt,dur,vel,phr,res,pout_cor,lse,lsu,lsf,dip,
     & mvd,snr
c
      character 
     & header*512,outarr(mcol)*20,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & version*64                 !current program version

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)

      warning_flag=0
c
c      pout_cor = 1.0 ! DW090331 this was what the pout_cor was before 090324
c It probably should have been more like -0.89hPa, but the Zeno battery
c was not plugged in, and the values were very noisy.
       pout_cor = 1.1 ! DW090331 this is the pout_cor value after 090324

      version=' create_sunrun     Version 1.3.2     25-Jun-2009     GCT'

      write(6,'(a)') version
      write(6,'(a,f5.2,a)')'Pout will be adjusted by ',pout_cor,' mb'
      write(6,'(a)')'Temperature and humidity defaults used '// 
     & 'for May 2004.'
      col1=' '
      call getendian(iend)  ! iend=+1 on Sun; iend=-1 on PC
c
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)       !Length of root
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

      header=
     &' Spectrum_File_Name                  Obj   tcorr'//
     &'   oblat   oblon   obalt   tins   pins   hins'//
     &'  tout  pout   hout    sia    fvsi   wspd   wdir'//
     &'   Nus    Nue      FSF      lasf     wavtkr'//
     &'   AIPL   TM'
      call substr(header, outarr, mcol,ncol)
      open(luns,file=logfile,status='old')
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
      open(lunt,file=outfile,status='unknown')
      write(lunt,*)3,ncol
      write(lunt,'(a)') version
      write(lunt,'(a)') header(:lnbc(header))
c
      object=2 ! Sun
      tcorr=0.0d0       ! Times are in UT
      fsf=0.99999999d0
      wavtkr=9900.0d0
      aipl=0.002        ! Airmass-Independent Path Length (km)
      tel_mag=1.0
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if(ext(1:1).eq.'g'.and.specname(11:11).ne.'s') then
             write(*,*) ext, specname(:22),
     &      'Warning: .gnd input file contains non-solar spectra'
         endif
         if(ext(1:1).eq.'l'.and.specname(11:11).ne.'l') then
            write(*,*) ext,specname(:22),
     &      'Warning: .lab input file contains non-lab spectra'
         endif

         if(specname(16:16).eq.'a') then      ! InGaAs detector
            nus=4000.0d0
            nue=11000.0d0
         elseif(specname(16:16).eq.'b') then  ! Si detector
            if(specname(15:15).eq.'0') then
               nus=11000.0d0  ! No Red Filter
               nue=13500.0d0  ! No Red Filter
            else
               nus=11000.0d0  ! Red Filter
               nue=15000.0d0  ! Red Filter
            endif
         elseif(specname(16:16).eq.'c') then   ! InSb detector
            nus=1800.0d0
            nue=4000.0d0 ! set to 4000 to avoid overlap with InGaAs
         else
            write(*,*) 'Spectrum Name Length =',ls-1
            write(*,*) 'Unknown spectrum type: '//string
            stop
         endif
c
c  find the spectral file, return the PATH to the spectrum
          dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
          call gindfile(dplist,specname,path)
          if(lnbc(path).le.0) then
            write(6,*) ' Not Found : ',specname
            stop
          endif

          call read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,
     &     id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lwn,foc,nip,dfr,
     &     pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,
     &     tins,pins,hins,tout,pout,hout,wspd,wdir,sia,sis,vdc,
     &     lst,lse,lsu,lsf,dip,mvd,snr)

c  Apply correction to measured surface pressure.
          pout = pout + pout_cor

c  Apply correction for timing errors from July 15 - October 24, 2010
          if(iy.eq.2010 .and. ((im.eq.7 .and. id.ge.15) .or.
     &      im.eq.8 .or. im.eq.9 .or.
     &      (im.eq.10 .and. id.le.24))) then
c            approximate length of the month in days
             lom = 31
c            doy 181 is June 30, 196 is July 15, 297 is October 24
             doy = 181 + (im-7)*lom + id
c            linearly increasing tcorr from 0 to -37 seconds
             tcorr = (doy-196)*37/(196-297)
c            write(*,*)'doy=',doy,'tcorr=',tcorr
          else
             tcorr = 0.0d0
          endif

c  There was a bug in the QNX6 driver that made all tout values positive
c  (i.e. abs(tout)) from September 16, 2013 to February 28, 2014. If
c  these values are ever needed, they should be corrected here.
          if(((iy.eq.2013 .and. ((im.eq.9 .and. id.ge.16) .or.
     &       im.ge.10)) .or.
     &       (iy.eq.2014 .and. im.lt.3)).and.warning_flag.eq.0) then
             write(*,*)'Warning! During the period from '//
     &         'September 16, 2013 '//
     &         'through February 28, 2014, outside '//
     &         'temperatures were recorded without '//
     &         'any negative signs (i.e. they''re all >=0)'
             warning_flag=1
          endif

c         fvsi=sis/sia
          if(sia .le. 0.0 .or. specname(11:11).eq.'l') then
              fvsi=-.9999  ! Missing Value
          else
              fvsi=sis/sia
              if(fvsi.ge.1.0)fvsi=1.
          endif
          call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &    oblon,obalt,tins,pins,hins,tout,pout,hout,sia,fvsi,wspd,wdir,
     &    nus,nue,fsf,lwn,wavtkr,aipl,tel_mag,istat)
c
      if(mod(ispe,1000).eq.0) write(*,*)ispe
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
