c  Program to create a GGG-compatible sunrun from an ascii list
c  of OPUS spectra.
c
c  To create the input file for this program run list_maker.
c
c  The sunrun program can then be immediately run from there.
c  It will search the local directory for the .gnd file, and
c  write the .gop file to the sunruns/ directory.
c
c input
c   1.  Spectrum name (up to 21 characters)
c   2. Source object (1=Moon; 2=sun)
c   3. Time correction (seconds)
c   4. Observation Latitude (deg)
c   5. Observation Longitude (deg)
c   6. Observation Altitude (km)
c   7. Internal Temperature (C)
c   8. Internal Pressure (mbar)
c   9. Internal humidity (%)
c  10. External Temperature (C)
c  11. External Pressure (mbar)
c  12. External Humidity (%) 
c  13. Solar Intensity Average (Arbitrary Units)
c  14. Fractional Variation in Solar Intensity (SIS/SIA)
c  15. Wind Speed (m/s)
c  16. Wind Direction (deg.)
c  17. Useful Starting frequency (cm-1)
c  18. Useful Ending Frequency (cm-1)
c  19. Frequency Scale Factor (dimensionless)
c  20. Laser Frequency (cm-1)
c  21. Suntracker Center Frequency (cm-1) 
c  22. Airmass Independent Path Length (km)
c  23. Telescope Magnification (dimensionless)
c output
c format(a1,a38,1x,i2,f8.0,f9.4,f10.4,f7.3,f6.1,f8.2,2f6.1,f8.2,f6.1,f7.1,f7.4,f6.1,f6.0,1x,2f7.0,f11.8,f11.3,f7.0,f7.3,f6.2)
      implicit none
      include '../../gfit/ggg_int_params.f'
      include '../params.f'

      integer*4 iy,im,id,hh,mm,ss,ms,pkl,prl,gfw,gbw,
     & fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nip,dfr,idum,
     & lr,lrt,ls,possp,istat,object,mcol,ncol,kfail,lst
      parameter (mcol=40)
c
      real*8 tins,pins,hins,tout,pout,hout,wspd,wdir,
     & wavtkr,oblat,oblon,obalt,lfl,hfl,foc,
     & fsf,tcorr,nus,nue,lwn,sia,sis,fvsi,aipl,tel_mag,
     & fxv,lxv,apt,dur,vel,phr,res,vdc,lse,lsu,lsf,
     & dip,mvd,snr
c
      character apf*2,dl*1,ext*3,
     & root*64,specname*(nchar),
     & col1*1,version*56,header*512,outarr(mcol)*20

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
c
      version=' create_sunrun     Version 1.33      2016-04-02     GCT'

      write(6,'(a)') version
      col1=' '
      call getendian(iend)  ! iend=+1 on Sun; iend=-1 on PC
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     Root specification:        DG000909
      call getenv('GGGPATH',root)
      dl='/'
      root=root(:lnbc(root))//dl
      lrt=lnbc(root)       !Length of root
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      lr=0
      do while(lr.eq.0)
         write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd): '
         read(*,'(a)') logfile
         lr=lnbc(logfile)
      end do
      ext=logfile(lr-2:lr)

      header=
     &' Spectrum_File_Name                  Obj   tcorr'//
     &'   oblat   oblon   obalt   tins   pins   hins'//
     &'  tout  pout   hout   sia    fvsi  wspd   wdir'//
     &'   Nus    Nue      FSF      lasf     wavtkr'//
     &'   AIPL   TM'
      call substr(header, outarr, mcol,ncol)
      open(luns,file=logfile,status='old')
      outfile=root(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
      open(lunt,file=outfile,status='unknown')
      write(lunt,*)3,ncol
      write(lunt,'(a)') version
      write(lunt,'(a)') header(:lnbc(header))
c
      object=2 ! Sun
      oblat=53.1036d0
      oblon=8.8495d0
      obalt=0.027d0
      tcorr=0.0d0       ! Times are in UT
      fsf=0.99999900
      wavtkr=9900.0d0
      aipl=0.004        ! Airmass-Independent Path Length (km)
      tel_mag=1.0
      kfail=0
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 
         tcorr = 0  ! initialise for each spectrum
         pout = -999. ! also initialise

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if(specname(9:10).eq.'.1' .or. specname(15:15).eq.'p') then
            nus=500.0d0
            nue=1600.0d0
         elseif(specname(9:10).eq.'.2' .or. specname(15:15).eq.'q') then
            nus=1700.0d0
            nue=2400.0d0
         elseif(specname(9:10).eq.'.3' .or. specname(15:15).eq.'r') then
            nus=1800.0d0
            nue=2800.0d0
         elseif(specname(9:10).eq.'.4' .or. specname(15:15).eq.'s') then
            nus=2300.0d0
            nue=3300.0d0
         elseif(specname(9:10).eq.'.5' .or. specname(15:15).eq.'t') then
            nus=2700.0d0
            nue=4000.0d0
         elseif(specname(9:10).eq.'.6' .or. specname(15:15).eq.'u') then
            nus=3800.0d0
            nue=4600.0d0
         elseif(specname(10:11).eq.'.7'.or.specname(9:10).eq.'.7') then
            nus=4000.0d0
            nue=11000.0d0
         elseif(specname(11:12).eq.'sa') then
            nus=4000.0d0
            nue=11000.0d0
         else
            nus = 0.0d0
            nue = 15000.0d0
         endif
 
c         if(specname(11:11).eq.'7') then      ! InGaAs detector
c            nus=4000.0d0
c            nue=11000.0d0
c            write(*,*) 'Spectrum Name Length=',ls-1
c            write(*,*) 'Unkown spectrum type: '//string
c         elseif(specname(12:16).eq.'ghost') then   ! InSb detector
c            nus=4000.0d0
c            nue=11000.0d0
c         else
c            write(*,*) 'Spectrum Name Length =',ls-1
c            write(*,*) 'Unknown spectrum type: '//string
c            stop
c         endif
c
c  find the spectral file, return the PATH to the spectrum
          dplist=root(:lrt)//'config'//dl//'data_part.lst'
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

          if(iy.eq.2005) then
            if(im.eq.3 .and. id.eq.31) then
              tcorr=-7200
            elseif(im.gt.3 .and. im.lt.10) then
              tcorr=-7200
            elseif(im.eq.10 .and. id.lt.19) then
              tcorr=-7200
            elseif(im.gt.9) then
              tcorr=-3600
            endif
          elseif(iy.eq.2006) then
            if(im.lt.4) then
              tcorr=-3600
            elseif(im.lt.10) then
              tcorr=-7200
            elseif(im.eq.10 .and. id.lt.16) then
              tcorr=-7200
            endif
          endif

!          if(specname(11:11).eq.'s') then

!          else
!  the if statement around the following call removes the need
!   to call this routine if the pressure is already in the spectrum
!   header. It is commented out so that pressure correction can be
!   applied consistently if necessary.
!          if(pout.le.0) then
            call get_bremen_log_vals(specname,iend,iy,im,id,
     &       hh,mm,ss,ms,dur,tout,pout,hout)
!          endif

           if(pins.gt.500.0d0) pins=pout

          if(vdc .ge. 0) then
              fvsi=vdc  ! Use DC igm variation
          elseif(sia .le. 0.0 .or. specname(11:11).eq.'l') then
              fvsi=1.0  ! Missing Value
          else
              fvsi=sis/sia
          endif
          if(fvsi.ge.1.0)fvsi=1.

          if(pout.gt.100) then
            call write_sunrun(lunt,col1,specname,object,
     &      tcorr,oblat,oblon,obalt,
     &      tins,pins,hins,tout,
     &      pout,hout,sia,fvsi,wspd,
     &      wdir,nus,nue,fsf,
     &      lwn,wavtkr,aipl,
     &      tel_mag,istat)
          else
            kfail=kfail+1
            ! skip
          endif

      if(mod(ispe,1000).eq.0) write(*,*)ispe
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      write(*,*) kfail, ' skipped'
      stop
      end
