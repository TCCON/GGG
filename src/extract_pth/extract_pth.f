c  Program extract_pth
c  Interpolates into a NCEP model files to determine
c  the T/P/H at the altitude of the site, and then output these
c  values along with the corresponding values from the runlog.
c  Useful for comparing the site weather station data with NCEP.

      implicit none
      integer*4 i,iwas,lunm,lunw,ncol,nlhead,lnbc,
     & jul,iyyyy,imm,idd,lrt,ldot,
     & lunr_rlg,             ! Logical unit number
     & istat,             ! status flag (0=success, 1=EOF)
     & iyr,               ! year
     & idoy,              ! day of year
     & ifirst,            ! index of first spectral point in disk file
     & ilast,             ! index of last spectral point in disk file
     & possp,             ! Length of attached header in bytes
     & bytepw,            ! Bytes per data word, usually 2 (I*2) or 4 (R*4)
     & iline              ! loop index for .mod file header

      real*8
     & frdoy, frwas,diff,
     & oblat, latwas,     ! observation latitude (deg).
     & oblon, lonwas,     ! observation longitude (deg).
     & obalt,             ! observation altitude (km)
     & asza,              ! astronomical solar zenith angle (unrefracted)
     & opd,               ! Optical path difference (cm) of interferogram
     & graw,              ! spacing of raw spectrum (cm-1) from GETINFO
     & azim,              ! Solar Azimuth Angle
     & osds,              ! Observer-Sun Doppler Stretch (ppm_
     & fovi,              ! Internal angular diameter of FOV (radians)
     & fovo,              ! External angular diameter of FOV (radians)
     & amal,              ! angular misalignment of interferometer
     & zoff,              ! Zero level offset (dimensionless fraction)
     & snr,               ! Signal-to-Noise Ratio (dimensionless)
     & tins,              ! Inside temperature
     & pins,              ! Inside pressure
     & hins,              ! Inside humidity
     & tout,              ! Outside temperature
     & pout,              ! Outside pressure
     & hout,              ! Outside humidity
     & sia,               ! Solar Intensity Average (arbitrary units)
     & fvsi,              ! Fractional Variation in Solar Intensity
     & wspd,              ! Wind Speed (m/s)
     & wdir,              ! Wind Direction (deg)
     & zpdtim,            ! Time of ZPD (UT hours)
     & zenoff,            ! Zenith angle pointing offset (deg)
     & aipl,              ! Airmass-Independent Path Length (km)
     & lasf,              ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr             ! suntracker frequency (active tracking)

      character
     & col1*1,            ! first column of runlog record
     & specname*57,       ! spectrum name
     & gggdir*256,
     & dl*1,
     & ext*3,
     & data_fmt_read_rl*256,
     & col_labels_rl*320,
     & apf*2,             ! apodization function (e.g. BX N2, etc)
     & modtype*4          ! which met model type ("NCEP" or "FPIT")

      parameter (lunr_rlg=14,lunm=15,lunw=16)
      character modname*80,modname_old*80,vmrname*80,ns*1,ew*1,
     &   version*50,runlog*80
     
      real*4 roc,ecc2,alat,gs,zz,pfact,ptrop,pmod,tmod,hmod,
     &   scht,log1pxox,wexler_svp,svp_wv_over_ice,tout_k,hout_vmr,
     &   p0,t0,z0,m0,h0,p1,t1,z1,m1,h1

      logical read_mod    ! controls whether a new model file needs read

      version=' extract_pth     Version 1.14   2013-03-01    GCT '

      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir

c -------- Get the runlog and model type from the user --------
c May either specify the runlog and y or n for FPIT on the command
c line or enter them manually. The runlog should just be the runlog
c file name, it should not include the path (automatically looks in
c the runlogs directory)
      if (iargc() .lt. 1) then
          write(*,*) ' Enter name of runlog'
          read(*,*) runlog
      else
          call getarg(1, runlog)
      end if

      if (iargc() .lt. 2) then
          write(*,9000) 
 9000 format(' Std. TCCON processing: 3-hourly model and VMR? (y/n): ',
     &       $)
          read(*,*) modtype
      else
          call getarg(2, modtype)
      end if
c -------------------------------------------------------------

c -------- Parse the runlog/model user input ------------------
c Find the index of the last '.' in the runlog to get the viewing
c geometry.
      ldot=lnbc(runlog)-3
      if(runlog(ldot+1:ldot+1).eq.'g') ext='gnd'
      if(runlog(ldot+1:ldot+1).eq.'b') ext='bal'
c     write(*,*)ldot,ext

      if (trim(modtype) .eq. 'y' .or. trim(modtype) .eq. 'Y') then
        modtype = 'FPIT'
      else
        modtype = 'NCEP'
      end if
c -------------------------------------------------------------
      iwas=0
      open(lunw,file='extract_pth.out', status='unknown')
      write(lunw,*)4,10
      write(lunw,*)version
c JLL: Mysterious newlines appear in `runlog` when compiled under
c ifort. trim() seems to remove them (and so keep the header the
c correct number of lines)
      write(lunw,'(a)')' Runlog= '//trim(runlog)
      write(lunw,'(a)')' spectrum  year  doy   uthour   tout      tmod
     &       pout     pmod       hout     hmod'
      
      open(lunr_rlg,file=gggdir(:lrt)//'runlogs'//dl//ext//dl//runlog,
     & status='old')
      call read_runlog_header(lunr_rlg,data_fmt_read_rl,col_labels_rl)
      frdoy=-999
      latwas=-999
      lonwas=-999
      read_mod = .false. 
      modname_old = ''

      do i=1,999999
1        continue
         frwas=frdoy
         latwas=oblat
         lonwas=oblon
         call read_runlog_data_record(lunr_rlg,data_fmt_read_rl,
     &    col1,specname,iyr,idoy,zpdtim,oblat,
     &    oblon,obalt,asza,zenoff,azim,osds,
     &    opd,fovi,fovo,amal,ifirst,ilast,
     &    graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

          if(istat.ne.0) exit   !  EOF
          if(col1.eq.':') go to 1

c-----------------------------------------------------------
c  The next 20 lines of code determine the name of the new model file.
          if(oblat.gt.0.0) then
            ns='N'
          else
            ns='S'
          endif
          if(oblon.gt.0.0) then
            ew='E'
          else
            ew='W'
          endif

c -------- Read the model file -----------------------------
          call julian(iyr,1,idoy,jul)

c              ** NCEP **
          if (modtype .eq. 'NCEP') then
            frdoy=jul+zpdtim/24+oblon/360
            diff=abs(frdoy-frwas)+
     &      0.1*abs(oblat-latwas)+0.1*abs(oblon-lonwas)
            if( diff.gt.0.5 ) then
               jul=int(frdoy)
               call caldat(jul,iyyyy,imm,idd)
               write(modname,'(a5,i4.4,2i2.2,a1,i2.2,a1,a1,i3.3,a1,a4)')
     &          'NCEP_',iyr,imm,idd,'_',nint(abs(oblat)),ns,
     &          '_',nint(abs(oblon)),ew,'.mod'
c               write(*,*)frdoy,oblat,oblon,diff,specname(:14),
c     &                   modname(:28)
               read_mod = .true.
            end if ! if need to read new .mod file
          elseif (modtype .eq. 'FPIT') then
            call next_fpit_files(jul, zpdtim, oblat, oblon, modname,
     &         vmrname)
            if (modname .ne. modname_old) then
               read_mod = .true.
c               write(*,*)frdoy,oblat,oblon,specname(:14),' ',modname
            end if
            modname_old = modname
          end if ! which model type

         if (read_mod) then
            open(unit=lunm,
     &   file=gggdir(:lrt)//'models'//dl//ext//dl//modname,status='old')
            read(lunm,*)nlhead,ncol
            read(lunm,*)roc,ecc2,alat,gs,zz,pfact,ptrop
            do iline=1,nlhead-2
                read(lunm,*)
            enddo
            read(lunm,*)p0,t0,z0,m0,h0
            read(lunm,*)p1,t1,z1,m1,h1
            close(lunm)
         endif
c -----------------------------------------------------------

         tmod=t0 + (t1-t0)*(obalt-z0)/(z1-z0)
         hmod=h0 + (h1-h0)*(obalt-z0)/(z1-z0)
         pmod=p0*exp(-(obalt-z0)/8)
         scht=8.3144*t0/m0/gs/log1pxox(tmod/t0-1)
         pmod=p0*exp(-(obalt-z0)/scht)
         tout_k=tout+273.16
c         hout_vmr=(hout/100)*svp_wv_over_ice(tout_k)/pout
         hout_vmr=(hout/100)*wexler_svp(tout_k)/pout
         write(lunw,'(a47,2i4,7f10.4)') specname,iyr,idoy,zpdtim,
     &    tout+273.16,tmod,pout,pmod,hout_vmr,hmod
      end do
c      write(*,*)i-1
      close(lunr_rlg)
      stop
      end

      function wexler_svp(t)
c Wexler 1976 J Res Natl Bur Stand A Phys Chem (doi:
c 10.6028/jres.080A.071) gives an update to SVP over liquid water 
c coinciding with an update to the temperature scale that happened 
c in 1968. Their Eq. 15 is used for the Vaisala calibration in 
c Miloshevich et al. 2004 (doi:10.1175/1520-0426(2004)021%3C1305:DAVOAT%3E2.0.CO;2),
c Eq. 6.
c 
c Inputs:
c   t: REAL*4 - temperature in Kelvin
c
c Outputs:
c   REAL*4 - saturation vapor pressure of water over liquid water in
c   mbar

      implicit none

      real*4 t, tmp, wexler_svp

      tmp = 2.858487*log(t) 
     & + -2.9912729e3*t**-2 + -6.0170128e3*t**-1
     & + 1.887643854e1      + -2.8354721e-2*t
     & + 1.7838301e-5*t**2  + -8.4150417e-10*t**3
     & + 4.4412543e-13*t**4

c exp(tmp) will give the SVP in Pa, convert to hPa i.e. mbar
      wexler_svp = exp(tmp) * 1e-2 

      end function



      function svp_wv_over_ice(temp)
c Uses the Goff-Gratch equation to calculate the saturation vapor
c pressure of water vapor over ice at a user-specified temperature.
c   Input:  temp (K)
c   Output: svp_wv_over_ice (mbar)
      real*4 temp,t0,tr,yy,svp_wv_over_ice
      t0=273.16   ! triple point temperature (K)
      tr=t0/temp
      yy=-9.09718*(tr-1)-3.56654*alog10(tr)+0.876793*(1-1/tr)
      svp_wv_over_ice=6.1173*10**yy
      return
      end

