      subroutine read_mkiv_header(runlab,path,iend,ifirst,ilast,possp,
     & bytepw,apf,delwav,opd,fovi,snr,oblat,oblon,obalt,
     & pout,tout,hout,asza,iy,im,id,gmt,wavtkr,tins,pins,hins,lasf)
c
c  Subroutine to extract information from MkIV headers.
c
c INPUTS:
c     runlab  C*(*)  Name of spectrum
c       path  C*(*)  Path to spectrum
c       iend  I*4  Endianess of host computer

c OUTPUTS:
c     Everything else

      implicit none

      real*8
     & delwav,opd,fovi,snr,oblat,oblon,obalt,
     & pout,tout,hout,asza,azim,gmt,wavtkr,tins,pins,hins,lasf,
     & wsfact,dstep

      character
     & runlab*(*), path*(*), apf*2

      integer*4
     & ifirst,ilast,iend,iy,im,id,possp,bytepw,i4fftpow,fftsiz

c-------------------------------------------------------------------
c  The following are the MkIV spectral header parameters
c
      character blank*1640,comments*70,detctr*4,detsn*4,filter*4,
     & ippd*12,ippnam*6,preamp*4,rund*2,runinfo*18,runloc*4,sigchn*4,
     & strt*12
c
      integer*2 fftpow,ialias,idecim,irun,iset,iyr,phase,sfct
c
      integer*4 pinl,pspl,spsv,sspp,totp
c
      real*4 alat,alon,altd,ccor,fovr,gains(8),hint,hext,ippver,
     & offsts(8),pint,pinv,pott,pext,psym,rzero,sampl,soaz,soze,
     & spwn,stnr,tint,text,wdir,wspd,zmst,zpdv,zsym,zpdtim
c
      real*8 lsem,pspv,sins,sinv,ssps,sspv,zpdl

c--------------------------------------------------------------------
c  The following table documents the structure of the MkIV spectral header
c
c   Start End  Variable   Data
c   Byte  Byte  Name      Type    Description
c--------------------------------------------------------------------
c     1    18   runinfo   c*18   ! Interferogram header information
c    19    30   strt      c*12   ! start time YYMMDDHHMMSS
c    31    34   runloc    c*4    ! Run location (e.g. JPL, BAL, etc.)
c    35    36   iset      i*2    ! Set number (usually equal to day of year)
c    37    38   irun      i*2    ! Run number (typically 000 to 999)
c    39    40   rund      c*2    ! Run direction (F=forward, R=reverse)
c    41    42   ialias    i*2    ! Alias (always zero for MkIV spectra)
c    43    46   pint      r*4    ! Internal Pressure (mbar) inside instrument
c    47    50   tint      r*4    ! Internal Temperature (C) inside instrument
c    51    54   hint      r*4    ! Internal Humidity (%) inside instrument
c    55    62   lsem      r*8    ! laser semi-frequency (cm-1)
c    63    66   altd      r*4    ! Altitude of observation (km)
c    67    70   alat      r*4    ! Latitude of instrument (deg.)
c    71    74   alon      r*4    ! Longitude of instrument (deg.)
c    75    78   soaz      r*4    ! Solar Azimuth angle (deg.)
c    79    82   soze      r*4    ! Solar Zenith angle (deg.)
c    83    86   text      r*4    ! External Temperature outside (C)
c    87    90   pext      r*4    ! External Pressure outside (mbar)
c    91    94   hext      r*4    ! External Humidity Outside (%)
c    95    98   wspd      r*4    ! Wind Speed (Kts)
c    99   102   wdir      r*4    ! Wind Direction (deg.)
c   103   106   pinv      r*4    ! Peak interferogram value (DN)
c   107   110   pinl      i*4    ! Peak Interferogram location (index)
c   111   114   psym      r*4    ! Interferogram Symmetry (uncorrected)
c   115   122   zpdl      r*8    ! ZPD location (fractional index)
c   123   126   zpdv      r*4    ! Interferogram value at ZPD (DN)
c   127   130   zsym      r*4    ! Interferogram symmetry (after phase correction)
c   131   134   zmst      r*4    ! Time delay (s) between I'gram start and ZPD
c   135   138   sampl     r*4    ! Interferogram samples per fringe (=1 for MkIV)
c   139   142   totp      i*4    ! Total # of I'gram points on long side of ZPD
c   143   146   detctr    c*4    ! Detector serial number
c   147   150   filter    c*4    ! Filter serial number
c   151   154   preamp    c*4    ! Preamplifier serial number
c   155   158   sigchn    c*4    ! Signal chain serial number
c   159   162   fovr      r*4    ! Field of view diameter (mrad)
c   163   170   sinv      r*8    ! Sum of interferogram values (after phase correction)
c   171   178   sins      r*8    ! Sum of quares of I'gram values
c   179   190   ippd      c*12   ! Version number of IPP which created spectrum
c   191   192   fftpow    i*2    ! Size of FFT (log base 2 thereof)
c   193   194   idecim    i*2    ! Decimation factor
c   195   198   sspp      i*4    ! Starting spectral Point
c   199   202   spsv      i*4    ! Number of spectral points saved.
c   203   206   spwn      r*4    ! Mean spectral signal in chosen window
c   207   210   stnr      r*4    ! Spectral signal-to-noise ratio
c   211   218   sspv      r*8    ! Sum of spectral values
c   219   226   ssps      r*8    ! Sum of squares of spectral values.
c   227   228   phase     i*2    ! Subjective phase grade
c   229   232   pspl      i*4    ! Location (index) of peak spectral value
c   233   240   pspv      r*8    ! Peak spectral value
c   241   244   sfct      i*2    ! Scale factor (usually 15000) used to multiply spectral values
c   243   246   ippver    r*4    ! Version number of IPP program which created spectrum
c   247   316   comments  c*70   ! Comments
c   317   318   iyr       i*2    ! Year (full 4 digit value)
c   319   322   ccor      r*4    ! Cross correlation coefficient (dimensionless) 
c   323   326   pott      r*4    ! Potential temperature (K) outside instrument.
c   327   330   zpdtim    r*4    ! ZPD time (fractional hours)
c   331   334   detsn     c*4    ! detector serial number
c   335   340   ippnam    c*6    ! Version of ipp program that created spectrum
c   341   344   rzero     r*4    ! Non-linearity coefficient used to correct interferogram
c   345   376   gains(8)  r*4    ! Array of signal chain electrical gains.
c   377   408   offsts(8) r*4    ! Array of signal chain electrical offsets.
c   409  2048   blank     c*1640 ! Unused part of header.
c---------------------------------------------------------------------------
c  For MkIV, spectra can be big- or little-endian. Therefore, it is
c  necessary to look at ALIAS value in header to decide whether
c  byte reversal is required.
      iy=10*ichar(runlab(4:4))+ichar(runlab(5:5))-528
      if(iy.gt.80) then
        iy=1900+iy
      else
        iy=2000+iy
      endif
      im=1
      id=100*ichar(runlab(6:6))+10*ichar(runlab(7:7))+
     &   ichar(runlab(8:8))-5328
      open(19,file=path,form='unformatted',access='direct',
     &  status='old',recl=2048)
c      read(19,rec=1)i2hedr
      read(19,rec=1)runinfo,strt,runloc,iset,irun,rund,ialias,pint,
     & tint,hint,lsem,altd,alat,alon,soaz,soze,text,pext,hext,
     & wspd,wdir,pinv,pinl,psym,zpdl,zpdv,zsym,zmst,sampl,totp,detctr,
     & filter,preamp,sigchn,fovr,sinv,sins,ippd,fftpow,idecim,sspp,
     & spsv,spwn,stnr,sspv,ssps,phase,pspl,pspv,sfct,ippver,comments,
     & iyr,ccor,pott,zpdtim,detsn,ippnam,rzero,gains,offsts,blank
      close(19)

c IALIAS should be a small integer (e.g. 1,2,3,4) unless byte-reversed
      bytepw=2*iend
      if(ialias.ge.256) then  
        bytepw=-bytepw
        call rbyte(tint,4,1)
        call rbyte(pint,4,1)
        call rbyte(hint,4,1)
        call rbyte(totp,4,1)
        call rbyte(spsv,4,1)
        call rbyte(sspp,4,1)
        call rbyte(ialias,2,1)
        call rbyte(fftpow,2,1)
        call rbyte(idecim,2,1)
        call rbyte(lsem,8,1)
        call rbyte(soze,4,1)
        call rbyte(altd,4,1)
        call rbyte(alat,4,1)
        call rbyte(alon,4,1)
        call rbyte(soaz,4,1)
        call rbyte(soze,4,1)
        call rbyte(zpdtim,4,1)
        call rbyte(text,4,1)
        call rbyte(pext,4,1)
        call rbyte(hext,4,1)
        call rbyte(stnr,4,1)
      endif
c      write(*,*)'zpdtim=',strt,alat,alon,soaz,soze,text,pext

      tout=dble(text)
      pout=dble(pext)
      hout=dble(hext)
      tins=dble(tint)
      pins=dble(pint)
      hins=dble(hint)
      oblat=dble(alat)
      oblon=dble(alon)
      obalt=dble(altd)
      asza=dble(soze)
      azim=dble(soaz)
      snr =dble(stnr)
      gmt=dble(zpdtim)

      wavtkr=9900.d0
      possp=2048
      apf='BX'
      wsfact=1.00000D0
      if(lsem.le.0.0d0) lsem=7899.0015d0
      lasf=2.d0*lsem  ! the value in the header is the semi-frequency
      dstep=dble(idecim)
      i4fftpow=fftpow
      fftsiz=2**i4fftpow
      if(runlab(2:3).eq.'hg') then
        fovi=.0043
      elseif(runlab(2:3).eq.'in') then
        fovi=.0036
      else
        write(6,*)runlab,' ??'
        write(6,*)'This cannot be M4'
      endif
      ifirst=sspp+fftsiz*(ialias-1)-1
      ilast=ifirst+spsv-1
c  Next, account for the possibility of the spectral point spacing
c  being non-uniform, due to the presence of air inside the instrument.
c  This means that DELWAV is a function of wavenumber.
c  xnu is the spectral point address corresponding to a wavenumber nu
c  delwav is the average spectral point spacing between 0 and lasf/4 cm-1.
      delwav=lasf*wsfact/dstep/fftsiz/2
      opd=totp/lasf/wsfact
cc  Apply air-to-vacuum  & FOV corrections
c      if(lasf.lt.9999.) write(6,*)'LASF =',lasf
c      vbar=0.5*delwav*(ifirst+ilast)
c      delwav=delwav*riair(sngl(lasf),tint,pint,hint)/
c     $riair(vbar,tint,pint,hint)
      delwav=delwav*(1.D0+(fovi**2)/16)  ! FOV correction
      return
      end
