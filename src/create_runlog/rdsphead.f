      Subroutine rdsphead(spfmt,specname,path,ifirst,ilast,possp,bytepw,
     & apf,delwav,opd,fovi,snr,iy,im,id,gmt,lasf,wavtkr,lst,lse,lsu)

c
c  Subroutine to extract information from spectrum headers.
c
c INPUTS:
c       spfmt   Spectral header format type
c       specname  Name of spectrum file
c       path    Path name to spectrum

c OUTPUTS:
c       IFIRST  The index of the first spectral point on disk
c       ILAST   The index of the last spectral point on disk
c       POSSP   Number of bytes before IFIRST'th point (i.e. header length)
c       BYTEPW  Number of bytes per data word (=2 for AT,M4; =4 for OP,GR)
c       APF     Apodization that has been applied to spectrum
c       DELWAV  The wavenumber spacing of the raw spectral points
c       OPD     Optical path difference (cm).
c       FOVI    Angular diameter of internal field of view (radians)
c       SNR     Signal to noise ratio
c       OBLAT   Latitude of observer (deg)
c       OBLON   Longitude of observer (deg)
c       OBALT   Observation Altitude (km)
c       POUT    Observation Pressure (mbar)
c       TOUT    Observation Temperature (C)
c       HOUT    Observation Humidity (%)
c       ASZA    Astronomical Solar Zenith Angle (degrees)
c       IY      Year (4 digit)
c       IM      Month
c       ID      Day of the month (001 to 31)
c       GMT     UT time in decimal hours
c       WAVTKR  The frequency (cm-1) at which the suntracker operates.
c       TINS    Instrument Temperature (C) inside interferometer
c       PINS    Instrument Pressure (mbar) inside interferometer
c       HINS    Instrument Humidity (%)    inside interferometer
c
c Supports several different spectral formats including:
c     MkIV (1024 byte binary header)
c     ATMOS (no spectral header, read the runlog)
c     Bruker_Giessen (1280 byte binary header)
c     Bruker_OPUS ( variable length binary header)
c     FITS ( ASCII header, 80 bytes per record, variable # of records)
c     GRAMS (512 byte binary header)
c     JA (ASCII Jungfraujoch format)
c
      implicit none
      include "../ggg_const_params.f"
      include "params.f"

      character specname*(*),spfmt*2,path*(*),amon*3,apf*2,
     & hst*80,runtype*1,date*8,c16*2
      integer*4 dtype,hedlen,luns,
     & idwas,hh,mm,ss,ms,
     & ktype,nip,nnn,isr,dfr,pkl,prl
c      integer*2 mrs
      parameter (hedlen=512,luns=19)
c      integer*2 i2val(mrs)

      integer*4
     & jj,
     & fftsiz,
     $ nfwdscans,posnrun,i,ic,krec,iflag,gfw,gbw,
     & lnbc,npoints,irec,reclen,
     $ iend,nscan,iscan,nintp,ncenter,lst

      real*8 delwav,wsfact,dstep,nus,nue,nubar,
     $ apt,dur,lfl,hfl,oblat,oblon,obalt,
     & wspd,wdir,sia,sis,vdc,lse,lsu,
     $ startf,stopf,res,phr,vel,foc,d2r,tgmt,
     & zpdtim,t_s,dt
      parameter (d2r=dpi/180)

      real*8 fovi,gmt,gmtstart,gmtwas,lasf,
     & opd,scandelt,sampfreq,
     $ snr,zoff,resn,wavtkr,dopp,
     & pins,tins,hins,pout,tout,hout,
     & asza,ed,efl,aptdiam(0:15),sza1,sza2

      save aptdiam
      save hh,mm,ss,idwas,gmtwas,nip
      character*(hedlen) chhedr
      character*1024 logblock, bombin
      integer*2 i2hedr(hedlen/2), i2bombin(512)
      integer*4 i4hedr(hedlen/4), i4bombin(256)
      real*4    r4hedr(hedlen/4), r4bombin(256)
      real*8    r8hedr(hedlen/8), r8bombin(128)
      real*4 sigflo(8), sigfhi(8)
      integer*4 logoffset, charoff, ifil, ivac
c    & ios
      logical*1 pcat

c      equivalence (i2val,i4val,r8val,cval)  ! for OPUS
      equivalence (chhedr,i2hedr,i4hedr,r4hedr,r8hedr)
      equivalence (bombin,i2bombin,i4bombin,r4bombin,r8bombin)
c     DG 4/7/2002  - define actual bandpass for NDSC filters to avoid GFIT
c     fitting noise regions at edge of bandpass.
      DATA sigflo /3950.0, 2960.0, 2350.0, 2050.0, 1850.0,
     &            700.0,  980.0,  700.0/
      DATA sigfhi /4450.0, 3500.0, 3120.0, 2600.0, 2250.0,
     &            1350.0, 1350.0, 1025.0/


c      write(*,*)'rdsphead: ',path
      object=2  ! 1=moon; 2=sun
      call getendian(iend)  ! iend=+1 on Sun; iend=-1 on PC
c  Initialize variables to zero so that they will not be remembered from the
c  previous spectrum in the event that a spectrum header could not be read.
      ifirst=0
      ilast=0
      possp=0
      bytepw=0
      delwav=0
      opd=0.0d0
      fovi=0.0d0
      snr=0.0d0
      apf='XX'
      iy=0
      id=0
      gmt=0
      wavtkr=9900.0d0

c================================================================
c  For MkIV, spectra can be big- or little-endian. Therefore, it is
c  necessary to look at ALIAS value in header to decide whether
c  byte reversal is required.
      if(spfmt.eq.'m4') then
      call read_mkiv_header(specname,path,iend,ifirst,ilast,possp,
     & bytepw,apf,delwav,opd,fovi,snr,oblat,oblon,obalt,
     & pout,tout,hout,asza,iy,im,id,gmt,wavtkr,tins,pins,hins,lasf)
c=====================================================================
c  GRAMS format for Wollongong solar BOMEM spectra collected by SUNRUN.
c  DG Jan03
      elseif(spfmt.eq.'D8') then
        call read_DA8_header(specname,path,iend,ifirst,ilast,possp,
     &   bytepw,apf,delwav,opd,fovi,snr,pout,tout,hout,asza,iy,im,id,
     &   gmt,wavtkr,tins,pins,hins,lasf)
c=====================================================================
      elseif(spfmt.eq.'at') then
        iy=1900+10*ichar(specname(2:2))+ichar(specname(3:3))-528
        im=10*ichar(specname(4:4))+ichar(specname(5:5))-528
        id=10*ichar(specname(6:6))+ichar(specname(7:7))-528

c  ATMOS spectral header information is stored in a RUNLOG
c  Get the RUNLOG record and extract pertinent info
       open(19,file='/runlogs/atmos_runlog',access='direct',recl=hedlen)
        read(19,rec=posnrun(19,specname,16,79109)) i2hedr
        close(19)
        dopp=r4hedr(104)
        wavtkr=9900.
        bytepw=2
        possp=0
        apf='BX'
        fftsiz=i4hedr(57)
        lasf=2.d0*r8hedr(61)
        wsfact=r8hedr(62)
        dstep=dble(r4hedr(71))
c        fovo=dble(r4hedr(101))/1000.  ! convert milli-radians to radians
c        fovi=fovo*dble(r4hedr(100))
        fovi=dble(r4hedr(100))*dble(r4hedr(101))/1000. 
        pins=dble(r4hedr(74))
        tins=dble(r4hedr(73))
        hins=dble(r4hedr(72))
c        tanght=dble(r4hedr(83))
c
c        oblat=r4hedr(84)
c        oblon=r4hedr(85)
c        obalt=r4hedr(81)
c        asza=r4hedr(82)
c        tout=r4hedr(97)
c        pout=r4hedr(98)
c        hout=r4hedr(99)
c        oblat=34.3819
c        oblon=-117.6776
c        obalt=2.258
c
        snr=dfloat(i4hedr(64))/r4hedr(80)
        gmt=dble(r4hedr(77)/3600.d0)
        ifirst=i4hedr(59)+fftsiz*(i4hedr(51)-1)-1
        ilast=ifirst+i4hedr(60)-1
        delwav=lasf*wsfact*(1.D0+dopp)/dstep/fftsiz/2
        opd=dfloat(i4hedr(55))/lasf/wsfact

c==================================================================
      elseif (spfmt.eq.'fi') then    ! FITS format
c        oblat=31.958d0   ! Kitt Peak
c        oblon=-111.595d0 ! Kitt Peak
c        obalt=2.092d0     ! Kitt Peak
        hh=0  ! avoids NaN if time not present in header
        mm=0
        ss=0
        nscan=2
        fovi=0.002d0
        bytepw=4
        apf='BX'
        resn=0.02  ! In case OPD value is not found in spectrum header
        snr=500.0d0
        reclen=80
        irec=0
        iflag=0
        hst=' '
        open(luns,file=path,form='unformatted',status='old',
     &  access='direct',recl=reclen)
        do while (iflag .eq. 0 )
           do krec=1,36
              irec=irec+1
              read(luns,rec=irec) hst       ! header string
              if    (hst(:9).eq.'WSTART  =') then
                 read(hst(10:),*)nus
              elseif(hst(:9).eq.'WSTOP   =') then
                 read(hst(10:),*)nue
              elseif(hst(:9).eq.'NPO     =') then
                 read(hst(10:),*)npoints
              elseif(hst(:9).eq.'DELW    =') then
                 read(hst(10:),*)delwav
              elseif(hst(:9).eq.'RESOLUTN=') then
                 read(hst(10:),*)resn
              elseif(hst(:9).eq.'ID      =') then
c                 write(*,'(a)') path(:lnbc(path)),hst(10:)
              elseif(hst(:9).eq.'DAY     =') then
                 read(hst(10:),*) date
                 read(date,'(i2,1x,i2,1x,i2)')im,id,iy 
              elseif(hst(:9).eq.'ZENSTRT =') then
                 read(hst(10:),*)sza1 
              elseif(hst(:9).eq.'ZENSTOP =') then
                 read(hst(10:),*)sza2 
              elseif(hst(:9).eq.'NINT    =') then
                 read(hst(10:),*)nintp
              elseif(hst(:9).eq.'NCENTER =') then
                 read(hst(10:),*)ncenter
              elseif(hst(:9).eq.'SAMPFREQ=') then
                 read(hst(10:),*)sampfreq
              elseif(hst(:9).eq.'TIMESTR =') then
                 read(hst(10:),*)date
                 ic=index(date,':')
                 if(ic.eq.2) read(date,'(i1,1x,i2,1x,i2)')hh,mm,ss
                 if(ic.eq.3) read(date,'(i2,1x,i2,1x,i2)')hh,mm,ss
              elseif(hst(:9).eq.'NSCAN   =') then
                 read(hst(10:),*)nscan
              elseif(hst(:9).eq.'ELTIM   =') then
                 read(hst(10:),*)dur
              elseif(hst(:9).eq.'PATM    =') then
                 read(hst(10:),*)pout
              elseif(hst(:9).eq.'TATM    =') then
                 read(hst(10:),*)tout
              elseif(hst(:9).eq.'REFWAVNO=') then
                 read(hst(10:),*)lasf
              elseif(hst(:9).eq.'END      ') then
                 iflag=1
              endif
           end do  ! krec=1,36
        end do  ! while iflag eq 0
        possp=80*irec
        if(iy.lt.50) then
           iy=iy+2000
        else
           iy=iy+1900     
        endif
        opd=0.5/resn
        pout=1013.25*pout/760
c        write(*,*)nus,nue,tout,pout
c        pins=0.
c        tins=25.
c        hins=10.
        ifirst=nus/delwav+0.5        ! subscript of first point in file
        if(abs(nus/delwav-ifirst) .gt. .04) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nus,nus/delwav
        endif
        ilast=nue/delwav+0.5        ! subscript of first point in file
        if(abs(nue/delwav-ilast) .gt. .04) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nue,nue/delwav
        endif
c
c  Compute mean ASZA as the average of the ZPD SZAs of the individual interferograms.
        tgmt=0.0d0
        zpdtim=hh+(mm+(ss+float(ncenter)/sampfreq)/60.)/60. 
        scandelt=(dur-dfloat(nintp)/sampfreq)/(nscan-1)
        do iscan=1,nscan
          tgmt=tgmt+zpdtim
          zpdtim=zpdtim+scandelt/3600
        end do
        gmt=tgmt/nscan ! airmass-weighted average GMT
c

c  Read the lab spectra taken on 10/29/2003 at Kitt Peak, using 100m 
c  path length cell at 3 Bar pressure, by Brown and Miller.
      
      elseif (spfmt.eq.'kp') then  !yzh 20031104
c
c  Read Header Block.
        oblat=31.958d0
        oblon=-111.595d0
        obalt=2.092d0
        bytepw=4
        apf='BX'
        snr=500.0d0
        irec=0
        reclen=100
        iflag=0
        hst=' '
        path=path(:lnbc(path)-3)//'hdr'   ! These spectra have seperate headers
        open(luns,file=path,status='old') 
        write(*,*) '.hdr file opened!'  !yzh temp
        do while (iflag .eq. 0 )
           do krec=1,99
              irec=irec+1
              read(luns,'(a)') hst       ! header string
              if    (hst(:9).eq.'id      =') then
c                 write(*,'(a)') hst(10:)
              elseif(hst(:9).eq.'day     =') then
                 read(hst(10:),*) date
                 read(date,'(i2,1x,i2,1x,i2)')im,id,iy
                 iy=iy+2000
              elseif(hst(:9).eq.'wstart  =') then
                 read(hst(10:),'(4x,f12.7)')nus
              elseif(hst(:9).eq.'wstop   =') then
                 read(hst(10:),'(4x,f12.7)')nue
              elseif(hst(:9).eq.'npo     =') then
                 read(hst(10:),'(4x,i6)')npoints
              elseif(hst(:9).eq.'delw    =') then
                 read(hst(10:),'(3x,f13.12)')delwav
              elseif(hst(:9).eq.'resolutn=') then
                 read(hst(10:),'(3x,f13.12)')resn
              elseif(hst(:9).eq.'nscan   =') then
                 read(hst(10:),'(2x,i4)')nscan
              elseif(hst(:9).eq.'timestr =') then
                 read(hst(10:),*)date
                 ic=index(date,':')
                 if(ic.eq.2) read(date,'(i1,1x,i2,1x,i2)')hh,mm,ss
                 if(ic.eq.3) read(date,'(i2,1x,i2,1x,i2)')hh,mm,ss
              elseif(hst(:9).eq.'eltim   =') then
                 read(hst(10:),'(2x,f8.3)')dur
              elseif(hst(:9).eq.'apin    =') then
                 read(hst(10:),'(2x,f4.2)')fovi
                 fovi=fovi/1000.0
              elseif(hst(:9).eq.'refwavno=') then
                 read(hst(10:),'(3x,f11.6)')lasf
              elseif(hst(:9).eq.'sampfreq=') then
                 read(hst(10:),'(3x,f11.6)')sampfreq
              elseif(hst(:9).eq.'nint    =') then
                 read(hst(10:),'(3x,i8)')nintp
              elseif(hst(:9).eq.'ncenter =') then
                 read(hst(10:),'(3x,i8)')ncenter
              elseif(hst(:9).eq.'END      ') then
                 iflag=1
              endif
           end do  ! krec=1,100
        end do  ! while iflag eq 0
        write(*,*) 'Finish reading .hdr file'  !yzh temp
        possp=0
        opd=0.5/resn
c        pout=800.0
c        pins=0.
c        tins=25.
c        hins=10.
        ifirst=aint(nus/delwav+0.5)        ! subscript of first point in file
        if(abs(nus/delwav-ifirst) .gt. .04) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nus,nus/delwav
        endif
        ilast=aint(nue/delwav+0.5)        ! subscript of first point in file
        if(abs(nue/delwav-ilast) .gt. .04) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nue,nue/delwav
        endif
c
        gmt=7.d0+hh+(mm+(ss+dur/2.)/60.d0)/60.d0     !yzh 20031104
c
c============================================================
c      elseif (spfmt.eq.'kp') then
c        open(19,file='/runlogs/kitt.log',access='direct',recl=400)
c  For some unknown reason the first four bytes of the KP header are nulls.
c        filter=char(0)//char(0)//char(0)//char(0)
c        read(19,rec=posnrun(19,filter//specname,16,205))specnamel,
c     &  comments,comments,
c     &  comments,strt,fftd,zstrt,zstop,strt,resn,detector,detector,pout,
c     &  tout,hout,pins,tins,snr,cspare1,detector,nus,nue,delwav
c        close(19)
c
c        nus=r8hedr(46)
c        nue=r8hedr(47)
c        delwav=r8hedr(48)
c        write(6,*)nus,nue,delwav
c        fovi=0.0
c        obalt=2.092
c        oblat=33.0
c        oblon=-100.
c        bytepw=4
c        snr=1200.
c        ifirst=nus/delwav+0.5        ! subscript of first point in file
c        ilast=nue/delwav+0.5        ! subscript of first point in file
c=======================================================================
      elseif (spfmt.eq.'ja') then  !jj ascii format with block markers removed
        open(19,file=path,status='old')
        read(19,'(///,8x,i2,1x,a3,i5)')id,amon,iy
        read(19,'(//,17x,f7.3)')       gmt
        read(19,'(23x,f9.4)')          nus
        read(19,'(23x,f9.4)')          nue
        read(19,'(14x,f5.1)')          resn
        read(19,'(22x,f5.1,29x,f3.0)') ed,efl
        read(19,'(24x,f7.3)')          asza
        read(19,'(24x,f6.0,//////)')   snr
        read(19,*)jj,nus,delwav
        close(19)
        asza=-asza       ! zenith angles are refracted
        gmt=gmt-1.d0     ! convert local time to UT
        if(amon.eq.'JAN')im=1
        if(amon.eq.'FEB')im=2
        if(amon.eq.'MAR')im=3
        if(amon.eq.'APR')im=4
        if(amon.eq.'MAY')im=5
        if(amon.eq.'JUN')im=6
        if(amon.eq.'JUL')im=7
        if(amon.eq.'AUG')im=8
        if(amon.eq.'SEP')im=9
        if(amon.eq.'OCT')im=10
        if(amon.eq.'NOV')im=11
        if(amon.eq.'DEC')im=12

        possp=1600  !  bytes in ascii header (20 lines x 80 bytes)
        fovi=0.1*ed/efl    ! internal FOV radius in milli-radians
        resn=0.001*resn    ! convert mK to cm-1
        opd=0.9/resn
c        obalt=3.58d0
c        pout=660.
c        tout=-5.
c        oblat=46.55
c        oblon=7.98
        bytepw=5          !  denotes an ascii file
        apf='BX'
        ifirst=nus/delwav+0.5        ! subscript of first point in file
        if(abs(nus/delwav-ifirst) .gt. .04) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nus,nus/delwav
        endif
        ilast=nue/delwav+0.5        ! subscript of first point in file
        if(abs(nue/delwav-ilast) .gt. .04) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nue,nue/delwav
        endif
c===================================================================
c  Bruker Giessen format spectra are always little-endian and therefore
c  BYTEPW is -ve and their headers require byte-reversing on the Sun (iend=+1)
      elseif(spfmt.eq.'br') then
        possp=1280
        bytepw=-4
        apf='BX'
        open(19,file=PATH,form='unformatted',status='old',
     &  access='direct',recl=hedlen) 
        read(19,rec=1)i4hedr
        if(iend.gt.0) call rbyte(i4hedr(1),4,hedlen/4)

c*JFB* begin change 23-Aug-94
c*     Starting with ACQUISIT V.2.6.0, released 13-Nov-92, the time is
c*     derived from the DOS time instead of the file creation time. This
c*     is more precise (1/18th of a second), but has the side effect of
c*     giving the year in full (century included). This change went unnoticed
c*     until now! It is actually a desirable feature, as the year 2000 should
c*     be correctly handled, if DOS does not crash on that date.
        iy=i4hedr(74)/4096
        if(iy.lt.1900) then
          iy=iy+1900     ! fix runs recorded before ACQUISIT V.2.6.0
        endif
c*JFB* end change 23-Aug-94

        im=mod(i4hedr(74)/64,64)
        id=mod(i4hedr(74),64)

        gmt= dble(i4hedr(75)/262144                 ! hours
     &  + float(mod(i4hedr(75)/4096,64))/60    ! minutes
     &  + float(mod(i4hedr(75),4096))/36000)    ! tenths of seconds

c*JFB* begin change 23-Aug-94
        gmt=gmt-1.d0  !  For historical reasons, ISSJ runs at UT+1    (!!!)
c*JFB* end change 23-Aug-94
        wavtkr=9900.
c        pins=1000.0
c        tins=20.0
c        hins=20.0
c        oblat=46.55
c        oblon=7.98
c        tout=-5.0
c        pout=660.0
c        hout=50.0
c        obalt=3.58d0
        snr=1000.0d0
        nus=dfloat(i4hedr(70))+dfloat(i4hedr(71))/2**24
        nue=dfloat(i4hedr(72))+dfloat(i4hedr(73))/2**24

c*JFB* begin change 23-Aug-94
c*     Table of field stop diameters, in millimeters, for ISSJ Bruker
        data aptdiam /0.6, 0.9, 1.1, 1.2, 1.45, 1.55, 1.75, 2.0,
     &                2.5, 3.2, 4.05, 5.05, 6.4, 8.0, 10.0, 15.0/
c*     Use header entry to find internal field-of-view.
c*     The table 'aptdiam' should really be read from 'instrum.des'
c*     but this will do for now...
        foc=418.  ! mm
        fovi=aptdiam(i4hedr(90))/foc
c*JFB* end change 23-Aug-94
        delwav=(nue-nus)/i4hedr(66)
        ifirst=int(nus/delwav+0.5)
        if(abs(nus/delwav-ifirst) .gt. .1) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nus,i4hedr(66),nus/delwav,ifirst
        endif
        ilast=int(nue/delwav+0.5)-1
        if(abs(nue/delwav-1-ilast).gt..1) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nue,i4hedr(66),nue/delwav,ilast
        endif
        read(19,rec=2)i4hedr
        if(iend.gt.0) call rbyte(i4hedr(5),4,2)
        lasf=dfloat(i4hedr(5))+dfloat(i4hedr(6))/2**24
        read(19,rec=3)i4hedr
        if(iend.gt.0) call rbyte(i4hedr(7),4,2)
        resn=dfloat(i4hedr(7))+dfloat(i4hedr(8))/2**24
        opd=0.5/resn
        close(19)
c===================================================================
c  FTUVS, using Bruker Giessen format: always little-endian and therefore
c  BYTEPW is -ve and their headers require byte-reversing on the Sun (iend=+1)
      elseif(spfmt.eq.'uv') then
        possp=1280
        bytepw=-4
        apf='BX'
        open(19,file=PATH,form='unformatted',status='old',
     &  access='direct',recl=hedlen) 
        read(19,rec=1)i4hedr
        if(iend.gt.0) call rbyte(i4hedr(1),4,hedlen/4)

        iy=i4hedr(74)/4096
        if(iy.lt.1900) then
          iy=iy+1900     ! fix runs recorded before ACQ V.2.6.0
        endif

c       iy=2001          ! abominable kludge

        im=mod(i4hedr(74)/64,64)
        id=mod(i4hedr(74),64)

        gmt= dble(i4hedr(75)/262144            ! hours
     &  + float(mod(i4hedr(75)/4096,64))/60    ! minutes
     &  + float(mod(i4hedr(75),4096))/36000)   ! tenths of seconds

c       gmt=gmt + 6.95  !  KLUDGE: -1h3' to get to local time, +8h to GMT

        wavtkr=15380.
c        pins=769.4
c        tins=20.3
c        hins=15.0
c        oblat=34.382
c        oblon=-117.68
c        obalt=2.286d0
c        tout=12.1
c        pout=769.4
c        hout=41.0
        snr=200.0d0
        nus=dfloat(i4hedr(70))+dfloat(i4hedr(71))/2**24
        nue=dfloat(i4hedr(72))+dfloat(i4hedr(73))/2**24
        nubar=(nus+nue)/2
        ktype=index(specname,'.dat')-1
        if(ktype.le.0) then
           write(*,*) 'Warning: Non-conforming FTUVS spectrum name'
           write(*,*) 'No time correction for run ', specname
        else
           runtype=specname(ktype:ktype)
           if(     runtype .eq. 'r' ) then  ! Reverse Scan
               gmt=gmt+0.00224d0*dble(nubar-1750)/3600
           elseif( runtype .eq. 'f' ) then  ! Forward scan
               gmt=gmt
           else                             ! F/R Averaged
               gmt=gmt+0.00112d0*dble(nubar-1750)/3600
           endif
        endif

        foc=203.3        ! Collimator focal length, in mm
        fovi=1.12/foc    ! Field stop of 1.12 mm in diameter
        delwav=(nue-nus)/i4hedr(66)
        ifirst=int(nus/delwav+0.5)
        if(abs(nus/delwav-ifirst) .gt. .1) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nus,i4hedr(66),nus/delwav,ifirst
        endif
        ilast=int(nue/delwav+0.5)-1
        if(abs(nue/delwav-1-ilast).gt..1) then
          write(6,*)'RDSPHEAD Warning: Non-integral spectral grid'
          write(6,*)nue,i4hedr(66),nue/delwav,ilast
        endif
        read(19,rec=2)i4hedr
        if(iend.gt.0) call rbyte(i4hedr(5),4,2)
        lasf=dfloat(i4hedr(5))+dfloat(i4hedr(6))/2**24
        read(19,rec=3)i4hedr
        if(iend.gt.0) call rbyte(i4hedr(7),4,2)
        resn=dfloat(i4hedr(7))+dfloat(i4hedr(8))/2**24
        opd=0.5/resn
        close(19)
c        write(*,*)nus,nue,delwav,resn,opd
c  Correct for ACQ laser wavenumber
        delwav=delwav*(15798.00364D0/lasf)
c  Correct for laser off-axis (to avoid laser back-feed)
        delwav=delwav*(1.D0-(0.00283**2)/2)
c===================================================================
c  OPUS spectra are always little-endian. Therefore, test IEND to see
c  whether byte-reversal is required.
      elseif(spfmt.eq.'op') then
        opd=0.0
        foc=1.E+18        ! default value (avoids zerodivide later)
        bytepw=-4
        zoff=0.
        snr=1000.0d0
c
        dtype=1031 ! data_type (1031=spectrum)
        c16=specname(16:17)
        if(c16.eq.'A.' .or. c16.eq.'B.' .or. c16.eq.'X.') dtype=2055 ! interferogram

       call read_opus_header(path,iend,dtype,npoints,startf,stopf,iy,im,
     &  id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lasf,foc,nip,dfr,
     &  pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,tins,pins,hins,
     &  tout,pout,hout,wspd,wdir,sia,sis,vdc,lst,lse,lsu)

c        call rdopushead(path,iend,dtype,npoints,startf,stopf,iy,im,id,
c     &  hh,mm,ss,ms,apt,nscans,dur,vel,apf,phr,res,lasf,foc,nip,dfr,
c     &  pkl,prl,gfw,gbw,lfl,hfl,possp)
c        write(*,'(a,7i6)')'iy,im,id,hh,mm,ss,ms=',iy,im,id,hh,mm,ss,ms
c
        delwav=(stopf-startf)/(npoints-1)
        ifirst=nint(startf/delwav)
        ilast=nint(stopf/delwav)
        delwav=15798.0138d0*delwav/lasf  ! Some OPUS Bruker data has wrong laser freq (e.g. 15798.10)
c        write(*,'(a,2f16.9,i6,f16.11)')'startf,stopf,npoints= ',
c     &  startf,stopf,npoints,delwav
        opd=0.9/res
        fovi=apt/foc
c        write(*,*)path,apt,foc,fovi
        delwav=delwav*(1.D0+(fovi**2)/16)  ! FOV correction
c        write(*,'(2a,2f15.8,i8,f16.12,f16.6,f8.5)')
c     &  'start,stop,npts,del,lasf= ',
c     &  specname,startf,stopf,npoints,delwav,lasf,fovi
        if(dtype.eq.2055) then   ! Interferogram
           delwav=1.d0
           ifirst=1
           ilast=nip
        endif
        gmtstart=dble(hh+(mm+(ss+ms/1000.)/60.)/60.)  ! Scan start time
c        write(*,*)'  gmtstart  Delta_t(s) =',specname,gmtstart,
c     &  nint(((id-idwas)*24.0d0+gmtstart-gmtwas)*3600)
        idwas=id
        gmtwas=gmtstart
c  Some old spectra have GFW=0 and GBW=0, so can't tell whether
c  it was a FWD or REV scan. So adopt the following kludge.
        if(gfw.eq.0 .and. gbw.eq.0) then
          if(pkl.gt.prl) gfw=1
          if(pkl.le.prl) gbw=1
        endif

        if(hfl.eq.0 .and. lfl.eq.0) hfl=lasf
        if(nip.le.0) then
           if(pkl.le.0 .or. prl.le.0) then
              nip=2*0.9*(hfl-lfl)*(1/res+1/phr)/dfr
           else
              nip=pkl+prl  ! Primary Sites (slice-ipp)
c              write(*,*)'pkl,prl,nip=',pkl,prl,nip
              if(gfw.le.0) pkl=0
              if(gbw.le.0) prl=0
           endif
        elseif(gfw.eq.gbw) then
           nip=nip/2  ! For FWD+REV, OPUS defines NIP=2*(PKL+PRL)
        endif
        nnn=nip
        if(gfw.eq.gbw) nnn=nip*2

        isr=nint(2000*vel*(hfl-lfl)/lasf)  ! Igram Sampling Rate (Hz)
        if(dfr.gt.1) isr=isr/dfr           ! Digital Filter Reduction
        t_s=dfloat(nip)/isr                ! Scan Time
c
c Reverse REV scans that OPUS has flipped to make appear as FWD scan
c        write(*,*) gfw, gbw, pkl, prl, nip, isr, t_s, dur
        if(gbw.gt.0 .and. prl.lt.nip/2) prl=nip-prl

        dt=dur/2+t_s*(dfloat(pkl+prl)/nnn-0.5d0)
c        write(*,*) gfw, gbw, pkl, prl, nip, isr, t_s, dur, dt
        gmt=gmtstart+dt/3600.
c        write(*,*)' gmtstart: ',dt,dur,pkl,prl,nnn

c===================================================================
c      elseif(spfmt.eq.'kp') then
c        openr(luns,file=path(lnbc(path)-3)//'hdr',
c     & form='unformatted',status='old',access=
c     & 'direct',recl=reclen)
c        read(luns,*)
c        read(luns,'(a9,i10)') inftype
c        read(luns,'(a9,3(1x,i2))') im,id,iy
c	iy=iy+2000
c        read(luns,'(a9,i18)') scanser
c	read(luns,'(a9,f17.7)') startf
c        read(luns,'(a9,f17.7)') stopf
c	read(luns,'(a9,i10)') npoints
c        read(luns,*)  !data type
c	read(luns,'(a9,f17.13)') delwav
c        read(luns,'(a9,f17.13)') res
c	read(luns,'(a9,2x,i2)') fboffs
c	read(luns,*)
c	read(luns,*)
c	read(luns,'(a9,5x,e12.4)' wavcorr
c	read(luns,'(a9,5x,e12.4)' rdsclfct
c	read(luns,'(a9,5x,e12.4)' noiselev
c	read(luns,*)
c	read(luns,*)
c	read(luns,*)
      else
        stop 'RDSPHEAD: unknown spectral format'
      endif
c============================================================
      return
      end

      integer function posnrun(unit,specname,nchr,nrun)
c  Positions in the (already opened) runlog to spectrum RUNLAB
      integer unit,nchr,jhi,jlo,jj,nrun
      character specname*(*),specnamel*20
c
c  First try a bisecting search (see function LOCATE in Numerical Recipes)
c  This only takes about log2(NRUN) iterations if the list is ordered.
      jlo=0
      jhi=nrun+1
      do while(jhi-jlo.gt.1)
         jj=(jlo+jhi)/2                ! Bisect remaining range of positions
         read(unit,rec=jj) specnamel(:nchr)
         if(specnamel.le.specname) then
            jlo=jj
         else
            jhi=jj
         endif
      end do
c
c  JLO points to where the run should be in the list.
c  If jlo=0  or  jlo=nrun+1 then it definately did not find the run
      if(jlo.gt.0 .and. jlo.le.nrun) then
         posnrun=jlo
         read(unit,rec=jlo) specnamel(:nchr)
         if(specnamel.eq.specname) return   ! success
      endif
c
c  If the ordered search failed then it can mean that either the run
c  is not in the runlog, or the runlog is not sorted properly.
c  So try a complete systematic search (very slow)
      do posnrun=1,nrun
        read(unit,rec=posnrun) specnamel(:nchr)
        if(specnamel.eq.specname) return   ! success
      end do
      write(6,*)'POSNRUN:', specname,' not present in runlog '
      write(6,*)'May need to increase parameter NRUN currently = ',nrun
      stop
      end
