      subroutine read_DA8_header
     &  (specname,path,iend,ifirst,ilast,possp,bytepw,apf,delwav,opd,
     &   fovi,snr,pout,tout,hout,asza,iy,im,id,gmt,wavtkr,tins,pins,
     &   hins,lasf) 

c     Read DA8-GRAMS format for Wollongong solar BOMEM spectra collected by SUNRUN.
c     DG Jan03
c     Revised Jun2007 as subroutine for new sunrun/runlog split

      implicit none

      character   specname*(*),path*(*),apf*2
      integer*4   i,hedlen,hh,mm,ss,
     &            ifirst,ilast,possp,bytepw,
     &            iy,im,id,jd,
     &            nfwdscans,npoints,iend, ispeed
      parameter   (hedlen=512)

      real*8      delwav,nus,nue,fovi,asza,gmt,lasf,
     &            opd,snr,resn,wavtkr,pins,tins,hins,pout,tout,hout

      character*(hedlen) chhedr
      integer*2   i2hedr(hedlen/2), i2bombin(512)
      integer*4   i4hedr(hedlen/4), i4bombin(256)
      real*4      r4hedr(hedlen/4), r4bombin(256)
      real*8      r8hedr(hedlen/8), r8bombin(128)
      real*4      sigflo(8), sigfhi(8), speed(16)
      character*1024 logblock, bombin
      integer*4   logoffset, charoff, ifil, ivac, iostat
      logical*1   pcat


      equivalence (chhedr,i2hedr,i4hedr,r4hedr,r8hedr)
      equivalence (bombin,i2bombin,i4bombin,r4bombin,r8bombin)

c     DG 4/7/2002  - define actual bandpass for NDSC filters to avoid GFIT
c     fitting noise regions at edge of bandpass.
      DATA sigflo /3950.0, 2800.0, 2350.0, 2050.0, 1850.0,
     &            700.0,  980.0,  700.0/
      DATA sigfhi /4450.0, 3500.0, 3120.0, 2600.0, 2250.0,
     &            1350.0, 1350.0, 1025.0/
      DATA speed /0.01, 0.02, 0.03, 0.05, 0.07, 0.10,
     &            0.15, 0.20, 0.30, 0.50, 0.70, 1.00,
     &            1.50, 2.00, 3.00, 4.6/

      tins=-999.0d0
      tout=-999.0d0
      hins=-999.0d0
      hout=-999.0d0
      iend=-1
c      idum=iend    ! Avoid compiler warning (unused variable)
c      rdum=hout    ! Avoid compiler warning (unused variable)
c      rdum=tout    ! Avoid compiler warning (unused variable)
c      rdum=hins    ! Avoid compiler warning (unused variable)
c      rdum=tins    ! Avoid compiler warning (unused variable)
c      idum=mfilepath ! Avoid compiler warning (unused parameter)
c      idum=mauxcol ! Avoid compiler warning (unused parameter)
c      idum=mcolvav ! Avoid compiler warning (unused parameter)
c      idum=mgas    ! Avoid compiler warning (unused parameter)
c      idum=mlev    ! Avoid compiler warning (unused parameter)
c      idum=mrow_qc ! Avoid compiler warning (unused parameter)
c      idum=mspeci  ! Avoid compiler warning (unused parameter)
c      idum=mvmode  ! Avoid compiler warning (unused parameter)
c      idum=ncell   ! Avoid compiler warning (unused parameter)

      open(19,file=path,form='unformatted',access='direct',
     &status='old',recl=hedlen)   !!!!! needs compiler option /assume:byterecl
c     Read Grams header 512 bytes
      read(19,rec=1)i4hedr
      nus=r8hedr(2)           !start freq.
      nue=r8hedr(3)           !end freq.
      npoints=i4hedr(2)
c     Read in Bomem binary block 1024 bytes at offset flogoff+64
      logoffset=i4hedr(63)+64                               !i4hedr(63)=flogoff
      read(19,rec=int(logoffset/hedlen)+1,iostat=iostat)chhedr
      if(iostat.lt.0)goto 100
      charoff=logoffset-int(logoffset/hedlen)*hedlen
      bombin(1:hedlen-charoff)=chhedr(charoff+1:hedlen)
      read(19,rec=int(logoffset/hedlen)+2,iostat=iostat)chhedr
      if(iostat.lt.0)goto 100
      bombin(hedlen-charoff+1:2*hedlen-charoff)=chhedr
      read(19,rec=int(logoffset/hedlen)+3,iostat=iostat)chhedr
      bombin(2*hedlen-charoff+1:1024)=chhedr
      if(iostat.lt.0)goto 100
c     Read in Log text block 1024 bytes at offset flogoff+64+1024
      logoffset=logoffset+1024
      read(19,rec=int(logoffset/hedlen)+1,iostat=iostat)chhedr
      if(iostat.lt.0)goto 100
      charoff=logoffset-int(logoffset/hedlen)*hedlen
      logblock(1:hedlen-charoff)=chhedr(charoff+1:hedlen)
      read(19,rec=int(logoffset/hedlen)+2,iostat=iostat)chhedr
      if(iostat.lt.0)goto 100
      logblock(hedlen-charoff+1:2*hedlen-charoff)=chhedr
      read(19,rec=int(logoffset/hedlen)+3,iostat=iostat)chhedr
      logblock(2*hedlen-charoff+1:1024)=chhedr
      if(iostat.lt.0)goto 100
100   close(19)

c     Get filter and date from file name
      read(specname(1:8),'(i1,i2,z1,i2)')ifil,iy,im,id
      if(iy.gt.90)then
         iy=iy+1900
      else
         iy=iy+2000
      endif

c     Get stuff from Bomem binary block
      nfwdscans=i4bombin(24)
      ivac=i2bombin(274)
      ispeed=i2bombin(262)
      pout=i2bombin(319)/3276.8*260.0/5+800.0-0.1
      if (pout.le.950.0.or.pout.gt.1040.0)pout=1013.0
c     Note ivac is incorrect in PCAT1.2c - always=0

c     Get stuff from text log block (same for both PCDA and PCAT)
      i=index(logblock,'Start time = ')+13
      if(i.eq.13) then
         write(*,*) specname
         stop 'Cannot find "Start time"'
      endif
      if(ichar(logblock(i+7:i+7)).eq.13)then         !13=CR, single digit seconds
         read(logblock(i:i+6),'(i2,x,i2,x,i1)')hh,mm,ss
      else
         read(logblock(i:i+7),'(i2,x,i2,x,i2)')hh,mm,ss
      endif
      i=index(logblock,'SZA = ')+6
      if(i.eq.6) stop 'Cannot find "SZA"'
      read(logblock(i:i+5),*)asza

c     Get other stuff from logblock, PCAT/PCDA dependent
      i=index(logblock,'PCAT')
      if(i.ne.0)then                 !PCAT spectrum, read aperture from log text
         pcat=.true.
         i=index(logblock,'Resoln = ')+9
         if(i.eq.9) stop 'Cannot find "Resoln"'
         read(logblock(i:i+9),*)resn
         i=index(logblock,'Solar_Aperture  = ')+18
         read(logblock(i:i+4),'(f5.1)')fovi
      else
         pcat=.false.
         i=index(logblock,'Res= ')+4
         if(i.eq.4) stop 'Cannot find "Res"'
         read(logblock(i:i+10),*)resn
         call julian(iy,im,id,jd)         !PCDA spectrum, determine aperture from date
         if(jd.lt.2450253)then                !18-Jun-1996
            if(ifil.le.5)then
               fovi=0.6
            else
               fovi=0.7
            endif
         elseif(jd.ge.2450253.and.jd.lt.2450527)then         !19-Mar-1996
            if(ifil.le.5)then
               fovi=0.6
            else
               fovi=1.0
            endif
         elseif(jd.ge.2450527)then                           !19-Mar-1997
            if(ifil.le.5)then
               fovi=0.8
            else
               fovi=1.0
            endif
         endif
      endif

      fovi=fovi/325.  !GFIT uses fov=angular diamter, not radius
      wavtkr=9900. ! solar tracker is pointed by a near-IR Si photodiode
      bytepw=-3    ! PC format I*4
      possp=544    ! 512 byte main header + 32 byte sub-file header
      apf='BX'     ! Boxcar apodisation
      opd=1./resn
      delwav=(nue-nus)/(npoints-1)
      delwav=delwav*15798.0138D0/lasf         !laser wavenumber correction
      delwav=delwav*(1.0D0+fovi*fovi/16.0D0)  !fov correction

      if(ivac.le.1)then
         pins=0.1
      else
         pins=pout
      endif

c     DG 4/7/2002
c     Following are the first and last points in the file
c     Replace by narrower limits to avoid GFIT trying to fit noise regions at edge of bandpass
c        ifirst=nint(nus/delwav)
c        ilast=ifirst+npoints-1
      ifirst=nint(sigflo(ifil)/delwav)
      ilast=nint(sigfhi(ifil)/delwav)
      possp=possp+4*(ifirst-nint(nus/delwav))
c!!!!!!!!!!!
      snr= 2000.0d0

c     DG070603 - Calculate the mean zpd time +mean of first and last zpdtime)
      ss=ss+(nfwdscans-1)*opd*(1.0/speed(ispeed)+0.5)/4.0
      gmt=dble(hh+(mm+ss/60.d0)/60.d0)
c
      return
      end
