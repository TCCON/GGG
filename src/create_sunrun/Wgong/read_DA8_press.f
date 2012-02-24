      subroutine read_DA8_press
     &     (specname,path,pout,pins) 
c     Read DA8-GRAMS format for Wollongong solar BOMEM spectra collected by SUNRUN.
c     DG Jan03
c     Revised Jun2007 as subroutine for new sunrun/runlog split

      implicit none
      include "../ggg_const_params.f"

      character   specname*(*),spfmt*2,path*(*),apf*2
      integer*4   i,hedlen,luns,iy,im,id,hh,mm,ss,
     &            possp,bytepw,ifirst,ilast,jd,nfwdscans,npoints,iend
      parameter   (hedlen=512,luns=19)

      real*8      d2r,delwav,nus,nue,fovi,asza,gmt,gmtstart,lasf,
     &            opd,snr,resn,wavtkr,pins,tins,hins,pout,tout,hout
      parameter   (d2r=dpi/180)

      character*(hedlen) chhedr
      integer*2   i2hedr(hedlen/2), i2bombin(512)
      integer*4   i4hedr(hedlen/4), i4bombin(256)
      real*4      r4hedr(hedlen/4), r4bombin(256)
      real*8      r8hedr(hedlen/8), r8bombin(128)
      real*4      sigflo(8), sigfhi(8)
      character*1024 logblock, bombin
      integer*4   logoffset, charoff, ifil, ivac, iostat
      logical*1   pcat


      equivalence (chhedr,i2hedr,i4hedr,r4hedr,r8hedr)	
      equivalence (bombin,i2bombin,i4bombin,r4bombin,r8bombin)

c     DG 4/7/2002  - define actual bandpass for NDSC filters to avoid GFIT
c     fitting noise regions at edge of bandpass.
      DATA sigflo /3950.0, 2960.0, 2350.0, 2050.0, 1850.0,
     &            700.0,  980.0,  700.0/
      DATA sigfhi /4450.0, 3500.0, 3120.0, 2600.0, 2250.0,
     &            1350.0, 1350.0, 1025.0/

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

c     Get stuff from Bomem binary block
      ivac=i2bombin(274)
      pout=i2bombin(319)/3276.8*260.0/5+800.0-0.1
      if (pout.le.950.0.or.pout.gt.1040.0)pout=1013.0

c     Note ivac is incorrect in PCAT1.2c - always=0
      if(ivac.le.1)then
          pins=0.1
      else
          pins=pout
      endif
c
      return
	end
