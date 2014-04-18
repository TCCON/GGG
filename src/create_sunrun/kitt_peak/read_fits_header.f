      subroutine read_fits_header(path,iend,dtype,nsp,nus,nue,iy,im,
     & id,hh,mm,ss,apt,dur,vel,apf,phr,res,lasf,nip,pkl,prl,
     & possp,oblat,oblon,obalt,tins,pins,hins,tout,pout,hout)
c
c  Subroutine to extract information from FITS headers.
c
c INPUTS:
c       path  C*(*)  Path to spectrum
c       iend  I*4  Endianess of host computer
c       dtype I*4  Data Type (spectrum=1031; interferogram=2055)

c OUTPUTS:
c       nsp,   I*4  ! Number of Spectral Points
c       nus,   R*8  ! First X-value
c       nue,   R*8  ! Last X-value
c       iy,    I*4  ! Year (4 digit)
c       im,    I*4  ! Month
c       id,    I*4  ! Day
c       hh,    I*4  ! Hour
c       mm,    I*4  ! Minute
c       ss,    I*4  ! Seconds
c       apt,   R*8  ! Aperture diameter (mm)
c       vel,   R*8  ! Bruker "Velocity" (laser fringe rate in kHz)
c       apf,   C*2  ! Apodization function (e.g. "BX")
c       phr,   R*8  ! Phase Resolution (cm-1) = 0.9/OPD_short
c       res,   R*8  ! Spectral Resolution (cm-1) = 0.9/OPD_long
c       lasf,  R*8  ! Laser WaveNumber (cm-1 ) = 15798.
c       nip,   I*4  ! Numner of Interferogram Points
c       pkl,   I*4  ! Interferogram Peak (ZPD) location: forward scan
c       prl,   I*4  ! Interferogram Peak (ZPD) location: reverse scan
c       possp  I*4  ! Position Of Starting Spectral Point (bytes)
c
      implicit none
      character path*(*),apf*2,hst*80,date*8,comment*40
      integer*4 dtype,hedlen,luns,iend,
     & npoints,ncenter,ic
      integer*4 nsp,iy,im,id,hh,mm,ss,nip,pkl,prl,possp
      parameter (hedlen=512,luns=19)

      real*8 vel,phr,res,apt,dur,oblat,oblon,obalt,
     & tins,pins,hins,tout,pout,hout,fovi,fovo,
     & amal,resn,snr,nus,nue,delwav,sza1,sza2,sampfreq,lasf,
     & tgmt,zpdtim,scandelt,opd,gmt

      integer*4 reclen,bytepw,irec,iflag,
     & iscan,nscan,krec,nintp,ifirst,ilast

c==================================================================
      if(1.eq.2)write(*,*)iend,dtype,nsp,apt,vel,phr,res,nip,pkl,prl
        oblat=31.958d0   ! Kitt Peak
        oblon=-111.595d0 ! Kitt Peak
        obalt=2.092d0     ! Kitt Peak
        fovi=0.002d0
        fovo=fovi/50d0
        amal=0.0d0
        bytepw=4
        apf='BX'
        resn=0.02  ! In case OPD value is not found in spectrum header
        snr=500.0d0
        pout=0.0d0 ! Initial value
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
                 read(hst(10:),'(a)') comment
              elseif(hst(:9).eq.'DAY     =') then
                 read(hst(10:),*) date
                 if(date(1:1).eq.' ') read(hst(10:),'(3x,a)') date
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
        iy=iy+1900     
        opd=0.5/resn
        if (pout.ne.0) then
           pout=1013.25*pout/760
        else
           pout=9999.99
        endif
        hout=40
        pins=0.
        tins=25.
        hins=10.
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

        write(*,'(a,f9.1,2x,a)') path(46:58),
     &  90*(sza1+sza2)/3.14159,comment
      return
      end
