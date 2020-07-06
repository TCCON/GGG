      subroutine jetspe(specpath,resn,graw,ifirst,ilast,possp,bytepw,
     & nus,nue,apo_m,interp,foff,res,yobs,mip,nip,nustrt,delwav,status)
c  Reads a portion of spectrum from disk, apodizes, interpolates, and resamples
c  the spectrum as requested, then places the result at the beginning of YOBS.
c  Note that YOBS is also used as workspace and so the elements of YOBS above
c  NMP might contain garbage.
c
c INPUTS:
c SPECPATH  C**  spectrum path
c      OPD  R*8  Maximum OPD (cm)
c     GRAW  R*8  Raw spectral point spacing of measured spectrum (cm-1)
c   IFIRST  I*4
c    ILAST  I*4 
c    POSSP  I*4
c   BYTEPW  I*4
C      NUS: R*8  Starting frequency of window (cm-1)
C      NUE: R*8  Ending frequency of window (cm-1)
C     FOFF: R*4  Frequency offset (cm-1) by which to resample spectrum
C      APO: I*4  Desired apodization to be applied to the data
C   INTERP: I*4  Interpolation factor; ratio output/input points
C      RES: R*4  Spectral resolution (0.5/OPD) at which spectrum should appear
C                If RES < actual resolution spectrum will appear at actual resn
C      MMP: I*4  Maximum dimension of array YOBS
C
C OUTPUTS:
C     YOBS: R*4  Array of spectral values
C      NMP: I*4  The number of spectral points found in the requested region
C   nustrt: R*8  Frequency (cm-1) of the first point in array YOBS
C   DELWAV: R*8  The frequency spacing in cm-1 between returned points
C   STATUS: I*4  Error return code
C      RES: R*4  Spectral resolution (0.5/OPD) at which spectrum appears
C
C   status=-1  Requested spectral portion exceeds frequency limits of disk file
C         =-2  NPTS =< 0   Zero length read attempted (e.g. NUS=NUE)
C         =-3  NPTS < NSF  Cannot convolve with ILS (e.g. NUS=NUE)
C         =-4  NPTS > MMP  YOBS(MMP) cannot hold all raw spectral values; increase MMP.
C         =-5  YOBS(MMP) cannot hold all interpolated values; increase MMP.
C          -6  NII > MII   Increase size of parameter MII
C
C  STATUS= 0  no error, everything worked OK !
C
C  Error return codes:
C        =-2  No overlap between requested interval and disk file.
C             This can also happen if MMP=0, or if NUS=NUE
C======================================================================
      implicit none

      character specpath*(*)

      INTEGER*4 apo_m,possp,nip,nmp,STATUS,msi,mip,nexpl,nexpr,
c     & j,
     & k1,m1,m2,i1,i2,bytepw,ifirst,ilast,interp,iskip,iabpw,
     & nhw,nsf,nscycle,j1,j2,ii
      parameter (msi=50753,nscycle=25)  ! max dimension of vsinc function 
c
      REAL*8 dzero,fr,resnog,resn,rect,vbar,hwid,dd,sh
      REAL*8 nus,nue,graw,delwav,nustrt
c
      REAL*4 vsinc(msi),yobs(mip),foff,res,dum
      parameter (dzero=0.0d0)

      if(nus.ge.nue) stop 'NUS >= NUE'
      status=0
      rect=0.0d0
      resn=dmax1(resn,dble(res))
      dd=nscycle*resn  ! half-wifth of ILS in cm-1
c      if(resn.lt.graw) resn=graw
      resnog=dabs(resn/graw)
c      write(*,*)'jetspe: specpath=',specpath(:80)
c      write(*,*)'jetspe: resn,graw,resnog=',resn,graw,resnog
      nhw=nint(nscycle*resnog)
      nsf=2*iabs(nhw)+1
c-------------------------------------------------------------------------
c  Determine the indexes of the first (m1) and last (m2) raw spectral
c  points in the requested spectral region (nus to nue)
      if(graw.gt.0) then
         m1=1+int(nus/graw)
         m2=int(nue/graw)
      else    ! if graw is -ve, switch m1 & m2
         m2=-1+int(nus/graw)
         m1=int(nue/graw)
      endif
c  Note that it is always true that M2>M1 irrespective of whether they are +ve or -ve
c-------------------------------------------------------------------------
      vbar=0.5d0*graw*(ilast+ifirst)
      hwid=0.5d0*dabs(graw*(ilast-ifirst))
c      write(*,*) 'jetspe: ',ifirst,ilast,graw,vbar,hwid

c  Check that the lowest required measurement frequency is in the spectrum file
      if( nus-dd .lt. vbar-hwid ) then
c         write(*,*) 'Lower window limit < disk file:',nus-dd,vbar-hwid
         status=-1
         return
      endif 
c
c  Check that the highest required measurement frequency is in the spectrum file
      if( nue+dd .gt. vbar+hwid ) then
c         write(*,*) 'Upper window limit > disk file:',nue+dd,vbar+hwid
         status=-1
         return
      endif 
      nmp=m2-m1+nsf  ! number of measured points to be read from disk
c-------------------------------------------------------------------------
c  Check that MPTS > 0 to avoid attempting a zero length read
      if(nmp .le. 0) then
         status=-2                  ! zero length read attempted
         write(*,89)specpath,vbar-hwid,vbar+hwid
89       format(a48,' only encompasses',f9.3,' to ',f9.3,' cm-1')
         return
      endif
c-------------------------------------------------------------------------
c  Check that MPTS > NS to avoid a futile read
      if(nmp .lt. nsf) then
         write(*,*)' Insufficient overlap to fill ILS'
         status=-3                  ! attempted read will not even fill ILS
         return
      endif
c-------------------------------------------------------------------------
c  Check that original spectra values will not overflow YOBS(MMP)
      if( nmp .gt. mip ) then
         write(*,*)'JETSPE warning: increase MMP to',nmp
         status=-4
         return
      endif
c-------------------------------------------------------------------------
c  Check that interpolated values will not overflow YOBS(MMP)
      if( 1+interp*(nmp-nsf) .gt. mip ) then
         write(*,*)'JETSPE warning: increase MMP to ',1+interp*(nmp-nsf)
         status=-5
         return
      endif
c-------------------------------------------------------------------------
c  Fetch the raw spectral points from disk file. Note that if INTERP > 1
c  the raw spectral values are not placed at the beginning of YOBS. 
c  Instead they fill addresses starting at K1. This allows the convolution
c  to be performed "in place" without prematurely overwriting any points.
      k1=1+(interp-1)*(nmp-nsf)
      iabpw=iabs(bytepw)
      if(iabpw.eq.3) iabpw=4
      if(iabpw.eq.6) iabpw=8  ! ACE binary [Real, Imag]
      if(iabpw.eq.11) then
         iskip=m1-iabs(nhw)-ifirst
      else
         iskip=possp+iabpw*(m1-iabs(nhw)-ifirst)
      endif
      if(iskip.lt.0) then 
         write(*,*)'jetspe: IFIRST,ILAST,M1,NHW = ',IFIRST,ILAST,M1,NHW
         stop 'JETSPE: iskip < 0    IFIRST > ILAST ?'
      endif
c      write(*,*)'jetspe: bpw, m1, nhw, ifirst=',iabpw, m1, nhw, ifirst
c      write(*,*)'jetspe: possp,nmp,nsf,k1,iskip=',possp,nmp,nsf,k1,iskip
      call fetch(specpath,bytepw,iskip,yobs(k1),nmp)
      if(graw.lt.0.0) then
         j1=k1
         j2=k1+nmp-1
         do ii=1,nmp/2
            dum=yobs(j1)
            yobs(j1)=yobs(j2)
            yobs(j2)=dum
            j1=j1+1
            j2=j2-1
         end do
      endif
c-------------------------------------------------------------------------
c  i1/i2 are the interpolated or output spectrum starting/ending indices
      delwav=dabs(graw)/interp
      i1=1+int(nus/delwav)
      i2=int(nue/delwav)
      nip=i2+1-i1       !  number of interpolated spectral po
      nustrt=i1*delwav
c      write(*,*)nus,nue
c      write(*,*)i1,i2,nip,delwav,i1*delwav,i2*delwav

c      i1=1+int(interp*nus/dabs(graw))
c      i2=int(interp*nue/dabs(graw))
c      nip=i2+1-i1   ! Number of Interpolated Points
c      delwav=dabs(graw)/interp
c      nustrt=i1*delwav
c-------------------------------------------------------------------------
c  Calculate vsinc function and normalize each subset to unity independently
c  RECTOG is zero since we are dealing with measured spectra.
      fr=dble(foff)/delwav
      if(1+2*interp*iabs(nhw).gt.msi) then
         write(*,*)' JETSPE: Increase MII to',1+2*interp*iabs(nhw)
         status=-7
         return
      endif
c      write(*,*)'jetspe: interp,nhw=',interp,nhw
      call compute_ils(apo_m,nhw,interp,interp*resnog,dzero,fr,
     & vsinc)

C  Perform the convolution that does the shifting, interpolation & apodization
c  Convolve spectral values with pre-normalized operator A (sinc function)
      if(m1.gt.0) then
         sh=dfloat(i1-interp*m1)/interp
      else
         sh=dfloat(i1+interp*m2)/interp
      endif


c     write(*,*) 'JETSPE: Calling regrid2:',m1,nmp,k1,interp,nhw
      call regrid2(m1,nmp,yobs(k1),nhw,vsinc,interp,1.d0/interp,
     & m1+interp*nhw,nip,yobs,nexpl,nexpr)
c      write(*,*) 'JETSPE: Called regrid2:'
      if(nexpl.ne.0) write(*,*)'Warning: JETSPE: NEXPL=',nexpl
      if(nexpr.ne.0) write(*,*)'Warning: JETSPE: NEXPR=',nexpr
      return
      end
