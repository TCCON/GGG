      subroutine jetspe(specpath,opd,graw,ifirst,ilast,possp,bytepw,
     & nus,nue,apo_m,interp,foff,res,
     & yobs,mmp,nmp,nustrt,delwav,status)
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

      INTEGER*4 apo_m,possp,nmp,npts,STATUS,mii,MMP,KINTPC,
     & k1,m1,m2,i1,i2,k,bytepw,IFIRST,ilast,INTERP,iskip,iabpw,
     & nhw,nsf,nii,nele,nscycle
      parameter (mii=50753,nscycle=25)  ! max dimension of slit function 
c
      REAL*8 dzero,fr,resnog,resn,rect,opd,vbar,hwid,dd
      REAL*8 nus,nue,graw,delwav,nustrt
c
      REAL*4 slit(mii),yobs(mmp),foff,res,unity,tot
      parameter (dzero=0.0d0,unity=1.0)
c
      if(nus.ge.nue) stop 'NUS >= NUE'
      status=0
      rect=0.0d0
      resn=0.5d0/opd
      resn=dmax1(resn,dble(res))
      dd=nscycle*resn  ! half-wifth of ILS in cm-1
c      if(resn.lt.graw) resn=graw
      resnog=resn/graw
c      write(*,*)opd,resn,graw,resnog
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
      npts=M2-M1+nsf  ! number of points to be read from disk
c-------------------------------------------------------------------------
c  Check that MPTS > 0 to avoid attempting a zero length read
      if(npts .le. 0) then
         status=-2                  ! zero length read attempted
         write(*,89)specpath,vbar-hwid,vbar+hwid
89       format(a48,' only encompasses',f9.3,' to ',f9.3,' cm-1')
         return
      endif
c-------------------------------------------------------------------------
c  Check that MPTS > NS to avoid a futile read
      if(npts .lt. nsf) then
         write(*,*)' Insufficient overlap to fill ILS'
         status=-3                  ! attempted read will not even fill ILS
         return
      endif
c-------------------------------------------------------------------------
c  Check that original spectra values will not overflow YOBS(MMP)
      if( npts .gt. mmp ) then
         write(*,*)'JETSPE warning: increase MMP to',npts
         status=-4
         return
      endif
c-------------------------------------------------------------------------
c  Check that interpolated values will not overflow YOBS(MMP)
      if( 1+interp*(npts-nsf) .gt. mmp ) then
      write(*,*)'JETSPE warning: increase MMP to',1+interp*(npts-nsf)
         status=-5
         return
      endif
c-------------------------------------------------------------------------
c  Fetch the raw spectral points from disk file. Note that if INTERP > 1
c  the raw spectral values are not placed at the beginning of YOBS. 
c  Instead they fill addresses starting at K1. This allows the convolution
c  to be performed "in place" without prematurely overwriting any points.
      k1=1+(interp-1)*(npts-nsf)
      iabpw=iabs(bytepw)
      if(iabpw.eq.3) iabpw=4
      iskip=possp+iabpw*(m1-iabs(nhw)-ifirst)
      call fetch(specpath,bytepw,iskip,yobs(k1),npts)
      if(graw.lt.0.0)
     & call vswap(yobs(k1),1,yobs(k1+npts-1),-1,npts/2)
c-------------------------------------------------------------------------
      i1=1+int(interp*nus/dabs(graw))
      i2=int(interp*nue/dabs(graw))
      nmp=i2+1-i1
c      write(*,*)graw,delwav
c      write(*,*)nus,nue,i1,i2,m1,m2
      delwav=dabs(graw)/interp
      nustrt=i1*delwav
c-------------------------------------------------------------------------
c  Calculate slit function and normalize each subset to unity independently
c  RECTOG is zero since we are dealing with measured spectra.
      fr=dble(foff)/delwav
      nii=1+2*interp*iabs(nhw)
      if(nii.gt.mii) then
         write(*,*)' JETSPE: Increase parameter MII to',nii
         status=-7
         return
      endif
      call profzl(apo_m,nii,interp*dabs(resnog),dzero,fr,slit)
      nele=nii/interp
      do k=interp,1,-1
         call vdot(slit(k),interp,unity,0,tot,nele)
         call vmul(slit(k),interp,1./tot,0,slit(k),interp,nele)
      end do   ! k=1,interp
      slit(nii)=slit(nii)/tot
c
C  Perform the convolution that does the shifting, interpolation & apodization
c  Convolve spectral values with pre-normalized operator A (sinc function)
c      write(*,*)'JETSPE NEWDEC..',apo_m,resnog,npts,nmp,nhw,nii,k1
c      write(*,*) (slit(k),k=1,nii)
c      if(apo_m.gt.0) then
         call newdec(yobs(k1),npts,slit,nii,interp,1.d0/interp,0.0d0,
     &   yobs,nmp)
c      else
c         call vmov(yobs(k1+nhw),1,yobs,1,nmp)
c      endif
      return
      end
