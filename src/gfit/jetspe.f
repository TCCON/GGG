      subroutine jetspe(specpath,opd,graw,ifirst,ilast,possp,bytepw,
     & nus,nue,apo_m,interp,foff,res,slit,mii,nscycle,
     & yobs,mmp,nmp,nustrt,delwav,status)
c  Reads a portion of spectrum from disk, apodizes, interpolates, and resamples
c  the spectrum as requested, then places the result at the beginning of YOBS.
c  Note that YOBS is also used as workspace and so the elements of YOBS above
c  NMP might contain garbage.
c
c INPUTS:
c SPECPATH  C**  spectrum path
c      OPD  R*8  Maximum OPD
c     GRAW  R*8  Raw spectral point spacing
c   IFIRST  I*4
c    ILAST  I*4 
c    POSSP  I*4
c   BYTEPW  I*4
C      NUS: R*8  Starting frequency in cm-1
C      NUE: R*8  Ending frequency in cm-1
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
C  Warnings (i.e. requested spectral portion could only partially be read)
C  status =7  Parameter MSF too small: Convolution operator was truncated 
C         =5  Requested spectral portion ends after last disk value
C         =4  Requested spectral portion starts before first disk value
C         =2  YOBS(MMP) cannot hold all interpolated values; increase MMP.
C         =1  YOBS(MMP) cannot hold all raw spectral values; increase MMP.
C
C  STATUS= 0  no error, everything worked OK !
C
C  Fatal error return codes:
C  STATUS=-1  Spectrum file could not be found
C        =-2  No overlap between requested interval and disk file.
C             This can also happen if MMP=0, or if NUS=NUE
C        =-3  Requested interval is too close to edge of disk file to perform
C             convolution. Can only be returned completely with APO=0, INTERP=0
C======================================================================
      implicit none
      character specpath*(*)

      INTEGER*4 apo_m,possp,nmp,mpts,STATUS,mii,MMP,KINTPC,
     & k1,M1,M2,I1,I2,k,bytepw,IFIRST,ilast,INTERP,iskip,iabpw,
     & nhw,nsf,nii,nele,nscycle,j
c      parameter (mii=50753,nscycle=18)  ! max dimension of slit function 
c
      REAL*8 dzero,fr,resnog,resn,rect,opd
      REAL*8 nus,nue,graw,delwav,nustrt
c
      REAL*4 slit(mii),yobs(mmp),foff,res,unity,tot
      parameter (dzero=0.0d0,unity=1.0)
c
      status=0
      rect=0.0d0
      resn=0.5d0/opd
      resn=dmax1(resn,dble(res))
c      if(resn.lt.graw) resn=graw
      resnog=resn/graw
c      write(*,*)opd,resn,graw,resnog
      nhw=nint(nscycle*resnog)
      nsf=2*nhw+1
c-------------------------------------------------------------------------
c  Determine the indexes of the first (m1) and last (m2) raw spectral
c  points in the requested spectral region (nus to nue)
c     KINT replaced by KINTPC to avoid conflict with instrisic function KINT     !DG000831
      m1=1+KINTPC(nus/graw)
      m2=KINTPC(nue/graw)
c-------------------------------------------------------------------------
c  Check that the spectral interval m1-nhw to m2+nhw is present on disk file
      if(m1-nhw .lt. ifirst) then
         status=4
         m1=nhw+ifirst
      endif 
c
      if(m2+nhw .gt. ilast) then
         status=5
         m2=ilast-nhw
      endif 
      mpts=M2-M1+nsf  ! number of points to be read from disk
c-------------------------------------------------------------------------
c  Check that MPTS > 0 to avoid attempting a zero length read
      if(mpts .le. 0) then
         status=-2                  ! zero length read attempted
         write(*,89)specpath,ifirst*graw,ilast*graw
89       format(a48,' only encompasses',f9.3,' to ',f9.3,' cm-1')
         return
      endif
c-------------------------------------------------------------------------
c  Check that MPTS > NS to avoid a futile read
      if(mpts .lt. nsf) then
         write(*,*)' Insufficient overlap to fill ILS'
         status=-3                  ! attempted read will not even fill ILS
         return
      endif
c-------------------------------------------------------------------------
c  Check that original spectra values will not overflow YOBS(MMP)
      if( mpts .gt. mmp ) then
         write(*,*)'HETSPE warning: increase MMP to',mpts
         status=1
         mpts=mmp
         m2=m1+mpts-nsf
      endif
c-------------------------------------------------------------------------
c  Check that interpolated values will not overflow YOBS(MMP)
      if( 1+interp*(mpts-nsf) .gt. mmp ) then
      write(*,*)'IETSPE warning: increase MMP to',1+interp*(mpts-nsf)
         status=2
         mpts=nsf+(mmp-1)/interp
         m2=m1+mpts-nsf
      endif
c-------------------------------------------------------------------------
c  Fetch the raw spectral points from disk file. Note that if INTERP > 1
c  the raw spectral values are not placed at the beginning of YOBS. 
c  Instead they fill addresses starting at K1. This allows the convolution
c  to be performed "in place" without prematurely overwriting any points.
      k1=1+(interp-1)*(mpts-nsf)
      iabpw=iabs(bytepw)
      if(iabpw.eq.3) iabpw=4
      iskip=m1-nhw-ifirst+possp/iabpw
      call fetch(specpath,bytepw,iskip,yobs(k1),mpts)
c      write(*,*)(yobs(j),j=k1,k1+mpts-1)
c-------------------------------------------------------------------------
      i1=1+KINTPC(interp*nus/graw)
      if(i1.le.interp*(m1-1)) i1=interp*m1
      i2=KINTPC(interp*nue/graw)
      if(i2.ge.interp*(m2+1)) i2=interp*m2
      nmp=i2+1-i1
c      write(*,*)graw,delwav
c      write(*,*)nus,nue,i1,i2,m1,m2
      delwav=graw/interp
      nustrt=i1*delwav
c-------------------------------------------------------------------------
c  Calculate slit function and normalize each subset to unity independently
c  RECTOG is zero since we are dealing with measured spectra.
      fr=dble(foff)/delwav
      nii=1+2*interp*nhw
      if(nii.gt.mii) write(*,*)' FETSPE: Increase parameter MII to',nii
      call profzl(apo_m,nii,interp*resnog,dzero,fr,slit)
      nele=nii/interp
      do k=interp,1,-1
         call vdot(slit(k),interp,unity,0,tot,nele)
         call vmul(slit(k),interp,1./tot,0,slit(k),interp,nele)
      end do   ! k=1,interp
      slit(nii)=slit(nii)/tot
c
C  Perform the convolution that does the shifting, interpolation & apodization
c  Convolve spectral values with pre-normalized operator A (sinc function)
c      write(*,*)'IETSPE calling NEWDEC..',apo_m,resnog,mpts,nmp,nhw,nii
c      write(*,*) (slit(j),j=1,nii)
      call newdec(yobs(k1),mpts,slit,nii,interp,1.d0/interp,0.d0,
     & yobs,nmp)
      return
      end
C
c   KINT replaced by KINTPC to avoid instrisic function KINT    !DG000831
      function KINTPC(x)
      integer*4 kintPC
      real*8 x
      KINTPC=idint(x)
      if (x.le.0.0D0) KINTPC=KINTPC-1
      return
      end
