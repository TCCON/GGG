c  avg_ker.f   GCT 2005-03-03
c  Program to compute averaging kernels from the zXXXXXXX
c  output files written by the ak.f subroutine called by GFIT.
c
c  Current version is a simplified version of aker.f.
c  It was simplified to try to troubleshoot some problems.
c  This new version no longer has the capability to handle
c  multiple target gases, previously used to kluge profile
c  retrievals.
c
c  The following improvements still need to be made:
c  1) Augment the A and B matrices in the equation Ax=B
c   to include the a priori information. This will require
c   expanding the first dimensions of arrays A and B from
c   MMP to MMP+MFP. Since the constraint on the VF of the
c   first target gas is extremely weak, I don't expect that
c   this will make much difference, which is why I haven't
c   been motivated to fix it.
c
c  2) This routine currently assumes that the CL, CT and FS
c   are always fitted. But this may not be true. And if ZO
c   is retrieved, this case needs to be handled. So we
c   need a more flexible way of handling CL, CT, FS, & ZO.
c
c  3) A thought:  The subroutine AK.F  and FM.F are very
c  similar.  If FM were equipped with the option of
c  calculating the single-level PDs, then we wouldn't need
c  the AK subroutine at all. [Done, 2007]
c
c  4) Another thought: It might be more convenient to
c  implement the averaging kernel calculation inside GFIT,
c  rather than as a stand-alone program. This is because all
c  the a priori information is already there.
c  I hesitated doing this in the past because it requires
c  an extra array of dimension (MMP+MFP,MLEV) to hold
c  the single-level PDs. The current version of GFIT has
c  MMP=360000, MLEV=150, so that's 216 Mbyte of memory for
c  an array that will typically be unused.
c
c  The other option would be to compute the averaging kernels
c  layer-by-layer, so that the extra array would be (360000,1).
c  But this would be much slower.
c---------------------------------------------------
c
c  AVG_KER computes the averaging kernels by solving the matrix equation
c     A.x=b
c  using my favorite subroutine HFTI
c
c  On input, matrix A(NMP,NFP) contains the PD's of all retrieved
c  quantities as a function of frequency.  The NFP retrieved parameters
c  include the VSF of the various target gases and the CL, CT, FS, ZO,
c  (but not channel fringes).
c
c  On input, matrix B(NMP,NLEV) contains the single-level PDs of the
c  first target gas. So these are the effect on the calculated spectrum
c  of scaling the vmr at a particular level.  Subroutine HFTI solves
c  for the multiple right-hand sides (multiple levels) in a single call.
c
c  On output, B(NFP,NLEV) contains the impacts on the NFP fitted parameters
c  of adjusting the first target gas at the NLEV altitude levels.
c  The averaging kernels of the first target gas are to be found in
c  the first row of B, i.e. B(1,NLEV), which are then normalized using
c  the effective pressure thicknesses of the layers.

c  NMP is the number of spectral points in the fitted window
c  NFP is the number of fitted/retrieved parameters
c  NLEV is the number of atmospheric levels.
c
c  2009-08-31  Added the ZO (zero offset) option to the calculation of the kernels.
c  This option required corresponding changes to the GFIT code, and so is not
c  backward-compatible.
c
      implicit none
      integer*4 lunr,luns,lunw,mmp,nmp,imp,mfp,
     & ntg,itg,mlev,nlev,ilev,nfp,fbc,fnbc,lnbc,lf,i,ispe,nlhead
      parameter (lunr=12,luns=13,lunw=14,mmp=32000,mfp=8,
     & mlev=100)
      integer*4 krank,ip(mfp)
      real*4 a(mmp,mfp),b(mmp,mlev),work(mfp),rnorm(mlev),
     & tau,psc(mlev),pres(mlev),ps,pwas,tsc,
     & ak1,ak2,ak,akwas,tak,tb
      character colfile*128,filename*128,akpath*128,spectrum*80,
     & version*60

      version=' avg_ker   Version 1.2.0    2009-08-31   GCT'

      tau=1.e-7

c  Find out what window the AK's are to be calculated for.
      write(*,*)'Enter name of most recent .col file'
      read(*,'(a)') colfile
      open(luns,file=colfile,status='old')
      read(luns,*) nlhead
      do i=2,nlhead-4
         read(luns,'(34x,a)') akpath
      end do
      do i=nlhead-3,nlhead
         read(luns,*)
      end do

      write(*,*)'   Spectrum Path/Name              '//
     &'         p-averaged AK   c-averaged AK'
c  Main loop over spectra.
      do ispe=1,99999
      read(luns,'(a)',end=99) spectrum
      filename=akpath(:lnbc(akpath))//spectrum(fnbc(spectrum):)
      lf=fbc(filename)-1
c
c  Read total column partial differentials
      open (lunr, file=filename(:lf), status='old')
      read(lunr,*)nmp,ntg
      if(nmp.gt.mmp) then
         write(*,*)'NMP, MMP = ',nmp,mmp
         stop 'increase MMP'
      endif
      nfp=ntg+4
      if(nfp.gt.mfp) stop 'increase MFP'
      do itg=1,ntg
         read(lunr,*) (a(imp,itg),imp=1,nmp)
      end do
c
c  Read single-level partial differentials of first target gas.
      read(lunr,*)nmp,nlev
      if(nmp.gt.mmp) stop 'increase MMP'
      if(nlev.gt.mlev) stop 'increase MLEV'
      do ilev=1,nlev
         read(lunr,*) (b(imp,ilev),imp=1,nmp) ! Single-Level PDs of First Target gas
      end do
c
      read(lunr,*) (a(imp,ntg+1),imp=1,nmp)  ! CL partial differential
      read(lunr,*) (a(imp,ntg+2),imp=1,nmp)  ! CT partial differential
      read(lunr,*) (a(imp,ntg+3),imp=1,nmp)  ! FS partial differential
      read(lunr,*) (a(imp,ntg+4),imp=1,nmp)  ! ZO partial differential
      read(lunr,*) (psc(ilev),ilev=1,nlev)   ! partial slant columns
      read(lunr,*)  ps                       ! surface pressure (atm)
      read(lunr,*) (pres(ilev),ilev=1,nlev)  ! pressure
      close(lunr)
      call vdot(psc,1,1.0,0,tsc,nlev)        ! total slant column
c
c  Write out the PD's in a form that can be easily plotted (e.g. xyplot)
c  This is for trouble-shooting/illustrative purposes only.
      open(lunw,file=filename(:lf)//'.wtf', status='unknown')
      write(lunw,*)2,1+nfp+nlev
      write(lunw,'(a18,999(9x,a1,i2.2))')
     & 'i  CL  CT  FS  ZO ',
     & ('T',itg,itg=1,ntg),
     & ('S',ilev,ilev=0,nlev-1)
      do imp=1,nmp
         write(lunw,'(i5,999(1pe12.4))') imp,
     &    (a(imp,itg),itg=ntg+1,ntg+4),  ! CL, CT, FS, ZO
     &    (a(imp,itg),itg=1,ntg),
     &    (b(imp,ilev),ilev=1,nlev)
      end do
      close(lunw)

c  Solve the equation A.x=b
      call shfti(a,mmp,nmp,nfp,b,mmp,nlev,tau,krank,rnorm,work,ip)
      if(krank.lt.nfp) write(*,*)' Rank Deficiency: ',krank,nfp
c
c  Write out the Averaging Kernels.
c  Also, calculate the pressure- and density-weighted kernels.
      open(lunw,file=filename(:lf)//'.aks', status='unknown')
      write(lunw,*)3,3
      write(lunw,*) version
      write(lunw,*)  'Level   AK   Pressure_(atm)'
      ak1=b(1,1)*tsc/psc(1)
      ak2=b(1,2)*tsc/psc(2)
c FIX ME: the following line assumes that Ps lies between p(1) and p(2)
      akwas=(ak1*(ps-pres(2))+ak2*(pres(1)-ps))/(pres(1)-pres(2))
      pwas=ps
      tak=0.0
      tb=0.0
      do ilev=1,nlev
         tb=tb+b(1,ilev)
         ak=b(1,ilev)*tsc/psc(ilev)
         write(lunw,*)ilev-1,ak,pres(ilev)
         tak=tak+0.5*(ak+akwas)*(pwas-pres(ilev))
         akwas=ak
         pwas=pres(ilev)
      end do
      close(lunw)
c      write(*,'(a,f9.6)')filename(:lf),tak/pres(1)
      write(*,'(a,2f15.5)')filename(:lf),tak/ps,tb
      end do ! ispe=1,999999
      write(*,*) 'Warning: Loop limit exceeded'
99    close(luns)
      write(*,*)'Number of Measured Points =',nmp
      write(*,*)'Number of Target Gases    =',ntg
      write(*,*)'Number of Model Levels    =',nlev
      stop
      end
