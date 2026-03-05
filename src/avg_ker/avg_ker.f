c  avg_ker.f   GCT 2005-03-03
c  Program to compute averaging kernels from the j_XXXXXXX
c  output files written by the fm.f and do_retrieval.f
c  subroutines called by gfit.f
c
c  Current version is a simplified version of aker.f.
c  It was simplified to try to troubleshoot some problems.
c  This new version no longer has the capability to handle
c  multiple target gases, previously used to kluge profile
c  retrievals.
c
c  The following improvements still need to be made:
c  1) Augment the A and B matrices in the equation Ax=B
c  to include the a priori information. This will require
c  expanding the first dimensions of arrays A and B from
c  MMP to MMP+MFP. Since the constraint on the VF of the
c  first target gas is extremely weak, I don't expect that
c  this will make much difference, which is why I haven't
c  been motivated to fix it. [Done, Mar 2017]
c
c  2) This routine currently assumes that the CL, CT, CC, FS, 
c  and ZO are always fitted. But this may not be true.
c
c  3) A thought:  The subroutine AK.F  and FM.F are very
c  similar.  If FM were equipped with the option of
c  calculating the single-level PDs, then we wouldn't need
c  the AK subroutine at all. [Done, 2007]
c
c  4) Another thought: It might be more convenient to
c  implement the averaging kernel calculation inside GFIT,
c  rather than as a stand-alone program. This is because
c  all the a priori information is already there.
c  This is true for single windows, but this wouldn't help
c  when computing the average kernel over multiple windows.
c  I hesitated doing this in the past because it requires
c  an extra array of dimension (MMP+MFP,MLEV) to hold
c  the single-level PDs, which are currently written out
c  to the ak file one by one. The current version of GFIT
c  has MMP=360000, MLEV=150, so that's 216 Mbyte of memory
c  for an array that will typically be unused.
c
c  The other option would be to compute the averaging kernels
c  layer-by-layer, so that the extra array would be (360000,1).
c  But this would be much slower.
c---------------------------------------------------
c
c  AVG_KER computes the averaging kernels by solving the matrix
c  equation
c     A.x=B
c
c  On input, matrix A(NMP,NFP) contains the PD's of all retrieved
c  quantities as a function of wavenumber index.  The NFP retrieved
c  parameters include the VSF of the various target gases and
c  CL, CT, CC, ... FS, SG, ZO, (but not channel fringes).
c
c  On input, matrix B(NMP,NLEV) contains the single-level PDs of the
c  first target gas. So these represent the effect on the calculated
c  spectrum of scaling the vmr at a particular level. Subroutine HFTI
c  solves for the multiple right-hand sides (multiple levels) in a
c  single call.
c
c  On output X(NFP,NLEV), which is stored in B, contains the impacts
c  on the NFP fitted parameters of adjusting the first target gas at
c  the NLEV altitude levels. The averaging kernels of the first target
c  gas are to be found in the first row of B, i.e. B(1,NLEV), which are
c  then normalized using the effective pressure thicknesses of the layers.

c  NMP is the number of spectral points in the fitted window
c  NFP is the number of fitted/retrieved parameters
c  NLEV is the number of atmospheric levels.
c
c  2009-08-31  Added the ZO (zero offset) option to the calculation of the kernels.
c  This option required corresponding changes to the GFIT code, and so is not
c  backward-compatible.
c
c  2017-03-15  Added ap_flag to facilitate comparison of kernels
c  computed with and without the a priori constraints.

      implicit none
      include "../gfit/ggg_int_params.f"
      
      integer*4 lunr_jac,lunw_aks,lunw_akall,
     & mmp,nmp,imp,mfp,j,ap_flag,lj,lrt,ls,lak,
     & ntg,jtg,nlev,ilev,jlev,nfp,lnbc,ispec,idum
      parameter (lunr_jac=12,lunw_aks=15,
     & lunw_akall=16,mmp=90000,mfp=24)
      integer*4 krank,ip(mfp)
      real*4 a(mmp+mfp,mfp),b(mmp+mfp,mlev),work(mfp),rnorm(mlev),
     & tau,psc(mlev),z(mlev),pres(mlev),ps_atm,pwas,tsc,
     & rmsocl,zmin,sza,airmass,
     & ak1,ak2,ak,akwas,tak,tb
      character jacfile*512,akpath*150,specname*128,
     & version*60,winfo*128,gggdir*128,dl*1

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol   ! Avoid compiler warning (unused parameter)
      idum=mcolvav   ! Avoid compiler warning (unused parameter)
      idum=mgas      ! Avoid compiler warning (unused parameter)
      idum=mlev      ! Avoid compiler warning (unused parameter)
      idum=mrow_qc   ! Avoid compiler warning (unused parameter)
      idum=mspeci    ! Avoid compiler warning (unused parameter)
      idum=mvmode    ! Avoid compiler warning (unused parameter)
      idum=ncell     ! Avoid compiler warning (unused parameter)
      idum=nchar     ! Avoid compiler warning (unused parameter)

      version=' avg_ker   Version 2.11-alpha     2021-02-19   GCT'
      write(*,*) version
      tau=1.e-7
      ap_flag=1   ! Include a priori constraints
      ap_flag=0   ! Ignore a priori constraints

      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir

c  Open .jac file and read header.
      if (iargc() == 0) then
         write(*,*)'Enter name of relevent .jac file'
         read(*,'(a)') jacfile
      elseif (iargc() == 1) then
         call getarg(1, jacfile)
      else
         stop 'Usage: $gggpath/bin/avg_ker jacfile'
      endif
      lj=lnbc(jacfile)
c  JLL 2021-03-04: I'm deliberately keeping this different from
c  Geoff's code - his assumes that the .jac files are always in
c  $GGGPATH/jac, mine will allow it to be anywhere.
      open(lunr_jac,file=jacfile,status='old')
      read(lunr_jac,'(a)')
      read(lunr_jac,'(a)') winfo

      write(*,*)'   Spectrum Path/Name            '//
     &' Ps       P-averaged_AK   C-averaged_AK'

      akpath=gggdir(:lrt)//'ak'//dl//'k'//char(48+ap_flag)
     &//jacfile(2:lj)
      lak=lnbc(akpath)
      open(lunw_akall,file=akpath(:lak)//'.all',status='unknown')
      write(lunw_akall,*) 3,6
      write(lunw_akall,*) version
      write(lunw_akall,*)' ispec zmin sza airmass z ak p'

c  Main loop over spectra.
      do ispec=1,99999
c
c  Loop over different windows of same gas.
c  Read total column partial differentials of first target gas.
         read(lunr_jac,'(a)',end=99)specname
         ls=lnbc(specname)
         write(*,*)'ls,specname = ',ls, specname(:ls)
         read(lunr_jac,*,end=99)zmin,sza,rmsocl,airmass
         write(*,*)'zmin,sza,rmsocl,airmass = ',zmin,sza,rmsocl,airmass
         read(lunr_jac,*)nmp,ntg,nfp,nlev
         write(*,*)'nmp,ntg,nfp = ',nmp,ntg,nfp
         if(nmp.gt.mmp) stop 'increase parameter MMP'
         if(nfp.gt.mfp)  stop 'increase parameter MFP'
         if(nlev.gt.mlev) stop 'increase parameter MLEV'

c  Read column Jacobians of target gases into A
c  Set A to zero for elements A(nmp+1,jtg) and beyond.
         do jtg=1,ntg
            read(lunr_jac,*) (a(imp,jtg),imp=1,nmp)
            call vmov(0.0,0,a(nmp+1,jtg),1,nfp)
         end do

c  Read single-level partial differentials of first target gas.
c  Set b to zero for elements above b(nmp,*)
         do ilev=1,nlev
            read(lunr_jac,*) (b(imp,ilev),imp=1,nmp) ! Single-Level PDs of 1st Target gas
            call vmov(0.0,0,b(nmp+1,ilev),1,nfp)
         end do
c
c  Read continuum (and FS, SG, ZO) PDs
         do j=ntg+1,nfp
            read(lunr_jac,*) (a(imp,j),imp=1,nmp)
            call vmov(0.0,0,a(nmp+1,j),1,nfp)
         enddo

c  Read other parameters, including a priori constraint
         read(lunr_jac,*) (psc(ilev),ilev=1,nlev)   ! partial slant columns
         read(lunr_jac,*) ps_atm                    ! surface pressure (atm)
         read(lunr_jac,*) (z(ilev),ilev=1,nlev)     ! altitudes of levels
         read(lunr_jac,*) (pres(ilev),ilev=1,nlev)  ! pressures of levels

c  Augment A and B matrices with NFP extra rows, even if we're not going
c  to use them (ap_flag=0) to maintain sync with the .jac file.
         read(lunr_jac,*) (a(nmp+j,j),j=1,nfp)      ! ynoise/apu(j)
         read(lunr_jac,*) (b(nmp+j,1),j=1,nfp)      ! (ax(j)-cx(j))*ynoise/apu(j)

c  Copy (ax(j)-cx(j))*ynoise/apu(j) from b(nmp+j,1) to b(nmp+j,ilev),j=1,nfp, ilev=1,nlev
         do ilev=2,nlev
            call vmov(b(nmp+1,1),1,b(nmp+1,ilev),1,nfp)
         end do

c  Sum partial slant columns (psc) to yield total slant column (tsc)
         call vdot(psc,1,1.0,0,tsc,nlev)          ! total slant column
c
c  Write out the PD's in a form that can be easily plotted (e.g. xyplot)
c  This is for trouble-shooting/illustrative purposes only.
c        if (index(winfo,' wtf ').gt.0) then
c           open(lunw_wtf,file=filename(:lf)//'.wtf', status='unknown')
c           write(lunw_wtf,*)2,1+nfp+nlev
c           write(lunw_wtf,'(a4,999(9x,a1,i2.2))')
c    &      ' i  ',
c    &      ('C',j,j=1,nfp-ntg),
c    &      ('T',itg,itg=1,ntg),
c    &      ('S',ilev,ilev=0,nlev-1)
c           do imp=1,nmp
c              write(lunw_wtf,'(i5,999(1pe12.4))') imp,
c    &         (a(imp,itg),itg=ntg+1,nfp),  ! CL, CT, CC, FS, ZO
c    &         (a(imp,itg),itg=1,ntg),
c    &         (b(imp,ilev),ilev=1,nlev)
c           end do
c           close(lunw_wtf)
c        endif

c  Solve equation A.x=b  for multiple RHS, representing different levels.
c  A contains the Jacobians for the target gas columns, CL, FS, ZO, etc
c  b contains the single-level Jacobians of the first target gas.
c  A is MMPxNFP, x is NFPxNLEV, b is MMPxNLEV
c  x is returned inside b, the original contents of which are destroyed.
c  A is also destroyed.
         call shfti(a,mmp+mfp,nmp+nfp*ap_flag,nfp,b,mmp+mfp,nlev,tau,
     &   krank,rnorm,work,ip)
         if(krank.lt.nfp) write(*,*)' Rank Deficiency: ',krank,nfp
c
c  Write out the Averaging Kernels.
c  Also, calculate the pressure- and density-weighted kernels.
         open(lunw_aks,file=akpath(:lak)//'_'//specname(:ls),
     &   status='unknown')
         write(lunw_aks,*)3,3
         write(lunw_aks,*) version
         write(lunw_aks,*)  'Altitude_(km)   AK   Pressure_(atm)'
         do ilev=1,nlev
            if(psc(ilev).gt.0.0) exit
         end do
         ak1=b(1,ilev)*tsc/psc(ilev)
         ak2=b(1,ilev+1)*tsc/psc(ilev+1)
         akwas=(ak1*(ps_atm-pres(ilev+1))+ak2*(pres(ilev)-ps_atm))/
     &   (pres(ilev)-pres(ilev+1))
         pwas=ps_atm
         tak=0.0
         tb=0.0
         do jlev=ilev,nlev
            tb=tb+b(1,jlev)
            ak=b(1,jlev)*tsc/psc(jlev)
            write(lunw_aks,'(f6.2,f9.4,1pe10.3)') z(jlev),ak,pres(jlev)
            write(lunw_akall,*) ispec,zmin,sza,airmass,z(jlev),ak,
     &      pres(jlev)
            tak=tak+0.5*(ak+akwas)*(pwas-pres(jlev))
            akwas=ak
            pwas=pres(jlev)
         end do
         close(lunw_aks)
c         write(*,'(a,f9.6)')filename(:lf),tak/pres(1)
         write(*,'(a24,4f15.5)')specname(:ls),ps_atm,tak/ps_atm,tb
      end do ! ispec=1,999999
      write(*,*) 'Warning: Loop limit exceeded'
99    close(lunr_jac)
      close(lunw_akall)
      write(*,*)'Number of Measured Points   =',nmp
      write(*,*)'Number of Target Gases      =',ntg
      write(*,*)'Number of Fitted Parameters =',nfp
      write(*,*)'Number of Model Levels      =',nlev
      stop
      end
