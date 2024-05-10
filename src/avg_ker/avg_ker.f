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
c  quantities as a function of frequency.  The NFP retrieved
c  parameters include the VSF of the various target gases and
c  CL, CT, FS, SG, ZO, (but not channel fringes).
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
      
      integer*4 lunr_pd,lunr_col,lunw_wtf,lunw_aks,lunw_akall,
     & mmp,nmp,imp,mfp,j,ap_flag,
     & ntg,itg,nlev,ilev,jlev,nfp,fbc,fnbc,lnbc,lf,i,ispec,nlhead,idum
      parameter (lunr_pd=12,lunr_col=13,lunw_wtf=14,lunw_aks=15,
     & lunw_akall=16,mmp=32000,mfp=16)
      integer*4 krank,ip(mfp),nit
      real*4 a(mmp+mfp,mfp),b(mmp+mfp,mlev),work(mfp),rnorm(mlev),
     & tau,psc(mlev),z(mlev),pres(mlev),ps_atm,pwas,tsc,
     & cl,tilt,cc,fqshift,
     & sg,zlo,rmsoclpc,zmin,airmass,sza,
     & ak1,ak2,ak,akwas,tak,tb
      character colfile*128,filename*128,akpath*128,spectrum*80,
     & version*60,winfo*128,col_fmt*128

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol ! Avoid compiler warning (unused parameter)
      idum=mcolvav ! Avoid compiler warning (unused parameter)
      idum=mgas    ! Avoid compiler warning (unused parameter)
      idum=mlev    ! Avoid compiler warning (unused parameter)
      idum=mrow_qc ! Avoid compiler warning (unused parameter)
      idum=mspeci  ! Avoid compiler warning (unused parameter)
      idum=mvmode  ! Avoid compiler warning (unused parameter)
      idum=ncell   ! Avoid compiler warning (unused parameter)
      idum=nchar   ! Avoid compiler warning (unused parameter)

      version=' avg_ker   Version 1.41     2018-12-20   GCT'
      write(*,*) version
      tau=1.e-7
      ap_flag=1   ! Include a priori constraints
      ap_flag=0   ! Ignore a priori constraints

c  Open .col file and read header.
      if (iargc() == 0) then
         write(*,*)'Enter name of relevent .col file'
         read(*,'(a)') colfile
      elseif (iargc() == 1) then
         call getarg(1, colfile)
      else
         stop 'Usage: $gggpath/bin/avg_ker colfile'
      endif
      open(lunr_col,file=colfile,status='old')
      read(lunr_col,*) nlhead
      do i=2,nlhead-5
         read(lunr_col,'(34x,a)') akpath
      end do
      read(lunr_col,*)
      read(lunr_col,*)
      read(lunr_col,'(34x,a)') col_fmt
      read(lunr_col,'(a)') winfo
      read(lunr_col,*) 

      write(*,*)'   Spectrum Path/Name              '//
     &'         Ps     P-averaged_AK   C-averaged_AK'

      open(lunw_akall,file=colfile(:lnbc(colfile)-3)//'akall',
     & status='unknown')
      write(lunw_akall,*) 2,7
      write(lunw_akall,*)' ispec zmin airmass sza  z ak p'

c  Main loop over spectra.
      do ispec=1,99999
         read(lunr_col,col_fmt,end=99) spectrum,nit,cl,tilt,cc,fqshift,
     &    sg,zlo,rmsoclpc,zmin,airmass
          sza=acos(1.03/airmass-0.03)*180/3.14159
         filename=akpath(:lnbc(akpath))//'_'//spectrum(fnbc(spectrum):)
         lf=fbc(filename)-1
c
c  Read total column partial differentials of first target gas.
         open (lunr_pd, file=filename(:lf), status='old')
         read(lunr_pd,*)nmp,ntg,nfp

c  Check that hard-wired array bounds won't be exceeded.
         if(nmp.gt.mmp) then
            write(*,*)'NMP, MMP = ',nmp,mmp
            stop 'increase parameter MMP'
         endif

         if(nfp.gt.mfp) then
            write(*,*)'NFP, MFP = ',nfp,mfp
            stop 'increase parameter MFP'
         endif

         do itg=1,ntg
            read(lunr_pd,*) (a(imp,itg),imp=1,nmp)
         end do
c
c  Read single-level partial differentials of first target gas.
         read(lunr_pd,*)nmp,nlev
         if(nmp.gt.mmp) stop 'increase MMP'
         if(nlev.gt.mlev) stop 'increase MLEV'
         do ilev=1,nlev
            read(lunr_pd,*) (b(imp,ilev),imp=1,nmp) ! Single-Level PDs of First Target gas
            call vmov(0.0,0,b(nmp+1,ilev),1,nfp)
         end do
c
c  Read continuum (and FS, SG, ZO) PDs
         do j=ntg+1,nfp
            read(lunr_pd,*) (a(imp,j),imp=1,nmp)
            call vmov(0.0,0,a(nmp+1,j),1,nfp)
         enddo

c  Read other parameters, including a priori constraint
         read(lunr_pd,*) (psc(ilev),ilev=1,nlev)   ! partial slant columns
         read(lunr_pd,*) ps_atm                    ! surface pressure (atm)
         read(lunr_pd,*) (z(ilev),ilev=1,nlev)     ! altitudes of levels
         read(lunr_pd,*) (pres(ilev),ilev=1,nlev)  ! pressures of levels
         if(ap_flag.ge.1) then  ! Augment A and B matrices with NFP extra rows
            read(lunr_pd,*) (a(nmp+j,j),j=1,nfp)      ! ynoise/apu(j)
            read(lunr_pd,*) (b(nmp+j,1),j=1,nfp)      ! (ax(j)-cx(j))*ynoise/apu(j)

c  Copy (ax(j)-cx(j))*ynoise/apu(j) from b(nmp+j,1) to b(nmp+j,ilev),j=1,nfp, ilev=1,nlev
            do ilev=2,nlev
               call vmov(b(nmp+1,1),1,b(nmp+1,ilev),1,nfp)
            end do
         endif
         close(lunr_pd)

c  Sum partial slant columns (psc) to yield total slant column (tsc)
         call vdot(psc,1,1.0,0,tsc,nlev)          ! total slant column
c
c  Write out the PD's in a form that can be easily plotted (e.g. xyplot)
c  This is for trouble-shooting/illustrative purposes only.
         if (index(winfo,' wtf ').gt.0) then
            open(lunw_wtf,file=filename(:lf)//'.wtf', status='unknown')
            write(lunw_wtf,*)2,1+nfp+nlev
            write(lunw_wtf,'(a4,999(9x,a1,i2.2))')
     &      ' i  ',
     &      ('C',j,j=1,nfp-ntg),
     &      ('T',itg,itg=1,ntg),
     &      ('S',ilev,ilev=0,nlev-1)
            do imp=1,nmp
               write(lunw_wtf,'(i5,999(1pe12.4))') imp,
     &         (a(imp,itg),itg=ntg+1,nfp),  ! CL, CT, CC, FS, ZO
     &         (a(imp,itg),itg=1,ntg),
     &         (b(imp,ilev),ilev=1,nlev)
            end do
            close(lunw_wtf)
         endif

c  Solve equation A.x=b  for multiple RHS, representing different levels.
c  A contains the Jacobians for the target gas columns, CL, FS, ZO, etc
c  b contains the single-level Jacobians of the first target gas.
c  A is MMPxNFP, x is NFPxNLEV, b is MMPxNLEV
         call shfti(a,mmp+mfp,nmp+nfp*ap_flag,nfp,b,mmp+mfp,nlev,tau,
     &   krank,rnorm,work,ip)
         if(krank.lt.nfp) write(*,*)' Rank Deficiency: ',krank,nfp
c
c  Write out the Averaging Kernels.
c  Also, calculate the pressure- and density-weighted kernels.
         open(lunw_aks,file=filename(:lf)//'.aks', status='unknown')
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
            write(lunw_aks,'(f5.1,f9.4,1pe10.3)') z(jlev),ak,pres(jlev)
            write(lunw_akall,*) ispec,zmin,airmass,sza,z(jlev),ak,
     &      pres(jlev)
            tak=tak+0.5*(ak+akwas)*(pwas-pres(jlev))
            akwas=ak
            pwas=pres(jlev)
         end do
         close(lunw_aks)
c         write(*,'(a,f9.6)')filename(:lf),tak/pres(1)
         write(*,'(a,4f15.5)')filename(:lf),ps_atm,tak/ps_atm,tb
      end do ! ispec=1,999999
      write(*,*) 'Warning: Loop limit exceeded'
99    close(lunr_col)
      close(lunw_akall)
      write(*,*)'Number of Measured Points =',nmp
      write(*,*)'Number of Target Gases    =',ntg
      write(*,*)'Number of Model Levels    =',nlev
      stop
      end
