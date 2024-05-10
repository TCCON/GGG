      subroutine fm3(lun_ak,slit,nhw,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & ldec,spts,spxv,
     & nfov,z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     & vac,splos,nlev,ncp,rdec,
     & cont_level,cont_tilt,cont_curv,xzo,
     & cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,tcalc,tpd,nmp)

c  Forward Model: Computes a calculated spectrum and its matrix of partial differentials
c
c Inputs:
c   LUN_AK              I*4  Logical Unit Number for writing to AK file
c   WINFO               C**  List of state vector variables to be adjusted.
c   SS
c   SLIT(NII)           R*4  Pre-computed ILS (oversampled by LDEC wrt SPVAC)
c   NII                 I*4  Size of ILS vector
c   LDEC                I*4  Over-sampling factor for ILS
c   SPTS(NCP)           R*4  Pre-computed Solar Pseudo-Transmittance Spectrum
c   SPXV(NCP,0:NTG+2)   R*4  Pre-computed A Priori limb opacities (slant-path x VAC)
C   IPCL                I*4  Pointer to Continuum Level in the state vector
C   IPCT                I*4  Pointer to Continuum Tilt  in the state vector
C   IPCC                I*4  Pointer to Continuum Curvature in the state vector
C   IPFS                I*4  Pointer to Frequency Shift in the state vector
C   IPZO                I*4  Pointer to Zero Offset     in the state vector
c   VAC(NCP,NLEV,0:NTG) R*4  Volume Absorption Coefficients (pre-computed)
c   SPLOS(NLEV)         R*4  LOS Slant Path distances
c   NLEV                I*4  Number of levels
c   NCP                 I*4  Number of primitive spectral grid points
c   RDEC                R*8  Ratio of primitive/observed spectral grids
c   SHSHS               R*8  Spectral shift
c   CX(nfp)             R*4  forward model state vector
c   NTG                 I*4  Number of Target Gases
c   NFP                 I*4  Number of Fitted Parameters (=NTG+5)
c   NMP                 I*4  Number of Measured Points
c
c Outputs:
c   CONT(nmp)           R*4  Continuum  = CL.[1+CT.R+CC.S]
c   CALC(nmp)           R*4  Calculated spectrum = CL.[Z.W+ZO]
c   SLPD(NMP,NLEV)      R*4  Single Level Partial Differentials
c   PD(NFP,NMP+NFP)     R*4  matrix of partial differentials dCALC/dCX
c
c
c Description:
c     CALC(i) = CONT(i,CX).{SLIT(CX(ipfs))*T(i,CX)} + CX(ipzo)
c  where "*" denotes convolution  and "." denotes multiplication
c  ipzo is a pointer to the Zero Offset element of the state vector
c  ipfs is a pointer to the Frequency Shift element of state vector
c
c     CONT(i,CX) = SUM_m[CX(m+ntg).DLPBF(i,m)]  m=1,NCT
c  is the continuum level
c  DLPBF are the Discrete Legendre Polynomial Basis Functions.
c  CX(ntg+m) (m=1,NCT) are the DLPBF continuum coefficients.
c  i is an index over the frequency grid of the measured spectrum
c  m is an index over DLP basis functions (NCT of them)
c
c     T(j,CX) = Exp{ - SUM_g [CX(g).SPXV(g,j)] }  g=1,NTG
c  is the monochromatic atmospheric transmittance
c  j is an index over the primitive/monochromatic frequency grid
c  g is an index over target gas
c
c     SPXV(g,j) = Sum_l [VAC(g,j,l).SP(l)]
c  is the limb opacity spectra of the different gases.
c  l is an index over the atmospheric levels.
c  SP(l) are the slant path distances for each level
c
c  CX is the state vector with the following contents:
c  CX(j) (j=1,ntg) is the VMR Scaling Factor applied to the j'th target gas
c  CX(ipcl)  = CL  is the Continuum Level
c  CX(ipcl+1)  = CT  is the Continuum Tilt
c  CX(ipcl+2)  = CC  is the Continuum Curvature
c  CX(ntg+m) = CC  is the m'th DLP coefficient
c  CX(ipfs)  = FS  is the Frequency Shift
c  CX(ipzo)  = ZO  is the Zero Offset
c
c  Partial Differentials (aka Jacobians)
c     CALC(i) = CONT(i,CX).{SLIT(CX(ipfs))*T(i,CX)} + CX(ipzo)
c
c    dCALC(i)/dCX(m) = DLPBF(i,m).{SLIT(CX(ipfs))*T(i,CX)}
c    dCALC(i)/dCX(g) = CONT(i,m).{SLIT(CX(ipfs))*dT(i,CX)/dCX(g))}
c
c     CL.dF/dCL     = PD(*,ipcl) = CL.[Z.W+ZO]
c        dF/dCT     = PD(*,ipcl+1) = CL.R.W
c        dF/dCC     = PD(*,ipcl+2) = CL.S.W
c        dF/dFS     = PD(*,ipfs) = CL.d{Z.W}/di
c        dF/dZO     = PD(*,ipzo) = CL
c        dF/dCX(j)  = PD(*,j)    = CL.Z.V(i,j)
c
c    spxv(*,0)      ! Non-target gases
c    spxv(*,1)      ! First target gases
c    spxv(*,2)      ! Second target gases
c    spxv(*,ntg)    ! Last target gases
c    spxv(*,ntg+1)  ! Workspace for total transmittances
c    spxv(*,ntg+2)  ! Workspace for PD's
c
c  VAC, SPLOS are needed here only to computed single-level PD's
c  used in the computation of averaging kernels and in GFIT2.

      implicit none

      real*4  zero,spi
      parameter(zero=0.0,spi=3.14159265)

      include "ggg_int_params.f"
      include "const_params.f"
      
      integer ncp,nmp,jmp,ntg,jtg,nfp,nhw,ldec,nterm,k,kk,jj,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,ncbf,jbf,
     & nexpl,nexpr,
     & rc,lun_ak,nlev,ilev,jva,jsp,nfov,ifov,idum
      real*8 rdec,sh,shshs,hh,fovo,wavtkr,obalt,rdum
      real*4 slit(1+2*ldec*nhw),cx(nfp),vac(ncp,nlev,0:ntg),splos(nlev),
     & z(nlev),t(nlev),p(nlev),
     & spts(ncp),
     & slpd(nmp,nlev),cont_level,cont_tilt,cont_curv,xzo,
     & rk,cont(nmp),calcul(nmp),pd(nmp+nfp,nfp),spxv(ncp,0:ntg+2),
     & tcalc(nmp),tpd(nmp+nfp,nfp),
     & solzen,fovr,roc,fbar,zmin,bend,bend0,frangl,wt,twt,sza

      idum=lun_ak    ! Prevent compiler warnings (unused variable)
      rdum=big       ! Prevent compiler warnings (unused variable)
      rdum=slpd(1,1) ! Prevent compiler warnings (unused variable)
      idum=mauxcol   ! Prevent compiler warnings (unused variable)
      idum=mcolvav   ! Prevent compiler warnings (unused variable)
      idum=mfilepath ! Prevent compiler warnings (unused variable)
      idum=mgas      ! Prevent compiler warnings (unused variable)
      idum=mlev      ! Prevent compiler warnings (unused variable)
      idum=mrow_qc   ! Prevent compiler warnings (unused variable)
      idum=mspeci    ! Prevent compiler warnings (unused variable)
      idum=mvmode    ! Prevent compiler warnings (unused variable)
      idum=ncell     ! Prevent compiler warnings (unused variable)
      idum=nchar     ! Prevent compiler warnings (unused variable)
      idum=iptg      ! Prevent compiler warnings (unused variable)
      idum=ipcf      ! Prevent compiler warnings (unused variable)
      idum=ipsg      ! Prevent compiler warnings (unused variable)
      rdum=tiny      ! Prevent compiler warnings (unused variable)

      shshs=0.d0
      if(ncbf.ge.1) cont_level=cx(ipcl)
      if(ncbf.ge.2) cont_tilt=cx(ipcl+1)
      if(ncbf.ge.3) cont_curv=cx(ipcl+2)
      if(ipzo.gt.0) xzo=cx(ipzo)
      if(ipfs.gt.0) then
         sh=rdec*(cx(ipfs)+shshs)
         hh=rdec*cx(ipfs)/(ifcsp+0.5*(ncp-1))
      else
         sh=rdec*shshs
         hh=0.0d0
      endif

c Zero TCALC & TPD arrary.
      call vmov(zero,0,tcalc,1,nmp)
      do jtg=1,nfp
         call vmov(zero,0,tpd(1,jtg),1,nmp+nfp)
      end do

      fovr=90.*sngl(fovo)/spi ! convert radians diameter to deg radius
      twt=0.0
      call tlpath(nlev-1,z(2),t(2),p(2),solzen,0.0,roc,sngl(obalt),
     $sngl(wavtkr),fbar,zmin,bend0,splos(2),rc)
      do ifov=1,nfov
         frangl=float(2*ifov-nfov-1)/nfov  ! varies between +/- (1-1/nfov)
         sza=-solzen+bend0+fovr*frangl
         wt=sqrt(1.0-frangl**2)
c         write(*,*)'Calling TLPATH....',solzen,solwas,fovr,roc,obalt,wavtkr
         call tlpath(nlev-1,z(2),t(2),p(2),sza,0.0,roc,sngl(obalt),
     $   sngl(wavtkr),fbar,zmin,bend,splos(2),rc)
c         write(*,*)ifov,wt,zmin,bend
         if(rc.ne.0) write(6,*)'Error in SLPATH:',rc
c         write(33,*)solzen,fovr,roc,sngl(obalt),zmin,bend

c  Multiple VAC by SPLOS
         jva=1
         jsp=1
         do jtg=0,ntg
            call vmov(zero,0,spxv(1,jtg),1,ncp)
            do ilev=1,nlev
               call vsma(vac(1,ilev,jtg),1,ckm2cm*splos(ilev),
     &         spxv(1,jtg),1,spxv(1,jtg),1,ncp)
               jva=jva+ncp
            end do
            jsp=jsp+ncp
         end do

c      write(*,*)'fm.f: sh=',sh,rdec,cx(ipfs),shshs
c  Compute primitive transmittance spectrum using CONT as work space.
c  Scale the limb opacities by CX and co-add to produce the total limb opacity
         call vmov(spxv(1,0),1,spxv(1,ntg+1),1,ncp)    ! non-target limb opacity
c         write(*,*)'fm:',0,spxv(1,0),spxv(1,ntg+1)
         do jtg=1,ntg         ! compute  SUM_g {CX(g).SPXV(j,g)}
            call vsma(spxv(1,jtg),1,cx(jtg),spxv(1,ntg+1),1,
     &      spxv(1,ntg+1),1,ncp)
c            write(*,*)'fm:',jtg,spxv(1,jtg),spxv(1,ntg+1)
         end do
c
         call vexp(spxv(1,ntg+1),1,spxv(1,ntg+1),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}

c  Multiply solar and atmospheric transmittance spectra
         call vmul(spts,1,spxv(1,ntg+1),1,spxv(1,ntg+1),1,ncp)    ! STS*T

c  Write out primitive-grid transmittance
c         write(primsp(8:10),'(i3.3)') lunp-40
c         open(lunp,file=primsp,status='unknown')
c         do j=1,ncp
c            write(lunp,*) j,spxv(j,ntg+1)
c         end do
c         close(lunp)


c         write(*,*)'SPTS:',spts(1),spts(nh),spts(ncp)
c         write(*,*)'CALC:',spxv(1,ntg+1),spxv(nh,ntg+1),spxv(ncp,ntg+1)
c         call vdot(spxv(1,ntg+1),1,spxv(1,ntg+1),1,sum2,ncp)
c         write(*,*)'fm: ss=',sqrt(sum2/ncp)
c         write(*,*) 'slit=',ldec,(slit(k),k=1,1+2*ldec*nhw)
c

c  Convolve hi-res spectrum with ILS. Store result temporarily in PD(1,1)
c         write(*,*)'fm3: calling regrid2',nhw
         call regrid2(ifcsp,ncp,spxv(1,ntg+1),nhw,slit,ldec,
     &   rdec*(1.0d0+hh),ifmsp,nmp,calcul,nexpl,nexpr)
         if(nexpl.ne.0) write(*,*) 'Warning: FM3: NEXPL==',nexpl
         if(nexpr.ne.0) write(*,*) 'Warning: FM3: NEXPR==',nexpr

c  Compute basis functions and store temporarily in PD(*,NTG+JBF)
         call vmov(cont_level,0,cont,1,nmp)
         if(ncbf.gt.0) call compute_dlpbf(nmp+nfp,nmp,ncbf,pd(1,ntg+1))

c  Compute CONT and associated PD's using new scheme.
c  Use PD(*,NTG+JBF) as temporary storage for Basis Functions
         do jbf=2,ncbf
            call vsma(pd(1,ntg+jbf),1,cx(ntg+jbf)*cont_level,cont,1,
     &      cont,1,nmp)
            call vmul(pd(1,ntg+jbf),1,calcul,1,pd(1,ntg+jbf),1,nmp)
            call vmul(pd(1,ntg+jbf),1,cont_level,0,pd(1,ntg+jbf),1,nmp)
         end do
         call vmul(calcul,1,cont,1,calcul,1,nmp)
         call vadd(calcul,1,cont_level*xzo,0,calcul,1,nmp)    ! calc = pd(1)+zoff
         if (ncbf.gt.0) call vmov(calcul(1),1,pd(1,ipcl),1,nmp)   ! pd(1) = calc
         if(ipzo.gt.0) call vmov(cont_level,0,pd(1,ipzo),1,nmp)   ! ipzo: ZOFF PD's

c  Compute target gas PD's
         if(lun_ak.gt.1) write(lun_ak,*) nmp,ntg,nfp
         do jtg=1,ntg
            call vmul(spxv(1,ntg+1),1,spxv(1,jtg),1,spxv(1,ntg+2),1,ncp)
            call regrid2(ifcsp,ncp,spxv(1,ntg+2),nhw,slit,ldec,
     &      rdec*(1.0d0+hh),ifmsp,nmp,pd(1,jtg),nexpl,nexpr)
            if(nexpl.ne.0) write(*,*)'Warning: REGRID: NEXPL=',nexpl
            if(nexpr.ne.0) write(*,*)'Warning: REGRID: NEXPR=',nexpr
            call vmul(pd(1,jtg),1,cont,1,pd(1,jtg),1,nmp)
            call vmov(zero,0,pd(nmp+1,jtg),1,nfp)           ! zero unused part of PD array
            if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,jtg),jmp=1,nmp)
         end do
c
c  Compute FS PDs
         if(ipfs.gt.0) then
            call vmov(zero,0,pd(1,ipfs),1,nmp+nfp)
            nterm=min0(4,nmp/2-1)
            do k=1,nterm     !  Apply triangular apodization to sinx/x operator
               jj=nmp-2*k
               rk=(nterm-k+1)*(-1)**k/float(k)/nterm
               call vsma(calcul(1),1,rk,pd(1+k,ipfs),1,
     &         pd(1+k,ipfs),1,jj)
               call vsma(calcul(1+2*k),1,-rk,pd(1+k,ipfs),1,
     &         pd(1+k,ipfs),1,jj)
            end do
         endif

         twt=twt+wt
         call vsma(calcul,1,wt,tcalc,1,tcalc,1,nmp)
         do jtg=1,nfp
            call vsma(pd(1,jtg),1,wt,tpd(1,jtg),1,tpd(1,jtg),1,nmp)
         end do

      end do  !  ifov=1,nfov
c
      call vdiv(tcalc,1,twt,0,calcul,1,nmp)
      do jtg=1,nfp
         call vdiv(tpd(1,jtg),1,twt,0,pd(1,jtg),1,nmp)
      end do

c  Zero the last NFP elements of each column of the PD Array (a priori).
c  If the parameter is not fitted, zero out the whole column.
      do kk=ntg+1,nfp
         call vmov(zero,0,pd(nmp+1,kk),1,nfp)       ! Zero last NFP elements
         if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,kk),jmp=1,nmp)   ! CL, CT, CC, FS, ZO Jacobians
      end do
c
      return
      end
