      subroutine fm(lun_ak,slit,nii,
     & iptg,ipcl,ipfs,ipzo,ipcf,
     & ldec,spts,spxv,
     & vac,splos,nlev,ncp,rdec,shshs,
     & cont_level, cont_tilt, cont_curv,xzo,
     & cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,nmp)

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
c   SHSHS               R*8  Spectral shift expressed as fracion of GINT
c   CX(nfp)             R*4  forward model state vector 
c   NTG                 I*4  Number of Target Gases
c   NFP                 I*4  Number of Fitted Parameters (=NTG+5)
c   NMP                 I*4  Number of Measured Points
c
c Outputs:
c   CONT(nmp)           R*4  Continuum  = CL.[1+CT.R+CC.S+...]
c   CALC(nmp)           R*4  Calculated spectrum = CL.[TRANS+ZO]
c   SLPD(NMP,NLEV)      R*4  Single Level Partial Differentials
c   PD(NFP,NMP+NFP)     R*4  matrix of partial differentials dCALC/dCX
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
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"
      include "const_params.f"
      
      integer ncp,nmp,ntg,jtg,nfp,nii,ldec,nterm,k,kk,jj,
     & iptg,ipcl,ipfs,ipzo,ipcf,ncbf,jbf,
     & lun_ak,jmp,nlev,ilev

      real*8 rdec,sh,shshs
      real*4 slit(nii),cx(nfp),vac(ncp,nlev,0:ntg),splos(nlev),
     & spts(ncp),
     & slpd(nmp,nlev),cont_level,cont_tilt,cont_curv,xzo,
     & rk,cont(nmp),calcul(nmp),pd(nmp+nfp,nfp),spxv(ncp,0:ntg+2)

      if(ncbf.ge.1) cont_level=cx(ipcl)
      if(ncbf.ge.2) cont_tilt=cx(ipcl+1)
      if(ncbf.ge.3) cont_curv=cx(ipcl+2)
      if(ipzo.gt.0) xzo=cx(ipzo)
      if(ipfs.gt.0) then
         sh=rdec*(cx(ipfs)+shshs)
      else
         sh=rdec*shshs
      endif
c      write(*,*)'fm.f: sh=',sh,rdec,cx(ipfs),shshs
c  Compute primitive transmittance spectrum using CONT as work space.
c  Scale the limb opacities by CX and co-add to produce the total limb opacity
      call vmov(spxv(1,0),1,spxv(1,ntg+1),1,ncp)    ! non-target limb opacity
c      write(*,*)'fm:',0,spxv(1,0),spxv(1,ntg+1)
      do jtg=1,ntg         ! compute  SUM_g {CX(g).SPXV(j,g)}
         call vsma(spxv(1,jtg),1,cx(jtg),spxv(1,ntg+1),1,
     &   spxv(1,ntg+1),1,ncp)
c         write(*,*)'fm:',jtg,spxv(1,jtg),spxv(1,ntg+1)
      end do
c
      call vexp(spxv(1,ntg+1),1,spxv(1,ntg+1),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}
      
c  Multiply solar and atmospheric transmittance spectra
      call vmul(spts,1,spxv(1,ntg+1),1,spxv(1,ntg+1),1,ncp)    ! STS*T

c  Write out primitive-grid transmittance
c       write(primsp(8:10),'(i3.3)') lunp-40
c       open(lunp,file=primsp,status='unknown')
c       do j=1,ncp
c          write(lunp,*) j,spxv(j,ntg+1)
c       end do
c       close(lunp)


c      write(*,*)'SPTS:',spts(1),spts(nh),spts(ncp)
c      write(*,*)'CALC:',spxv(1,ntg+1),spxv(nh,ntg+1),spxv(ncp,ntg+1)
c      call vdot(spxv(1,ntg+1),1,spxv(1,ntg+1),1,sum2,ncp)
c      write(*,*)'fm: ss=',sqrt(sum2/ncp)
c      write(*,*) 'slit=',ldec,(slit(k),k=1,nii)
c
      
c  Convolve hi-res transmittance spectrum with ILS.
c  Store result temporarily in PD(1,1)
      call newdec(spxv(1,ntg+1),ncp,slit,nii,ldec,rdec,sh,
     & pd(1,1),nmp)  !  ILS*T
c
c  Compute basis functions and store temporarily in PD(*,NTG+JBF)
      call vmov(cont_level,0,cont,1,nmp)
      if(ncbf.gt.0) call compute_dlpbf(nmp+nfp,nmp,ncbf,pd(1,ntg+1))
c  Compute CONT and associated PD's using new scheme.
c  Use PD(*,NTG+JBF) as temporary storage for Basis Functions
      do jbf=2,ncbf
         call vsma(pd(1,ntg+jbf),1,cx(ntg+jbf)*cont_level,cont,1,
     &   cont,1,nmp)
         call vmul(pd(1,ntg+jbf),1,pd(1,1),1,pd(1,ntg+jbf),1,nmp)
         call vmul(pd(1,ntg+jbf),1,cont_level,0,pd(1,ntg+jbf),1,nmp)
      end do
      call vmul(pd(1,1),1,cont,1,calcul,1,nmp)  ! calcul = (ILS*T) . cont
      call vadd(calcul,1,cont_level*xzo,0,calcul,1,nmp)    ! calc = calc + zoff*CL
      if (ncbf.gt.0) call vmov(calcul(1),1,pd(1,ipcl),1,nmp)   ! pd(1) = calcul
      if (ipzo.gt.0) call vmov(cont_level,0,pd(1,ipzo),1,nmp)  ! ipzo: ZOFF PD's

c  Compute target gas PD's
      if(lun_ak.gt.1) write(lun_ak,*) nmp,ntg,nfp
      do jtg=1,ntg
         call vmul(spxv(1,ntg+1),1,spxv(1,jtg),1,spxv(1,ntg+2),1,ncp)
         call newdec(spxv(1,ntg+2),ncp,slit,nii,ldec,rdec,sh,
     &   pd(1,jtg),nmp)
         call vmul(pd(1,jtg),1,cont,1,pd(1,jtg),1,nmp)
         call vmov(zero,0,pd(nmp+1,jtg),1,nfp)           ! zero unused part of PD array
         if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,jtg),jmp=1,nmp)
      end do

c  Calculate & Write single level partial differentials = v . df/dv
      if(lun_ak.gt.0) then
         if(lun_ak.gt.1)write(lun_ak,*)nmp,nlev-ncell

         do ilev=ncell+1,nlev
         do jtg=1,min0(ntg,1)   !  only do the first target gas
            call vmul(spxv(1,ntg+1),1,vac(1,ilev,jtg),1,
     &      spxv(1,ntg+2),1,ncp)
            call newdec(spxv(1,ntg+2),ncp,slit,nii,ldec,rdec,sh,
     &      slpd(1,ilev),nmp)
            call vmul(slpd(1,ilev),1,cont,1,slpd(1,ilev),1,nmp)
            call vmul(slpd(1,ilev),1,ckm2cm*splos(ilev),0,
     &      slpd(1,ilev),1,nmp)
            if(lun_ak.gt.1) write(lun_ak,*) (slpd(jmp,ilev),jmp=1,nmp)
         end do
         end do
      endif   ! lun_ak.gt.0

c  Compute FS PDs
      if(ipfs.gt.0) call lagrange_differentiate(nmp,calcul,pd(1,ipfs))
c      if(ipfs.gt.0) then
c         call vmov(zero,0,pd(1,ipfs),1,nmp+nfp)
c         nterm=min0(4,nmp/2-1)
c         do k=1,nterm     !  Apply triangular apodization to sinx/x operator
c            jj=nmp-2*k
c            rk=(nterm-k+1)*(-1)**k/float(k)/nterm
c            call vsma(calcul(1),1,rk,pd(1+k,ipfs),1,pd(1+k,ipfs),1,jj)
c            call vsma(calcul(1+2*k),1,-rk,pd(1+k,ipfs),1,
c     &      pd(1+k,ipfs),1,jj)
c         end do
c      endif

c  Zero the last NFP elements of each column of the PD Array (a priori).
c  If the parameter is not fitted, zero out the whole column.
      do kk=ntg+1,nfp
         call vmov(zero,0,pd(nmp+1,kk),1,nfp)       ! Zero last NFP elements
         if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,kk),jmp=1,nmp)   ! CL, CT, CC, FS, ZO Jacobians
      end do
      return
      end

      subroutine lagrange_differentiate(np,vin,vout)
c  Computes the gradient of the function that goes through
c  the points vin(np) at the points themselves.
c
c  Inputs:
c            NMP   I*4  Number of points in in/output vectors
c        VIN(NMP)  R*4  Input vector
c
c  Outputs:
c       VOUT(NMP)  R*4  Output vector
c
c  Does this by fitting a Lagrange polynomial, differentiating, 
c  then evaluating the differential at the equally-spaced abscissae.
c  Uses 3-point Lagrange interpolation for the edge points
c  Uses 4-point Lagrange interpolation for the penultimate points
c  Uses 5-point Lagrange interpolation for the rest of the window
c  Although the 5-point operator only has 4 terms, this is because
c  the coefficient of the VIN(i) term is zero.  It is still based
c  on a 5-point Lagrange Interpolation.

      integer*4 np,i
      real*4 vin(np),vout(np)

      if(np.ge.4) then
        vout(1)=(-3*vin(1)+4*vin(2)-vin(3))/2                  ! 3-point
        vout(2)=(-2*vin(1)-3*vin(2)+6*vin(3)-vin(4))/6         ! 4-point
        do i=1+2,np-2
          vout(i)=(vin(i-2)-8*vin(i-1)+8*vin(i+1)-vin(i+2))/12 ! 5-point
        end do
        vout(np-1)=(2*vin(np)+3*vin(np-1)-6*vin(np-2)+vin(np-3))/6 ! 4-p
        vout(np)=(3*vin(np)-4*vin(np-1)+vin(np-2))/2           ! 3-point
      elseif(np.eq.3) then
        vout(1)=(-3*vin(1)+4*vin(2)-vin(3))/2
        vout(2)=(-vin(1)+vin(3))/2
        vout(3)=(3*vin(1)-4*vin(2)+vin(3))/2
      elseif(np.eq.2) then
        vout(1)=(-vin(1)+vin(2))/2
        vout(2)=(-vin(1)+vin(2))/2
      elseif(np.eq.1) then
        vout(1)=0.0
      endif
      return
      end
