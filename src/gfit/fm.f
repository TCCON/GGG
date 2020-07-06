      subroutine fm(lun_ak,slit,nhw,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & ldec,spts,spxv,
     & vac,splos,nlev,ncp,rdec,
     & cont_level, cont_tilt, cont_curv,xzo,
     & cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,nmp)

c  Forward Model: Computes a calculated spectrum and its matrix of partial differentials
c
c Inputs:
c   LUN_AK              I*4  Logical Unit Number for writing to AK file
c   SLIT(NII)           R*4  Pre-computed ILS (oversampled by LDEC wrt SPVAC)
c   NII                 I*4  Size of ILS vector
c   IFCSP               I*4  Index of First Calculated Spectral Point
c   IFMSP               I*4  Index of First Measured Spectral Point
c   IPTG                I*4  Pointer to first target gas
C   IPCL                I*4  Pointer to Continuum Level in the state vector
C   IPFS                I*4  Pointer to Frequency Shift in the state vector
C   IPZO                I*4  Pointer to Zero Offset     in the state vector
C   IPCF                I*4  Pointer Channel Fringe in the state vector (placeholder)
c   LDEC                I*4  Over-sampling factor for ILS
c   SPTS(NCP)           R*4  Pre-computed Solar Pseudo-Transmittance Spectrum
c   SPXV(NCP,0:NTG+2)   R*4  Pre-computed A Priori limb opacities (slant-path x VAC)
c   VAC(NCP,NLEV,0:NTG) R*4  Volume Absorption Coefficients (pre-computed)
c   SPLOS(NLEV)         R*4  LOS Slant Path distances
c   NLEV                I*4  Number of levels
c   NCP                 I*4  Number of primitive spectral grid points
c   RDEC                R*8  Ratio of primitive/observed spectral grids
c   SHSHS               R*8  Spectral shift expressed as fracion of GINT
c   CONT_LEVEL          R*4  Continuum Level (in case not being retrieved)
c   CONT_TILT           R*4  Continuum Tilt (in case not being retrieved)
c   CONT_CURV           R*4  Continuum  Curvature (in case not....)
c   XZO                 R*4  Zero Level Offset (in case not....)
c   CX(nfp)             R*4  forward model state vector 
c   NTG                 I*4  Number of Target Gases
c   NFP                 I*4  Number of Fitted Parameters (=NTG+5)
c   NMP                 I*4  Number of Measured Points
c
c Outputs:
c   CONT(nmp)           R*4  Continuum  CONT = CL.[1+CT.R+CC.S+...]
c   CALCUL(nmp)         R*4  Calculated spectrum = CONT.[SPTS*TRAN*(1-ZO)+ZO]
c   SLPD(NMP,NLEV)      R*4  Single Level Partial Differentials
c   PD(NFP,NMP+NFP)     R*4  matrix of partial differentials dCALC/dCX
c
c Description:
c  The forward model is described by the function
c     F(vi,X) = SLIT*{ CONT(vj,X).[TRAN(vj,X).SPTS(vj) + X(ipzo)] }
c  where "*" denotes convolution  and "." denotes multiplication
c  ipzo is a pointer to the Zero Offset element of the state vector
c  ipfs is a pointer to the FS Stretch element of state vector
c  ipsg is a pointer to the SG Stretch element of state vector
c  SLIT is the ILS
c  CONT is the continuum
c  TRAN is the atmospheric transmittance
c  SPTS is the Solar Pseudo-Transmittance Spectrum
c
c  If structure in CONT is much broader than that in SLIT, then
c     CALC(i) = CONT(i,CX).SLIT*[SPTS(j).TRAN(j,CX) + CX(ipzo)]
c
c     CONT(i,CX) = SUM_m[CX(ntg+m).DLPBF(i,m)]  m=1,NCT
c  is the continuum level
c  DLPBF are the Discrete Legendre Polynomial Basis Functions.
c  CX(ntg+m) (m=1,NCT) are the DLPBF continuum coefficients.
c  i is an index over the frequency grid of the measured spectrum
c  m is an index over DLP basis functions (NCT of them)
c
c     TRAN(j,CX) = Exp{ - SUM_g [CX(g).SPXV(g,j)] }  g=1,NTG 
c  is the monochromatic atmospheric transmittance
c  j is an index over the primitive/monochromatic frequency grid
c  g is an index over target gas
c
c     SPXV(g,j) = Sum_l [VAC(g,j,l).SP(l)]
c  is the limb opacity spectra of the different gases.
c  l is an index over the atmospheric levels.
c  SP(l) are the slant path distances for each level
c  VAC are the Volume Absorption Coefficients.
c
c  CX is the state vector with the following contents:
c  CX(j) (j=1,ntg) is the VMR Scaling Factor applied to the j'th target gas
c  CX(ipcl)  = CL  is the Continuum Level
c  CX(ipcl+1)  = CT  is the Continuum Tilt
c  CX(ipcl+2)  = CC  is the Continuum Curvature
c  CX(ntg+m) = CC  is the m'th DLP coefficient
c  CX(ipfs)  = FS  is the FS Stretch
c  CX(ipsg)  = SG  is the SG Stretch
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
c        dF/dSG     = PD(*,ipsg) = CL.W.d{Z}/di 
c        dF/dZO     = PD(*,ipzo) = CL 
c        dF/dCX(j)  = PD(*,j)    = CL.Z.V(i,j)
c
c where
c   Z = SPTS    ! Solar Transmittance Spectrum
c   W = TRAN   ! Telluric Transmittance
c
c    spxv(*,0)      ! Non-target gases
c    spxv(*,1)      ! First target gases
c    spxv(*,2)      ! Second target gases
c    spxv(*,ntg)    ! Last target gases
c    spxv(*,ntg+1)  ! Workspace for total transmittances
c    spxv(*,ntg+2)  ! Workspace for solar PD's
c    spxv(*,ntg+3)  ! Workspace for re-sampled solar spectrum
c
c  VAC, SPLOS are needed here only to computed single-level PD's
c  used in the computation of averaging kernels and in GFIT2.

      implicit none

      real*4 zero
      parameter(zero=0.0)

      include "ggg_int_params.f"
      include "const_params.f"
      
      integer ncp,jcp,nmp,ntg,jtg,nfp,nhw,ldec,kk,
c     & jj,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,ncbf,jbf,
     & lun_ak,jmp,nlev,ilev,nexpl,nexpr,noff,idum

      real*8 rdec,xfs,xsg,xx,fx,rdum
      real*4 slit(1+2*nhw*ldec),cx(nfp),vac(ncp,nlev,0:ntg),splos(nlev),
     & spts(ncp),
     & slpd(nmp,nlev),cont_level,cont_tilt,cont_curv,xzo,
     & cont(nmp),tran(nmp),calcul(nmp),
     & pd(nmp+nfp,nfp),spxv(ncp,0:ntg+3)

      rdum=big       ! Prevent compiler warnings (unused variable)
      idum=mauxcol   ! Prevent compiler warnings (unused variable)
      idum=mcolvav   ! Prevent compiler warnings (unused variable)
      idum=mfilepath ! Prevent compiler warnings (unused variable)
      idum=mgas      ! Prevent compiler warnings (unused variable)
      idum=mlev      ! Prevent compiler warnings (unused variable)
      idum=mrow_qc   ! Prevent compiler warnings (unused variable)
      idum=mspeci    ! Prevent compiler warnings (unused variable)
      idum=mvmode    ! Prevent compiler warnings (unused variable)
      idum=nchar     ! Prevent compiler warnings (unused variable)
      idum=ipcf      ! Prevent compiler warnings (unused variable)
      rdum=tiny      ! Prevent compiler warnings (unused variable)
      idum=iptg      ! Prevent compiler warnings (unused variable


      if(ncbf.ge.1) cont_level=cx(ipcl)
      if(ncbf.ge.2) cont_tilt=cx(ipcl+1)
      if(ncbf.ge.3) cont_curv=cx(ipcl+2)
      if(ipzo.gt.0) xzo=cx(ipzo)
      if(ipfs.gt.0) then
         xfs=cx(ipfs)
      else
         xfs=0.0d0
      endif

c      write(*,*)'fm.f: beginning ',ipfs,xfs,rdec
c  Compute primitive transmittance spectrum.
c  Scale the limb opacities by CX and co-add to produce the total limb opacity
c      write(*,*)'fm: spxv(1,0)=',spxv(1,0)
      call vmov(spxv(1,0),1,spxv(1,ntg+1),1,ncp)    ! non-target limb opacity
c      write(*,*)'fm: spxv(1,ntg+1)=',spxv(1,ntg+1)
      do jtg=1,ntg         ! compute  SUM_g {CX(g).SPXV(j,g)}
         call vsma(spxv(1,jtg),1,cx(jtg),spxv(1,ntg+1),1,
     &  spxv(1,ntg+1),1,ncp)
c         write(*,*)'fm loop:',jtg,spxv(1,jtg),spxv(1,ntg+1)
      end do
c
      call vexp(spxv(1,ntg+1),1,spxv(1,ntg+1),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}

cc  Code to model earthshine for GOSAT, by reducing gas
cc  absorptions by a factor ES.
c      es=0.000  ! No earthshine
c      es=0.0067 ! GOSAT band 3 (5000 cm-1)
c      es=0.0140 ! GOSAT band 2 (6000 cm-1)
c      es=0.0105 ! GOSAT band 1 (O2 A-band)
c      es=1.000  ! Normal situation (SPXV=SPXV)
c      do jj=1,ncp
c         spxv(jj,ntg+1)=1.0-es*(1-spxv(jj,ntg+1)) 
c      end do

c SPXV(*,ntg+1) contains gas transmittances on primative grid.
      
c Compute SG shift Jacobians and convert to stretches.
      if(ipsg.gt.0) then
         xsg=cx(ipsg)
         call lagrange_differentiate(ncp,spts,spxv(1,ntg+2))    ! dS/dx
c SPXV(*,ntg+2) contains dS/dx per primative grid point spacing.
         xx=xsg*ifcsp
         do jcp=1,ncp
            noff=nint(xx)
            if(jcp+noff.lt.1) then
               noff=1-jcp
            elseif(jcp+noff.gt.ncp) then
               noff=ncp-jcp
            endif
            fx=xx-noff
c  Stretched SPTS
            spxv(jcp,ntg+3)=spts(jcp+noff)+sngl(fx)*spxv(jcp+noff,ntg+2)
            xx=xx+xsg
         end do
c SPXV(*,ntg+3) contains SG-stretched SPTS on primative grid.

c SPXV(*,ntg+2) contains dS/dx on primative grid.
         call vmul(spxv(1,ntg+2),1,spxv(1,ntg+1),1,spxv(1,ntg+2),1,ncp)
         call regrid2(ifcsp,ncp,spxv(1,ntg+2),nhw,slit,ldec,
     &   rdec*(1.0d0+xfs),ifmsp,nmp,pd(1,ipsg),nexpl,nexpr)
         if(nexpl.ne.0)write(*,*)'Warning: FM: NEXPL=',nexpl
         if(nexpr.ne.0)write(*,*)'Warning: FM: NEXPR=',nexpr
c  Convert shift per measured spectral point to stretch
         do jmp=1,nmp
            pd(jmp,ipsg)=sngl(rdec)*cont_level*pd(jmp,ipsg)*
     &      float(ifmsp+jmp)
         end do
      else
         call vmov(spts,1,spxv(1,ntg+3),1,ncp)
      endif

c  Multiply resampled solar and atmospheric transmittance spectra
c  storing result in SPXV(*,NTG+1), overwriting the transmittance
      call vmul(spxv(1,ntg+3),1,spxv(1,ntg+1),1,spxv(1,ntg+1),1,ncp)   ! T.SPTS
c SPXV(*,ntg+1) contains TRAN.SPTS on primative grid.

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
c      write(*,*) 'slit=',ldec,(slit(k),k=1,1+2*ldec*nhw)
c
      
c      write(*,*) rdec, cx(ipfs)

c  Convolve product of high-res solar and telluric transmittance spectra
c  with ILS, at same time stretching by 1+xfs. Store result in TRAN.
      call regrid2(ifcsp,ncp,spxv(1,ntg+1),nhw,slit,ldec,
     & rdec*(1.d0+xfs),ifmsp,nmp,tran,nexpl,nexpr)
      if(nexpl.ne.0) write(*,*)'Warning: FM: NEXPL=',nexpl
      if(nexpr.ne.0) write(*,*)'Warning: FM: NEXPR=',nexpr
 
c  Compute continuum basis functions and store temporarily in PD(*,NTG+JBF)
c      write(*,*)'nmp,nfp,ncbf=',nmp,nfp,ncbf
      call vmov(cont_level,0,cont,1,nmp)
      if(ncbf.gt.0) call compute_dlpbf(nmp+nfp,nmp,ncbf,pd(1,ntg+1))
c  Compute CONT and associated PD's using new scheme.
c      write(*,*)'fm: cx = ',cx
      do jbf=2,ncbf
c         write(*,*)'fm: jfb,pd = ',jbf,pd(1,ntg+jbf),pd(2,ntg+jbf)
         call vsma(pd(1,ntg+jbf),1,cx(ntg+jbf)*cont_level,cont,1,
     &   cont,1,nmp)
         call vmul(pd(1,ntg+jbf),1,tran,1,pd(1,ntg+jbf),1,nmp)
         call vmul(pd(1,ntg+jbf),1,cont_level,0,pd(1,ntg+jbf),1,nmp)
      end do
c      write(*,*)'fm: cont=', cont(1),cont(2),cont(nmp-1),cont(nmp)
c
c  Apply zero offset (XZO) to TRANS 
c    CALCUL = TRAN*(1-XZO)+XZO
c           = TRAN + XZO*(1-TRAN)
      call vsub(1.0,0,tran,1,calcul,1,nmp)  !  CALCUL = 1-TRAN
c      write(*,*)'fm: calcul1=',
c     & calcul(1),calcul(2),calcul(nmp-1),calcul(nmp)
      if (ipzo.gt.0) call vmul(calcul,1,cont_level,0,pd(1,ipzo),1,nmp)
      call vsma(calcul,1,xzo,tran,1,calcul,1,nmp) ! TRAN + XZO*(1-TRAN)

c      write(*,*)'fm: calcul2=',
c     & calcul(1),calcul(2),calcul(nmp-1),calcul(nmp)
c  Multiply XZO-adjusted transmittance by CONT
c    CALCUL = CONT*[TRAN*(1-XZO)+XZO]
      call vmul(calcul,1,cont,1,calcul,1,nmp)  ! calcul = (ILS*T) . cont
c      write(*,*)'fm: calcul3=',
c     & calcul(1),calcul(2),calcul(nmp-1),calcul(nmp)

c  Compute CL PD
      if (ncbf.gt.0) call vmov(calcul(1),1,pd(1,ipcl),1,nmp)   ! pd(1) = calcul

c  Compute ZO PD
c      if (ipzo.gt.0) call vmov(cont_level,0,pd(1,ipzo),1,nmp)  ! ipzo: ZOFF PD's

c  Compute target gas PD's
      if(lun_ak.gt.1) write(lun_ak,*) nmp,ntg,nfp
      do jtg=1,ntg
         call vmul(spxv(1,ntg+1),1,spxv(1,jtg),1,spxv(1,ntg+2),1,ncp)
         call regrid2(ifcsp,ncp,spxv(1,ntg+2),nhw,slit,ldec,
     &   rdec*(1.0d0+xfs),ifmsp,nmp,pd(1,jtg),nexpl,nexpr)
         if(nexpl.ne.0) write(*,*)'Warning: FM: NEXPL=',nexpl
         if(nexpr.ne.0) write(*,*)'Warning: FM: NEXPR=',nexpr
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
     &         spxv(1,ntg+2),1,ncp)
               call regrid2(ifcsp,ncp,spxv(1,ntg+2),nhw,slit,ldec,
     &         rdec*(1.0d0+xfs),ifmsp,nmp,slpd(1,ilev),nexpl,nexpr)
               if(nexpl.ne.0)write(*,*)'Warning: FM: NEXPL=',nexpl
               if(nexpr.ne.0)write(*,*)'Warning: FM: NEXPR=',nexpr
               call vmul(slpd(1,ilev),1,cont,1,slpd(1,ilev),1,nmp)
               call vmul(slpd(1,ilev),1,ckm2cm*splos(ilev),0,
     &         slpd(1,ilev),1,nmp)
               if(lun_ak.gt.1) write(lun_ak,*)(slpd(jmp,ilev),jmp=1,nmp)
            end do
         end do
      endif   ! lun_ak.gt.0

c  Compute FS PDs, then convert to Stretch PDs
      if(ipfs.gt.0) then
         call lagrange_differentiate(nmp,calcul,pd(1,ipfs))
         do jmp=1,nmp
            pd(jmp,ipfs)=pd(jmp,ipfs)*float(ifmsp+jmp)
         end do
      endif

c  Zero the last NFP elements of each column of the PD Array (a priori).
c  If the parameter is not fitted, zero out the whole column.
      do kk=ntg+1,nfp
         call vmov(zero,0,pd(nmp+1,kk),1,nfp)       ! Zero last NFP elements
         if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,kk),jmp=1,nmp)   ! CL, CT, CC, FS, ZO Jacobians
      end do
      return
      end
