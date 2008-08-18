      subroutine ak(luna,slit,nii,ldec,spxv,splos,vac,nlev,ncp,rdec,
     & x,ntg,pd,nmp)
c
c  Computes PD's for individual levels for use in calculation
c  of the averaging kernels.
c
c Inputs:
c   LUNA             I*4   Logical unit number for output file
c   SLIT(NII)        R*4   Pre-computed ILS
c   NII              I*4   Length of SLIT vector
c   LDEC             I*4   Oversampling factor of SLIT vector WRT VAC
c   SPXV(NCP,NTG+3)  R*4   Workspace array for SP(L)*VAC(L,K) & building PDs
c   SPLOS(NLEV)      R*4   Slant Path vector
c   VAC(NCP,NLEV,*)  R*4   Precomputed absorption coeficient array
c   NLEL             I*4   Number of atmospheric levels.
c   NCP              I*4   Number of Calculated spectral Points
c   RDEC             R*8   Measurement point spacing / primative point spacing
c   X(ntg+3)         R*4   forward model state vector X(ntg+3),
c   NNT              I*4   Group number for Non-Target gases
c   NTG              I*4   Number of Target Gases that were fitted
c
c Outputs:
c   PD(NMP.NTG)      R*4   Matrix of Partial Differentials
c     dF/dXj =  PD(i,j)     = X1.[1+X2.R+X4.sin(CF)].SLIT(X3)*(VAC(k,j).T)
c  X1.dF/dX1 =  PD(i,ntg+1) = X1.[1+X2.R+X4.sin(CF)].SLIT(X3)*T
c     dF/dX2 =  PD(i,ntg+2) = X1.R.SLIT(X3)*T
c     dF/dX3 =  PD(i,ntg+3) = X1.d{[1+X2.R+X4.sin(CF)].SLIT(X3)*T}/di
c
      implicit none
      integer*4 ncp,nmp,jmp,ntg,jtg,nii,ldec,nterm,k,
     & jj,n1,n2,n3,n4,nlev,jlev,luna
      real*4 slit(nii),x(ntg+3),vac(ncp,nlev,0:*),splos(nlev),
     & zero,unity,rk,pd(nmp+ntg+3,ntg+3),spxv(ncp,ntg+4),ckm2cm
      real*8 rdec,sh
      parameter (zero=0.0,unity=1.0,ckm2cm=1.0E+05)
c
      n1=ntg+1   ! Start address of Non-Target VACs
      n2=ntg+2   ! Start address of Solar Transmittance Spectrum
      n3=ntg+3   ! Workspace
      n4=ntg+4   ! Workspace
      sh=rdec*x(n3)
c
c  Compute primitive transmittance spectrum using spxv(1,n3) as work space.
      call vmov(spxv(1,n1),1,spxv(1,n3),1,ncp)    ! non-target VACs
      do jtg=1,ntg         ! compute  SUM{VAC(k,j).X(j)}
         call vsma(spxv(1,1+jtg),1,x(jtg),spxv(1,n3),1,spxv(1,n3),1,ncp)
      end do
      call vexp(spxv(1,n3),1,spxv(1,n3),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}
      call vmul(spxv(1,n2),1,spxv(1,n3),1,spxv(1,n3),1,ncp)    ! STS*T
      write(*,*)'ak : spxv=',spxv(1,n3),spxv(2,n3)
c
c  Compute PD's and convolve with ILS using PD(1,NTG+3) for workspace
c      call newdec(spxv(1,n3),rdec,slit,ldec,sh,nii,pd(1,n1),nmp) ! SLIT*T
      call newdec(spxv(1,n3),ncp,slit,nii,ldec,rdec,sh,pd(1,n1),nmp) 
      call vramp(pd(1,n2),1,nmp)                                 ! Ramp
      call vsma(pd(1,n2),1,x(n1)*x(n2),x(n1),0,pd(1,n3),1,nmp)   ! a.(1+b.R) 
      call vmul(pd(1,n2),1,pd(1,n1),1,pd(1,n2),1,nmp)            ! R.SLIT*T
      call vmul(pd(1,n2),1,x(n1),0,pd(1,n2),1,nmp)               ! a.R.SLIT*T
      call vmul(pd(1,n1),1,pd(1,n3),1,pd(1,n1),1,nmp)     ! a.(1+b.R).SLIT*T
c
c  Compute and write total column PD's
      write(luna,*) nmp,ntg
      do jtg=1,ntg
      call vmul(spxv(1,n3),1,spxv(1,1+jtg),1,spxv(1,n4),1,ncp)
      call newdec(spxv(1,n4),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp) 
      call vmul(pd(1,jtg),1,pd(1,n3),1,pd(1,jtg),1,nmp)
      call vmov(zero,0,pd(nmp+1,jtg),1,ntg+4)  ! zero unused part of PD array
      write(luna,*) (pd(jmp,jtg),jmp=1,nmp)  ! total column PD's
      end do
c
c  Write single level partial differentials = v . df/dv
      write(luna,*)nmp,nlev-1
      do jlev=2,nlev
      do jtg=1,1   !  only do the first target gas
      call vmul(spxv(1,n3),1,vac(1,jlev,jtg),1,spxv(1,n4),1,ncp)
      call newdec(spxv(1,n4),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp) 
      call vmul(pd(1,jtg),1,pd(1,n3),1,pd(1,jtg),1,nmp)
      call vmov(zero,0,pd(nmp+1,jtg),1,ntg+3)  ! zero unused part of PD array
      write(luna,*) (splos(jlev)*ckm2cm*pd(jmp,jtg),jmp=1,nmp)  ! Single-level PD's
      end do
      end do
c
      call vmov(zero,0,pd(nmp+1,n1),1,ntg+3) ! zero unused parts of PD array
      call vmov(zero,0,pd(nmp+1,n2),1,ntg+3) ! zero unused parts of PD array
c
      nterm=min0(4,nmp/2-1)
      call vmov(zero,0,pd(1,n3),1,nmp+ntg+3)
      do k=1,nterm     !  Apply triangular apodization to sinx/x operator
         jj=nmp-2*k
         rk=(nterm-k+1)*(-1)**k/float(k)/nterm
         call vsma(pd(1,n1),1,rk,pd(1+k,n3),1,pd(1+k,n3),1,jj)
         call vsma(pd(1+2*k,n1),1,-rk,pd(1+k,n3),1,pd(1+k,n3),1,jj)
      end do
      write(luna,*) (pd(jmp,n1),jmp=1,nmp) ! CL PD
      write(luna,*) (pd(jmp,n2),jmp=1,nmp) ! CT PD
      write(luna,*) (pd(jmp,n3),jmp=1,nmp) ! FS PD
c
      return
      end
