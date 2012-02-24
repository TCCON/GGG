      subroutine fm(lun_ak,winfo,slit,nii,ldec,spts,spxv,
     & vac,splos,nlev,ncp,rdec,shshs,cx,ntg,nfp,calc,pd,nmp)
c  Computes the forward model:
c    F(i,CX)=CX(n1).{ [1 + CX(n2).R + CX(n3).S].SLIT(CX(n4))*T(k) + CX(n5) }
c  and its matrix of partial differentials dF/dCX,
c  where 
c  n1=ntg+1
c  n2=ntg+2
c  n3=ntg+3
c  n4=ntg+4
c  n5=ntg+5
c  CX(j) (j=1,NTG) is the VMR Scaling Factor applied to the j'th target gas
c  CX(n1) is the Continuum Level (CL)
c  CX(n2) is the Continuum Tilt (CT)
c  CX(n3) is the Continuum Curvature (CC)
c  CX(n4) is the Frequency Shift (FS)
c  CX(n5) is the Zero Offset (ZO)
c
c  T(k)=exp[SUM_j {VAC(k,j).CX(j)}]
c  R = (i-ibar)/(nmp-1) = -0.5+(i-1)/(nmp-1)
c
c  S=a.R^2+b
c  a=3.(nmp-1)/(nmp-2)
c  b=-0.25.(nmp+1)/(nmp-2)
c
c  i is an index over the measured spectral points
c  k is an index over the primative spectral points
c  Let "." denote multiplication and "*" denote convolution
c  Let:  Z(i) = 1 + CX(n2).R(i) + CX(n3).S(i)
c  Let:  W(i) = SLIT(CX(n4))*T(k)
c  Let:  V(i,j) = SLIT(CX(n4))*(VAC(k,j).T(k)
c
c Inputs:
c   LUN_AK      I*4    Logical Unit Number
c   CX(nfp)     R*4    forward model state vector CX(ntg+3),
c   SPXV(NCP,NTG+3)    A Priori limb opacities (slant-path x VAC)
c   SPTS(NCP)          Solar Pseudo-Transmittance Spectrum
c   SLIT(NII)          pre-computed ILS (oversampled by LDEC wrt SPVAC)
c
c Outputs:
c        calc(CX)              = CX(n1).[Z.W+CX(n5)]
c CX(n1).dF/dCX(n1) = PD(*,n1) = CX(n1).[Z.W+CX(n5)]
c        dF/dCX(n2) = PD(*,n2) = CX(n1).R.W
c        dF/dCX(n3) = PD(*,n3) = CX(n1).S.W
c        dF/dCX(n4) = PD(*,n4) = CX(n1).d{Z.W}/di 
c        dF/dCX(n5) = PD(*,n5) = CX(n1) 
c        dF/dCX(j)  = PD(*,j)  = CX(n1).Z.V(i,j)
c
c  Note that CALC(CX) is identical to CX(n1).dF/dCX(n1) = PD(*,n1)
c
c Organization of X and SPXV arrays
c    cx(1)       ! oirst Target Gas
c    cx(2)       ! Second Target Gas
c    cx(i)       ! i'th target gas
c    cx(ntg)     ! Last Target Gas
c    cx(n1)   ! Continuum Level (CL)
c    cx(n2)   ! Continuum Tilt (CT)
c    cx(n3)   ! Continuum Curvature (CC)
c    cx(n4)   ! Frequency Shift (FS)
c    cx(n5)   ! Zero Offset (ZO)
c
c    spxv(*,0)   ! Non-target gases
c    spxv(*,1)   ! First target gases
c    spxv(*,2)   ! Second target gases
c    spxv(*,ntg) ! Last target gases
c    spxv(*,ntg+1)  ! Workspaca for total transmittances
c    spxv(*,ntg+2)  ! Workspace for PD's
c
      implicit none
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"
      include "const_params.f"
      
      integer ncp,nmp,ntg,jtg,nfp,nii,ldec,nterm,k,kk,jj,
     & n1,n2,n3,n4,n5,nnt,
c     & lunp,j,
     & lun_ak,jmp,nlev,ilev

      real*8 rdec,sh,shshs
      real*4 slit(nii),cx(nfp),vac(ncp,nlev,0:ntg),splos(nlev),
     & spts(ncp),sum2,
     & a,b,
     & rk,calc(nmp),pd(nmp+nfp,nfp),
     & spxv(ncp,0:ntg+2)

      parameter (nnt=5)
      character winfo*(*),ss(nnt)*4
c     & ,primsp*10

      data ss/' cl ',' ct ',' cc',' fs ',' zo '/
c      data lunp/39/

c      lunp=lunp+1
c      primsp='primsp.000'

      if(nfp-ntg.gt.nnt) stop ' Increase NNT:  nfp-ntg.gt.nnt '
c      nh=(1+ncp)/2
c      do jtg=0,ntg
c        call vdot(spxv(1,jtg),1,spxv(1,jtg),1,sum2,ncp)
c        write(*,*)'SPXV:',jtg,ntg,spxv(1,jtg),spxv(nh,jtg),spxv(ncp,jtg)
c      end do

      n1=ntg+1  ! 
      n2=ntg+2  ! 
      n3=ntg+3  !
      n4=ntg+4  !
      n5=ntg+5  !
      sh=rdec*(cx(n4)+shshs)
c      write(*,*)'fm.f: sh=',sh,rdec,cx(n4),shshs
c  Compute primitive transmittance spectrum using spxv(1,n4) as work space.
c  Scale the limb opacities by CX and co-add to produce the total limb opacity
      call vmov(spxv(1,0),1,spxv(1,n1),1,ncp)    ! non-target limb opacity
      do jtg=1,ntg         ! compute  SUM{VAC(k,j).X(j)}
      call vsma(spxv(1,jtg),1,cx(jtg),spxv(1,n1),1,spxv(1,n1),1,ncp)
      end do
c
      call vexp(spxv(1,n1),1,spxv(1,n1),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}
c      write(*,*)'CALC:',spxv(1,n1),spxv(nh,n1),spxv(ncp,n1)
      
c  Multiply solar and atmospheric transmittance spectra
      call vmul(spts,1,spxv(1,n1),1,spxv(1,n1),1,ncp)    ! STS*T

c  Write out primitive-grid transmittance
c       write(primsp(8:10),'(i3.3)') lunp-40
c       open(lunp,file=primsp,status='unknown')
c       do j=1,ncp
c          write(lunp,*) j,spxv(j,n1)
c       end do
c       close(lunp)


c      write(*,*)'SPTS:',spts(1),spts(nh),spts(ncp)
c      write(*,*)'CALC:',spxv(1,n1),spxv(nh,n1),spxv(ncp,n1)
c      call vdot(spxv(1,n1),1,spxv(1,n1),1,sum2,ncp)
c      write(*,*)'fm: ss=',sqrt(sum2/ncp)
c      write(*,*) 'slit=',ldec,(slit(k),k=1,nii)
c
c  Compute PD's and convolve with ILS. Can use PD(1,n4) for workspace
c  since the FS PD is computed later.
      call newdec(spxv(1,n1),ncp,slit,nii,ldec,rdec,sh,pd(1,n1),nmp)  ! n1: W=SLIT*T
      call vdot(pd(1,n1),1,pd(1,n1),1,sum2,nmp)
c      write(*,*) 'pd=',pd(1,n1),pd(2,n1),pd(3,n1),pd(nmp,n1)
c      write(*,*)'fm: sum=',sqrt(sum2/nmp)
      a=     3*float(nmp-1)/(nmp-2)
      b= -0.25*float(nmp+1)/(nmp-2)
      call vmov(zero,0,pd(1,n4),1,nmp)   !  Zero pd(*,n4) before using as workspace
      call vramp(pd(1,n2),1,nmp)                                   ! n2: Ramp = R
      call vmul(pd(1,n2),1,pd(1,n2),1,pd(1,n3),1,nmp)              ! n3: R^2
      call vsma(pd(1,n3),1,a,b,0,pd(1,n3),1,nmp)                   ! n3: S=a.R^2+b
      call vsma(pd(1,n2),1,cx(n1)*cx(n2),cx(n1),0,pd(1,n4),1,nmp)  ! n4: cl.(1+ct.R) 
      call vsma(pd(1,n3),1,cx(n1)*cx(n3),pd(1,n4),1,pd(1,n4),1,nmp)! n4: cl.[1+ct.R+cc.S]
      call vmul(pd(1,n3),1,cx(n1),0,pd(1,n3),1,nmp)                ! n3: Si = cx(n1)*S
      call vmul(pd(1,n2),1,pd(1,n1),1,pd(1,n2),1,nmp)              ! n2: R.W
      call vmul(pd(1,n2),1,cx(n1),0,pd(1,n2),1,nmp)                ! n2: cl.R.W
      call vmul(pd(1,n1),1,pd(1,n4),1,pd(1,n1),1,nmp)              ! n1: cl.[1+ct.R+cc.S].W
      call vadd(pd(1,n1),1,cx(n1)*cx(n5),0,calc(1),1,nmp)          ! f = pd(1)+zoff
      call vmov(calc(1),1,pd(1,n1),1,nmp)                          ! pd(1) = f
      call vmov(cx(n1),0,pd(1,n5),1,nmp)  !  Zero-Offset PD's
      
c
c      do i=1,nmp
c        r=-0.5+float(i-1)/(nmp-1)
c        s=a*r2**2+b
c        cont=cx(n1)*(1+cx(n2)*r+cx(n3)*s)
c        pd(i,n2)=cx(n1)*r2*pd(i,n1)
c        pd(i,n3)=cx(n1)*r3*pd(i,n1)
c        pd(i,n5)=cx(n1)
c      end do

c  Compute target gas PD's
      if(lun_ak.gt.0) write(lun_ak,*) nmp,ntg
      do jtg=1,ntg
         call vmul(spxv(1,n1),1,spxv(1,jtg),1,spxv(1,n2),1,ncp)
         call newdec(spxv(1,n2),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp)
         call vmul(pd(1,jtg),1,pd(1,n4),1,pd(1,jtg),1,nmp)
         call vmov(zero,0,pd(nmp+1,jtg),1,nfp)  ! zero unused part of PD array
         if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,jtg),jmp=1,nmp)
      end do

c  Calculate & Write single level partial differentials = v . df/dv
      if(lun_ak.gt.0) then
         write(lun_ak,*)nmp,nlev-ncell
         do ilev=ncell+1,nlev
         do jtg=1,min0(ntg,1)   !  only do the first target gas
         call vmul(spxv(1,n1),1,vac(1,ilev,jtg),1,spxv(1,n2),1,ncp)
         call newdec(spxv(1,n2),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp)
         call vmul(pd(1,jtg),1,pd(1,n4),1,pd(1,jtg),1,nmp)
         call vmov(zero,0,pd(nmp+1,jtg),1,ntg+3)  ! zero unused part of PD array
         write(lun_ak,*) (splos(ilev)*ckm2cm*pd(jmp,jtg),jmp=1,nmp)
         end do
         end do
      endif   ! lun_ak.gt.0

c
c      if(index(winfo,ss(3)).gt.0) then  !  compute "fs" PD's.
         call vmov(zero,0,pd(1,n4),1,nmp+nfp)
         nterm=min0(4,nmp/2-1)
         do k=1,nterm     !  Apply triangular apodization to sinx/x operator
            jj=nmp-2*k
            rk=(nterm-k+1)*(-1)**k/float(k)/nterm
            call vsma(pd(1,n1),1,rk,pd(1+k,n4),1,pd(1+k,n4),1,jj)
            call vsma(pd(1+2*k,n1),1,-rk,pd(1+k,n4),1,pd(1+k,n4),1,jj)
         end do
c      endif
c
cc  Write CT, CT, FS  partial differentials
c      if(lun_ak.gt.0) then
c         write(lun_ak,*) (pd(jmp,n1),jmp=1,nmp) ! CL PD
c         write(lun_ak,*) (pd(jmp,n2),jmp=1,nmp) ! CT PD
c         write(lun_ak,*) (pd(jmp,n4),jmp=1,nmp) ! FS PD
c         write(lun_ak,*) (pd(jmp,n5),jmp=1,nmp) ! ZO PD
c      endif
c
c  Zero the last NFP=NTG+5 elements of each column of the PD Array (a priori).
c  If the parameter is not fitted, zero out the whole column.
      do kk=1,nnt
c        write(*,*)'fm: ',kk,ss(kk),index(winfo,ss(kk))
        if(index(winfo,ss(kk)).gt.0) then
          call vmov(zero,0,pd(nmp+1,ntg+kk),1,nfp) ! Zero last NFP elements
        else
          call vmov(zero,0,pd(1,ntg+kk),1,nmp+nfp) ! Zero whole column
        endif
        if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,ntg+kk),jmp=1,nmp)   ! CL, CT, FS, ZO Jacobians
      end do
      return
      end
