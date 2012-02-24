      subroutine fm2(lun_ak,winfo,slit,nii,ldec,spts,spxv,
     & vac,splos,nlev,ncp,rdec,sssss,cx,ntg,f,slpd,pd,nmp)
c
c  Computes the forward model
c    F(i,x)=X1.[1+X2.R+X4.SIN(CF)].SLIT(X3)*T(k)
c  and its matrix of partial differentials dF/dX,
c  where CF=X5.R+X6,  R=-0.5+(i-1)/(nmp-1), and T=exp[SUM{VAC(k,j).X(j)}]
c
c Inputs:
c   LUN_AK      I*4   Logical Unit Number
c   CX(ntg+4)           forward model state vector CX(ntg+3),
c   SPXV(NCP,NTG+3)    slant-path x absorption coefficient dot products
c   SPTS(NCP )         Solar Pseudo-Transmittance Spectrum
c   SLIT(NII)          pre-computed ILS which is oversampled by LDEC wrt SPVAC
c
c Outputs:
c      F(X)                 = X1.{[1+X2.R+X5.sin(CF)].SLIT(X3)*T + X4}
c     dF/dXj =  PD(i,j)     = X1.[1+X2.R+X5.sin(CF)].SLIT(X3)*(VAC(k,j).T)
c  X1.dF/dX1 =  PD(i,ntg+1) = X1.[1+X2.R+X5.sin(CF)].SLIT(X3)*T
c     dF/dX2 =  PD(i,ntg+2) = X1.R.SLIT(X3)*T
c     dF/dX3 =  PD(i,ntg+3) = X1.d{[1+X2.R+X5.sin(CF)].SLIT(X3)*T}/di 
c     dF/dX4 =  PD(i,ntg+4) = X1  ! Zero Offset
c
c Organization of X and SPXV arrays
c     cx(1)    ! First Target Gas
c     cx(2)    ! Second Target Gas
c     cx(ntg)  ! Last Target Gas
c     cx(n1)   ! Continuum Level (CL)
c     cx(n2)   ! Continuum Tilt (CT)
c     cx(n3)   ! Frequency Shift (FS)
c     cx(n4)   ! Zero Offset (ZO)
c
c     spxv(*,0)   ! Non-target gases
c     spxv(*,1)   ! First target gases
c     spxv(*,2)   ! Second target gases
c     spxv(*,ntg) ! Last target gases
c     spxv(*,n1)  ! Workspace
c     spxv(*,n2)  ! Workspace
      implicit none
      include "../ggg_const_params.f"
      include "const_params.f"

      integer ncp,nmp,ntg,jtg,nii,ldec,nterm,k,kk,jj,n1,n2,n3,n4,
     & nlev,ilev,jsp,jva,jmp,lun_ak
      real*4 slit(nii),cx(ntg+4),vac(ncp,nlev,0:ntg),splos(nlev),
     & rk,f(nmp),pd(nmp+ntg+4,ntg+4),spxv(ncp,0:ntg+2),
     & slpd(nmp,nlev,ntg),
     & spts(ncp)
      real*8 rdec,sh,sssss
      character winfo*(*),ss(4)*4
      data ss/' cl ',' ct ',' fs ',' zo '/
c
      n1=ntg+1  !
      n2=ntg+2  !
      n3=ntg+3  !
      n4=ntg+4  !

c  Multiple VAC by SPLOS
      jva=1
      jsp=1
      do jtg=0,ntg
         call vmov(zero,0,spxv(1,jtg),1,ncp)
         do ilev=1,nlev
            call vsma(vac(1,ilev,jtg),1,ckm2cm*splos(ilev),
     &      spxv(1,jtg),1,spxv(1,jtg),1,ncp)
            jva=jva+ncp
         end do
         jsp=jsp+ncp
      end do

c      write(*,*)'fm: ',x(n3),sssss
      sh=rdec*(cx(n3)+sssss)
c  Compute primitive transmittance spectrum using spxv(1,n1) as work space.
      call vmov(spxv(1,0),1,spxv(1,n1),1,ncp)    ! non-target VACs
      do jtg=1,ntg         ! compute  SUM{VAC(k,j).X(j)}
      call vsma(spxv(1,jtg),1,cx(jtg),spxv(1,n1),1,spxv(1,n1),1,ncp)
      end do
      call vexp(spxv(1,n1),1,spxv(1,n1),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}
      call vmul(spts,1,spxv(1,n1),1,spxv(1,n1),1,ncp)    ! STS*T
c
c  Compute PD's and convolve with ILS using PD(1,NTG+3) for workspace
      call newdec(spxv(1,n1),ncp,slit,nii,ldec,rdec,sh,pd(1,n1),nmp) 
      call vramp(pd(1,n2),1,nmp)                                 ! Ramp
      call vsma(pd(1,n2),1,cx(n1)*cx(n2),cx(n1),0,pd(1,n3),1,nmp)   ! a.(1+b.R) 
      call vmul(pd(1,n2),1,pd(1,n1),1,pd(1,n2),1,nmp)            ! R.SLIT*T
      call vmul(pd(1,n2),1,cx(n1),0,pd(1,n2),1,nmp)               ! a.R.SLIT*T
      call vmul(pd(1,n1),1,pd(1,n3),1,pd(1,n1),1,nmp)     ! a.(1+b.R).SLIT*T
      call vadd(pd(1,n1),1,cx(n1)*cx(n4),0,f(1),1,nmp)      ! f = pd(1)+zoff
      call vmov(f(1),1,pd(1,n1),1,nmp)                      ! pd(1) = f
      call vmov(cx(n1),0,pd(1,n4),1,nmp)  !  Zero-Offset PD's
      
c  Compute target gas PD's
      if(lun_ak.gt.0) write(lun_ak,*) nmp,ntg
      do jtg=1,ntg
         call vmul(spxv(1,n1),1,spxv(1,jtg),1,spxv(1,n2),1,ncp)
         call newdec(spxv(1,n2),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp)
         call vmul(pd(1,jtg),1,pd(1,n3),1,pd(1,jtg),1,nmp)
         call vmov(zero,0,pd(nmp+1,jtg),1,ntg+4)  ! zero unused part of PD array
      if(lun_ak.gt.0) write(lun_ak,*) (pd(jmp,jtg),jmp=1,nmp)
      end do
c
c  Compute single level partial differentials = v . df/dv
c      write(*,*)'fm2: lun_ak, spxv=',lun_ak,spxv(1,n1),spxv(2,n1)
      if(lun_ak.gt.0) write(lun_ak,*)nmp,nlev-1
      do ilev=2,nlev
      do jtg=1,min0(ntg,1)   !  only do the first target gas (if there is one)
        call vmul(spxv(1,n1),1,vac(1,ilev,jtg),1,spxv(1,n2),1,ncp)
        call newdec(spxv(1,n2),ncp,slit,nii,ldec,rdec,sh,
     &  slpd(1,ilev,jtg),nmp)
        call vmul(slpd(1,ilev,jtg),1,pd(1,n3),1,slpd(1,ilev,jtg),1,nmp)
        call vmul(slpd(1,ilev,jtg),1,ckm2cm*splos(ilev),0,
     &  slpd(1,ilev,jtg),1,nmp)
        call vmov(zero,0,pd(nmp+1,jtg),1,ntg+3)  ! zero unused part of PD array
        if(lun_ak.gt.0) write(lun_ak,*) (slpd(jmp,ilev,jtg),jmp=1,nmp)! Single-level PD's
      end do
      end do

      call vmov(zero,0,pd(nmp+1,n1),1,ntg+3) ! zero unused parts of PD array
      call vmov(zero,0,pd(nmp+1,n2),1,ntg+3) ! zero unused parts of PD array

c      if(index(winfo,ss(3)).gt.0) then  !  compute "fs" PD's.
      call vmov(zero,0,pd(1,n3),1,nmp+ntg+4)
      nterm=min0(4,nmp/2-1)
      do k=1,nterm     !  Apply triangular apodization to sinx/x operator
         jj=nmp-2*k
         rk=(nterm-k+1)*(-1)**k/float(k)/nterm
         call vsma(pd(1,n1),1,rk,pd(1+k,n3),1,pd(1+k,n3),1,jj)
         call vsma(pd(1+2*k,n1),1,-rk,pd(1+k,n3),1,pd(1+k,n3),1,jj)
      end do
c      endif
c
c  Write CT, CT, FS  partial differentials
      if(lun_ak.gt.0) then
         write(lun_ak,*) (pd(jmp,n1),jmp=1,nmp) ! CL PD
         write(lun_ak,*) (pd(jmp,n2),jmp=1,nmp) ! CT PD
         write(lun_ak,*) (pd(jmp,n3),jmp=1,nmp) ! FS PD
      endif
c
c  Zero out the last NTG+4 elements of each column of the PD Array (a priori).
c  If the gas is not fitted, zero out the whole column.
      do kk=1,4
        if(index(winfo,ss(kk)).gt.0) then
          call vmov(zero,0,pd(nmp+1,ntg+kk),1,ntg+4) ! Zero last NTG+4 elements
        else
          call vmov(zero,0,pd(1,ntg+kk),1,nmp+ntg+4) ! Zero whole column
        endif
      end do
      return
      end
