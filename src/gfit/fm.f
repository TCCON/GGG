      subroutine fm(lun_ak,winfo,slit,nii,ldec,spts,spxv,dspdzxv,
     & vac,splos,nlev,ncp,rdec,sssss,cx,ntg,calc,pd,nmp)
c  Computes the forward model
c    F(i,x)=X1.[1+X2.R+X4.SIN(CF)].SLIT(X3)*T(k)
c  and its matrix of partial differentials dF/dX,
c  where CF=X5.R+X6,  R=-0.5+(i-1)/(nmp-1), and T=exp[SUM{VAC(k,j).X(j)}]
c
c Inputs:
c   LUN_AK      I*4   Logical Unit Number
c   CX(ntg+4)           forward model state vector CX(ntg+3),
c   SPXV(NCP,NTG+3)    A Priori limb opacities (slant-path x VAC dot productsi over levels)
c   SPTS(NCP)          Solar Pseudo-Transmittance Spectrum
c   SLIT(NII)          pre-computed ILS which is oversampled by LDEC wrt SPVAC
c
c Outputs:
c     calc(X)               = X1.{[1+X2.R+X5.sin(CF)].SLIT(X3)*T(Xj)+X4}
c     dF/dXj =  PD(i,j)     = X1.[1+X2.R+X5.sin(CF)].SLIT(X3)*(VAC(k,j).T(Xj))
c  X1.dF/dX1 =  PD(i,ntg+1) = X1.{[1+X2.R+X5.sin(CF)].SLIT(X3)*T(Xj)+X4}
c     dF/dX2 =  PD(i,ntg+2) = X1.R.SLIT(X3)*T(Xj)
c     dF/dX3 =  PD(i,ntg+3) = X1.d{[1+X2.R+X5.sin(CF)].SLIT(X3)*T(Xj)}/di 
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
c
      implicit none
      integer ncp,nmp,ntg,jtg,nii,ldec,nterm,k,kk,jj,n1,n2,n3,n4,
     & lun_ak,jmp,nlev,ilev,i
      real*4 slit(nii),cx(ntg+4),vac(ncp,nlev,0:ntg),splos(nlev),
     & spts(ncp),ckm2cm,sum2,
     & zero,unity,rk,calc(nmp),pd(nmp+ntg+4,ntg+4),
     & spxv(ncp,0:ntg+2),dspdzxv(ncp,0:ntg+2),d2,d4

      real*8 rdec,sh,sssss
      character winfo*(*),ss(4)*4
      parameter (zero=0.0,unity=1.0)
      parameter (ckm2cm=100000.0)
      data ss/' cl ',' ct ',' fs ',' zo '/

c      do jtg=0,ntg+2
c      call vdot(spxv(1,jtg),1,spxv(1,jtg),1,sum2,ncp)
c      write(*,*)'FM: sum SPXV:',jtg,sqrt(sum2/ncp)
c      call vdot(dspdzxv(1,jtg),1,dspdzxv(1,jtg),1,sum2,ncp)
c      write(*,*)'FM: sum dSPdzXV:',jtg,sqrt(sum2/ncp)
c      end do

      n1=ntg+1  ! 
      n2=ntg+2  ! 
      n3=ntg+3  !
      n4=ntg+4  !
      sh=rdec*(cx(n3)+sssss)
c      write(*,*)'fm.f: sh=',sh,rdec,cx(n3),sssss
c  Compute primitive transmittance spectrum using spxv(1,n3) as work space.
c  Scale the limb opacities by CX and co-add to produce the total limb opacity
      call vmov(spxv(1,0),1,spxv(1,n1),1,ncp)    ! non-target limb opacity
      call vmov(dspdzxv(1,0),1,dspdzxv(1,n1),1,ncp)    ! non-target limb opacity
      do jtg=1,ntg         ! compute  SUM{VAC(k,j).X(j)}
      call vsma(spxv(1,jtg),1,cx(jtg),spxv(1,n1),1,spxv(1,n1),1,ncp)
      call vsma(dspdzxv(1,jtg),1,cx(jtg),dspdzxv(1,n1),1,
     & dspdzxv(1,n1),1,ncp)
      end do
c
      call vexp(spxv(1,n1),1,spxv(1,n1),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}
c  Evaluate Bessel_I1(eta)
      do i=1,ncp
         d2=dspdzxv(i,n1)**2
         d4=d2*d2
c         dspdzxv(i,n1)=1+d2/8+d4/192+d2*d4/9216+d4*d4/737280
         dspdzxv(i,n1)=1+d2/8+d4/192+d2*d4/9216
      end do
      call vmul(dspdzxv(1,n1),1,spxv(1,n1),1,spxv(1,n1),1,ncp)  ! FOV correction
      call vmul(spts,1,spxv(1,n1),1,spxv(1,n1),1,ncp)    ! STS*T
      call vdot(spxv(1,n1),1,spxv(1,n1),1,sum2,ncp)
c      write(*,*)'fm: ss=',sqrt(sum2/ncp)
c      write(*,*) 'slit=',ldec,(slit(k),k=1,nii)
c
c  Compute PD's and convolve with ILS using PD(1,NTG+3) for workspace
      call newdec(spxv(1,n1),ncp,slit,nii,ldec,rdec,sh,pd(1,n1),nmp)  ! SLIT*T
      call vdot(pd(1,n1),1,pd(1,n1),1,sum2,nmp)
c      write(*,*) 'pd=',pd(1,n1),pd(2,n1),pd(3,n1),pd(nmp,n1)
c      write(*,*)'fm: sum=',sqrt(sum2/nmp)
      call vramp(pd(1,n2),1,nmp)                                    ! Ramp
      call vsma(pd(1,n2),1,cx(n1)*cx(n2),cx(n1),0,pd(1,n3),1,nmp)   ! a.(1+b.R) 
      call vmul(pd(1,n2),1,pd(1,n1),1,pd(1,n2),1,nmp)               ! R.SLIT*T
      call vmul(pd(1,n2),1,cx(n1),0,pd(1,n2),1,nmp)               ! a.R.SLIT*T
      call vmul(pd(1,n1),1,pd(1,n3),1,pd(1,n1),1,nmp)       ! a.(1+b.R).SLIT*T
      call vadd(pd(1,n1),1,cx(n1)*cx(n4),0,calc(1),1,nmp)   ! f = pd(1)+zoff
      call vmov(calc(1),1,pd(1,n1),1,nmp)                   ! pd(1) = f
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

c  Write single level partial differentials = v . df/dv
      if(lun_ak.gt.0) then
         write(lun_ak,*)nmp,nlev-1
         do ilev=2,nlev
         do jtg=1,min0(ntg,1)   !  only do the first target gas
         call vmul(spxv(1,n1),1,vac(1,ilev,jtg),1,spxv(1,n2),1,ncp)
         call newdec(spxv(1,n2),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp)
         call vmul(pd(1,jtg),1,pd(1,n3),1,pd(1,jtg),1,nmp)
         call vmov(zero,0,pd(nmp+1,jtg),1,ntg+3)  ! zero unused part of PD array
         write(lun_ak,*) (splos(ilev)*ckm2cm*pd(jmp,jtg),jmp=1,nmp)
         end do
         end do
      endif   ! lun_ak.gt.0

c
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
