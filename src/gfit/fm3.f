      subroutine fm3(lun_ak,winfo,slit,nii,ldec,spts,spxv,
     & z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,vac,
     & splos,nlev,ncp,rdec,sssss,cx,ntg,nfp,calc,pd,tcalc,tpd,nmp)
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
      include "../ggg_const_params.f"
      include "const_params.f"
      
      integer ncp,nmp,ntg,jtg,nfp,nii,ldec,nterm,k,kk,jj,
     & n1,n2,n3,n4,n5,
     & rc,lun_ak,nlev,ilev,jva,jsp,nfov,ifov,kfov
      real*4 slit(nii),cx(nfp),vac(ncp,nlev,0:ntg),splos(nlev),
     & z(nlev),t(nlev),p(nlev),
     & spts(ncp),
     & rk,calc(nmp),pd(nmp+nfp,nfp),
     & tcalc(nmp),tpd(nmp+nfp,nfp),
     & spxv(ncp,0:ntg+2),
     & solzen,fovr,roc,fbar,zmin,bend,bend0,frangl,wt,twt,sza

      real*8 rdec,sh,sssss,fovo,wavtkr,obalt
      character winfo*(*),ss(4)*4
      data ss/' cl ',' ct ',' fs ',' zo '/

c      do jtg=0,ntg+2
c      call vdot(spxv(1,jtg),1,spxv(1,jtg),1,sum2,ncp)
c      write(*,*)'FM: sum SPXV:',jtg,sqrt(sum2/ncp)
c      end do

      n1=ntg+1  ! 
      n2=ntg+2  ! 
      n3=ntg+3  !
      n4=ntg+4  !
      n5=ntg+5  !
      sh=rdec*(cx(n3)+sssss)

      call vmov(zero,0,tcalc,1,nmp)
      do jtg=1,nfp
      call vmov(zero,0,tpd(1,jtg),1,nmp+nfp)
      end do

      fovr=90.*sngl(fovo)/spi ! convert radians diameter to deg radius
      twt=0.0
      nfov=1
      kfov=index(winfo,'nfov=')
      if(kfov.gt.0) read(winfo(kfov+5:),'(i1)')nfov
      call tlpath(nlev-1,z(2),t(2),p(2),solzen,0.0,roc,sngl(obalt),
     $sngl(wavtkr),fbar,zmin,bend0,splos(2),rc)
      do ifov=1,nfov
       frangl=float(2*ifov-nfov-1)/nfov  ! varies between +/- (1-1/nfov)
       sza=-solzen+bend0+fovr*frangl
       wt=sqrt(1.0-frangl**2)
c      write(*,*)'Calling TLPATH....',solzen,solwas,fovr,roc,obalt,wavtkr
      call tlpath(nlev-1,z(2),t(2),p(2),sza,0.0,roc,sngl(obalt),
     $sngl(wavtkr),fbar,zmin,bend,splos(2),rc)
c      write(*,*)ifov,wt,zmin,bend
      if(rc.ne.0) write(6,*)'Error in SLPATH:',rc
c      write(33,*)solzen,fovr,roc,sngl(obalt),zmin,bend

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

c  Compute primitive transmittance spectrum using spxv(1,n3) as work space.
c  Scale the limb opacities by CX and co-add to produce the total limb opacity
      call vmov(spxv(1,0),1,spxv(1,n1),1,ncp)    ! non-target limb opacity
      do jtg=1,ntg         ! compute  SUM{VAC(k,j).X(j)}
      call vsma(spxv(1,jtg),1,cx(jtg),spxv(1,n1),1,spxv(1,n1),1,ncp)
      end do
c
      call vexp(spxv(1,n1),1,spxv(1,n1),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}
      call vmul(spts,1,spxv(1,n1),1,spxv(1,n1),1,ncp)    ! STS*T
c
c  Compute PD's and convolve with ILS using PD(1,NTG+3) for workspace
      call newdec(spxv(1,n1),ncp,slit,nii,ldec,rdec,sh,pd(1,n1),nmp)  ! SLIT*T
      call vramp(pd(1,n2),1,nmp)                                    ! Ramp
      call vsma(pd(1,n2),1,cx(n1)*cx(n2),cx(n1),0,pd(1,n3),1,nmp)   ! a.(1+b.R) 
      call vmul(pd(1,n2),1,pd(1,n1),1,pd(1,n2),1,nmp)               ! R.SLIT*T
      call vmul(pd(1,n2),1,cx(n1),0,pd(1,n2),1,nmp)               ! a.R.SLIT*T
      call vmul(pd(1,n1),1,pd(1,n3),1,pd(1,n1),1,nmp)       ! a.(1+b.R).SLIT*T
      call vadd(pd(1,n1),1,cx(n1)*cx(n4),0,calc(1),1,nmp)   ! f = pd(1)+zoff
      call vmov(calc(1),1,pd(1,n1),1,nmp)                   ! pd(1) = f
      call vmov(cx(n1),0,pd(1,n4),1,nmp)  !  Zero-Offset PD's
      
c  Compute target gas PD's
      do jtg=1,ntg
         call vmul(spxv(1,n1),1,spxv(1,jtg),1,spxv(1,n2),1,ncp)
         call newdec(spxv(1,n2),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp)
         call vmul(pd(1,jtg),1,pd(1,n3),1,pd(1,jtg),1,nmp)
         call vmov(zero,0,pd(nmp+1,jtg),1,nfp)  ! zero unused part of PD array
      end do
c
c      if(index(winfo,ss(3)).gt.0) then  !  compute "fs" PD's.
         call vmov(zero,0,pd(1,n3),1,nmp+nfp)
         nterm=min0(4,nmp/2-1)
         do k=1,nterm     !  Apply triangular apodization to sinx/x operator
            jj=nmp-2*k
            rk=(nterm-k+1)*(-1)**k/float(k)/nterm
            call vsma(pd(1,n1),1,rk,pd(1+k,n3),1,pd(1+k,n3),1,jj)
            call vsma(pd(1+2*k,n1),1,-rk,pd(1+k,n3),1,pd(1+k,n3),1,jj)
         end do
c      endif
c
      twt=twt+wt
      call vsma(calc,1,wt,tcalc,1,tcalc,1,nmp)
      do jtg=1,nfp
         call vsma(pd(1,jtg),1,wt,tpd(1,jtg),1,tpd(1,jtg),1,nmp)
      end do

      end do  !  ifov=-nfov,nfov
c
      call vdiv(tcalc,1,twt,0,calc,1,nmp)
      do jtg=1,nfp
      call vdiv(tpd(1,jtg),1,twt,0,pd(1,jtg),1,nmp)
      end do

c  Zero out the last NTG+4 elements of each column of the PD Array (a priori).
c  If the gas is not fitted, zero out the whole column.
      do kk=1,4
        if(index(winfo,ss(kk)).gt.0) then
          call vmov(zero,0,pd(nmp+1,ntg+kk),1,nfp) ! Zero last NTG+4 elements
        else
          call vmov(zero,0,pd(1,ntg+kk),1,nmp+nfp) ! Zero whole column
        endif
      end do
      return
      end
