      subroutine fm(winfo,slit,nii,ldec,spxv,ncp,rdec,x,ntg,f,pd,nmp)
c  Computes the forward model
c    F(i,x)=X1.[1+X2.R+X4.SIN(CF)].SLIT(X3)*T(k)
c  and its matrix of partial differentials dF/dX,
c  where CF=X5.R+X6,  R=-0.5+(i-1)/(nmp-1), and T=exp[SUM{VAC(k,j).X(j)}]
c
c Inputs:
c   X(ntg+4)           forward model state vector X(ntg+3),
c   SPXV(NCP,NTG+3)    slant-path x absorption coefficient dot products
c   SPXV(*,N2)         Solar Transmittance Spectrum
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
c The first column of SPXV contain the non-target VACs
c The next NTG columns of SPXV contain the target VACs
c The NTG+2'nd column contains the solar transmittance spectrum
c The NTG+3'nd and NTG+4'rd columns are workspace.
c
      implicit none
      integer ncp,nmp,ntg,jtg,nii,ldec,nterm,k,kk,jj,n1,n2,n3,n4
      real*4 slit(nii),x(ntg+4),sumr2,
     & zero,unity,rk,f(nmp),pd(nmp+ntg+4,ntg+4),spxv(ncp,ntg+4)
      real*8 rdec,sh
      character winfo*(*),ss(4)*4
      parameter (zero=0.0,unity=1.0)
      data ss/' cl ',' ct ',' fs ',' zo '/
c
      n1=ntg+1  ! Start address of Non-target VAC's
      n2=ntg+2  ! Start address of Solar Transmittance Spectrum
      n3=ntg+3  ! Workspace
      n4=ntg+4  ! Workspace
      sh=rdec*x(n3)
c  Compute primitive transmittance spectrum using spxv(1,n3) as work space.
      call vmov(spxv(1,1),1,spxv(1,n3),1,ncp)    ! non-target VACs
      do jtg=1,ntg         ! compute  SUM{VAC(k,j).X(j)}
c      call vdot(spxv(1,1+jtg),1,spxv(1,1+jtg),1,sumr2,ncp)
c      write(*,*)'FM: jtg,rms=',jtg,sumr2
      call vsma(spxv(1,1+jtg),1,x(jtg),spxv(1,n3),1,spxv(1,n3),1,ncp)
      end do
      call vexp(spxv(1,n3),1,spxv(1,n3),1,ncp) ! T=exp[SUM{VAC(k,j).X(j)}
      call vmul(spxv(1,n2),1,spxv(1,n3),1,spxv(1,n3),1,ncp)    ! STS*T
c
c  Compute PD's and convolve with ILS using PD(1,NTG+3) for workspace
      call newdec(spxv(1,n3),ncp,slit,nii,ldec,rdec,sh,pd(1,n1),nmp) 
      call vramp(pd(1,n2),1,nmp)                                 ! Ramp
      call vsma(pd(1,n2),1,x(n1)*x(n2),x(n1),0,pd(1,n3),1,nmp)   ! a.(1+b.R) 
      call vmul(pd(1,n2),1,pd(1,n1),1,pd(1,n2),1,nmp)            ! R.SLIT*T
      call vmul(pd(1,n2),1,x(n1),0,pd(1,n2),1,nmp)               ! a.R.SLIT*T
      call vmul(pd(1,n1),1,pd(1,n3),1,pd(1,n1),1,nmp)     ! a.(1+b.R).SLIT*T
      call vadd(pd(1,n1),1,x(n1)*x(n4),0,f(1),1,nmp)      ! f = pd(1)+zoff
      call vmov(x(n1),0,pd(1,n4),1,nmp)  !  Zero-Offset PD's
      
c  Compute target gas PD's
      do jtg=1,ntg
      call vmul(spxv(1,n3),1,spxv(1,1+jtg),1,spxv(1,n4),1,ncp)
      call newdec(spxv(1,n4),ncp,slit,nii,ldec,rdec,sh,pd(1,jtg),nmp) 
      call vmul(pd(1,jtg),1,pd(1,n3),1,pd(1,jtg),1,nmp)
      call vmov(zero,0,pd(nmp+1,jtg),1,ntg+4)  ! zero unused part of PD array
      end do
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
