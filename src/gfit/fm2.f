      subroutine fm2(lun_ak,slit,nii,
     & iptg,ipcl,ipfs,ipzo,ipcf,
     & ldec,spts,spxv,
     & vac,splos,nlev,ncp,rdec,sssss,
     & cont_level, cont_tilt, cont_curv,xzo,
     & cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,nmp,
     & prf_flag,sf_iter)

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

      real*4  zero
      parameter(zero=0.0)

      include "const_params.f"
      include "ggg_int_params.f"

      integer ncp,nmp,ntg,jtg,nfp,nii,ldec,
     & iptg,ipcl,ipfs,ipzo,ipcf,ncbf,
     & nlev,ilev,jsp,jva,lun_ak,nnt,ll
      parameter (nnt=5) ! Number of Non_Target-Gas parameters
      real*4 slit(nii),cx(nfp),vac(ncp,nlev,0:ntg),splos(nlev),
     & cont(nmp),calcul(nmp),pd(nmp+nfp,nfp),spxv(ncp,0:ntg+2),
     & temp(nlev),
     & cont_level, cont_tilt, cont_curv,xzo,
     & slpd(nmp,nlev),sf_iter(nlev,ntg),
     & spts(ncp)
      real*8 rdec,sssss
      logical prf_flag(ntg)


c  Multiply VAC by SPLOS, storing the produxt in SPXV.
c  This has to be done every iteration for the full profile retrieval,
c  versus once per spectrum for the simple scaling retrieval.
      jva=1
      jsp=1
      do jtg=1,ntg
         call vmov(zero,0,spxv(1,jtg),1,ncp)
         do ilev=1,nlev
            if(prf_flag(jtg))then
               temp(ilev)=splos(ilev)*sf_iter(ilev,jtg)
            else
               temp(ilev)=splos(ilev)
            endif
            call vsma(vac(1,ilev,jtg),1,ckm2cm*temp(ilev),
     &      spxv(1,jtg),1,spxv(1,jtg),1,ncp)
            jva=jva+ncp
         end do
         jsp=jsp+ncp
      end do 

c  If LL>0 then FM calculates single-level PDs and returns them in SLPD
c  If LL>1  FM also writes out SLPD into the already-open LUN_AK 
      if(prf_flag(1)) then
        ll=1
      elseif(lun_ak.gt.0) then
        ll=lun_ak
      else
        ll=0
      endif

c  Compute forward model.
c     call fm(ll,winfo,slit,nii,ldec,spts,spxv,
c    & vac,splos,nlev,ncp,rdec,sssss,cx,ntg,nfp,cont,f,slpd,pd,nmp)

      call fm(ll,slit,nii,
     & iptg,ipcl,ipfs,ipzo,ipcf,
     & ldec,spts,spxv,
     & vac,splos,nlev,ncp,rdec,sssss,
     & cont_level, cont_tilt, cont_curv,xzo,
     & cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,nmp)

      return
      end
