      subroutine do_retrieval(obsrvd,nmp,apx,apu,slit,nii,
     & iptg,ipcl,ipfs,ipzo,ipcf,xzo,
     & cont_level,cont_tilt,cont_curv,
     & ldec,rdec,spts,spxv,vac,splos,nlev,ncp,ntg,ncbf,nfp,snr,
     & corrld,sssss,debug,mit,nit,cont,calcul,rmsocl,cx,ex,
     & slpd,pd,ssnmp)

c  Adjusts the state vector (cx) to get the best fit to measured spectrum (y).
c  Performs an iterative solution to the optimal estimation equation
c  (K'.S-1.K+Sa-1)dx = (K'.Sy-1.(y-f(x))+Sa-1.(x-xa))              
c  where K is the matrix of partial differentials and K' is its transpose
c
c  Assumes that the noise in the spectrum is white, so that Sy, the
c  measurement covariance, is a diagonal matrix of constant value.
c
c  Further assumes that the state vector a priori covariance (Sa) is diagonal.
c
c   Inputs:
c     obsrvd(nmp)   R*4  Measured spectrum (y)
c     nmp           I*4  Number of measured points (in spectrum)
c     apx(nfp)      R*4  A priori state vector (xa)
c     apu(nfp)      R*4  A priori state vector uncertainties (=SQRT[diag{Sa}])
c     slit(nii)     R*4  The ILS (slit function) needed by FM
c     nii           I*4  length of the SLIT vector
c     ldec          I*4  Decimation of SLIT vector
c     rdec          R*8  Ratio:  GINT/GRID (used by FM)
c     spxv(ncp,ntg) R*4  Extinction: VAC integrated along slant path
c     ncp           I*4  Number of primative spectral points (FM)
c     ntg           I*4  Number of target gases
c     nfp           I*4  Number of fitted parameters
c     snr           R*8  Nominal signal-to-noise ratio of OBSRVD
c     corrld        R*4  Noise amplification factor (non-independence of spectra points) 
c     ssssss        R*8  Offset from start of primitive spectrum (fraction of GINT)
c     debug         L    Activates diagnostic write-statements when .true.
c     mit           I*4  Maximum permitted number of iterations
c
c   Outputs:
c     nit           I*4  Number of iterations performed
c     calcul(nmp)   R*4  Calculated spectrum (f(x)) at final iteration
c     rmsocl        R*4  RMS difference between OBSRVD & CALCUL divided by CONT
c     cx(nfp)       R*4  State vector (x)
c     ex(nfp)       R*4  State vector uncertainties
c     pd(nmp,ntg)   R*4  Individual gas transmittance spectra
c
c   Notes/Explanations:
c    if nit>mit, then it didn't converge
c    if nit=mit, then it converged on the last iteration


      implicit none
      include "../ggg_const_params.f"
      include "const_params.f"
      include "int_params.f"

      logical
     & debug,    !  activates debug write-statements when .true.
     & cf        !  Fits channel fringes when .true.

      integer*4
     & ncbf,
     & iptg,ipcl,ipfs,ipzo,ipcf,
     & nconv,kconv, ! Number of convergences
     & ncp,         ! Number of precomputed absorption frequencies in SPXV
     & krank,
     & ierr,
     & nlev,
     & nmp,nii,ldec,nmpfp,
     & ntg,jtg,
     & nfp,kfp,jfp,i,j,
     & k,
     & mit,nit,
     & jva,jpd,knn

      real*4
c     & sflag,
     & slit(nii),
     & cx(nfp),dx(mfp),ex(nfp),
     & obsrvd(nmp),calcul(nmp),cont(nmp),resids(mmp+nfp),
     & apx(nfp),apu(nfp),
     & wk(mfp),
     & rmsocl, rmswas,
     & thresh, cont_level, cont_tilt, cont_curv, xzo,
     & ti,cc,oc,sum_oc,sum_toc,sum_cc,sum_tcc,sum_ttcc,denom,
     & ynoise,
     & tau,
     & corrld,
     & rnfac,
     & cmin,cmax,omin,omax,
     & var,
     & slpd(nmp,nlev),
     & vac(ncp,nlev,0:ntg),splos(nlev),ssnmp(nmp),
     & pd((mmp+mfp)*mfp),spts(ncp),
     & spxv(ncp*(ntg+3)), 
     & fs,fr,tt,eps,esat(0:mfp),
     & rdum,tbar,
     & dxlimit,xfr,
     & sumr2

      real*8
     & sh,sssss,
     & snr,rdec

      integer*4
     & ip(mfp)

c      logical prf_flag(ntg)

c      character winfo*(*),ss(nfp)*4
      parameter (tau=6.e-06,eps=0.01)

      cont_level=1.0
      cont_tilt=0.0
      cont_curv=0.0
      nmpfp=nmp+nfp

c      do j=1,ntg
c        prf_flag(j)=.false.
c      end do

c  Check that static array dimensions are adequate
      if(nfp.gt.mfp) then
         write(*,*)'mfp,nfp=',mfp,nfp
         stop 'do_retrieval: Increase parameter MFP'
      endif
      if(nmp.gt.mmp) then
         write(*,*)'mmp,nmp=',mmp,nmp
         stop 'do_retrieval: Increase parameter MMP'
      endif

      if( ipcf .gt. 0 ) then
         nconv=3
         cf=.true.
      else
         nconv=2
         cf=.false.
      endif

      if (debug) write(*,*)'do_retrieval: ntg,ncbf,nfp=',ntg,ncbf,nfp
      rmswas=big
      kconv=nconv
      do nit=0,mit     ! Spectral fitting iteration loop

c  Limit frequency shift to 1.8*GINT
         if(ipfs.gt.0) then
         if(abs(cx(ipfs)).gt.1.8) then
            cx(ipfs)=sign(1.8,cx(ipfs))
c            write(6,*)' Warning: Limiting Frequency shift'
         endif
         endif
c  Limit TILT if it exceeds 1.0
c        if( abs(cx(ipct)) .ge. 1.0 ) then
c          if(debug) write(6,*)' Limiting TILT'
c          cx(ipct)=sign(1.0,cx(ipct))
c        endif

c  Calculate spectrum & PD's
         if(debug) write(*,*)'do_retrieval calling fm: cx=',cx

         call fm(0,slit,nii,
     &   iptg,ipcl,ipfs,ipzo,ipcf,
     &   ldec,spts,spxv,
     &   vac,splos,nlev,ncp,rdec,sssss,
     &   cont_level,cont_tilt,cont_curv,xzo,
     &   cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,nmp)
c     &   cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,nmp,prf_flag,sflag)

         if (debug) write(*,*)'calcul=',calcul(1),calcul(2),
     &   calcul(nmp-1),calcul(nmp)
         call vdot(calcul,1,unity,0,tbar,nmp)
         tbar=tbar/nmp/(abs(cont_level)+0.0001)
c  In the zeroth iteration, estimate CL & CT using the ratio of
c  the measured to calculated spectra.
         if(ncbf.gt.0) then  ! Estimate CL and CT
         if (nit.eq.0) then
            sum_oc=eps
            sum_cc=eps
            sum_toc=eps
            sum_tcc=eps
            sum_ttcc=0.0
            rnfac=sqrt(3*float(nmp-1)/(nmp+1))
            do i=1,nmp
               ti=(1.0-2*float(i-1)/(nmp-1))*rnfac
               cc=calcul(i)*calcul(i)
               oc=obsrvd(i)*calcul(i)
               sum_oc=sum_oc+oc
               sum_toc=sum_toc+oc*ti
               sum_cc=sum_cc+cc
               sum_tcc=sum_tcc+cc*ti
               sum_ttcc=sum_ttcc+cc*ti**2
            end do
            denom=sum_tcc**2-sum_cc*sum_ttcc
            if(denom.ne.0.0 .and. ncbf.ge.2) then
               cont_level=(sum_tcc*sum_toc-sum_oc*sum_ttcc)/denom
               cont_tilt=(sum_tcc*sum_oc-sum_cc*sum_toc)/denom/
     &         cont_level
            else
               cont_level=sum_oc/sum_cc
               cont_tilt=0.0
            endif
c            write(*,*)
c     &    'do_retrieval: cx(ipcl),cx(ipct) =',cont_level,cont_tilt
            if(debug) write(*,*)
     &    'do_retrieval: cx(ipcl),cx(ipct) =',cont_level,cx(ipcl+1)
c
c  Apply the derived CL and CT to the calculated spectra and PD's
            do i=1,nmp
               ti=(1.0-2*float(i-1)/(nmp-1))*rnfac
               cont(i)=cont_level*(1+cont_tilt*ti)
               calcul(i)=calcul(i)*cont(i)
            end do
            call vmul(pd,1,cont_level,0,pd,1,nfp*(nmp+nfp))
         endif    !  (nit.eq.0)
         endif    !  (ncbf.gt.0)

         if(ncbf.ge.1) cx(ipcl)=cont_level
         if(ncbf.ge.2) cx(ipcl+1)=cont_tilt
         if(ncbf.ge.3) cx(ipcl+2)=cont_curv

         if(debug) then
            cmin=calcul(1)
            cmax=calcul(1)
            omin=obsrvd(1)
            omax=obsrvd(1)
            do i=2,nmp
               if(calcul(i).lt.cmin) cmin=calcul(i)
               if(calcul(i).gt.cmax) cmax=calcul(i)
               if(obsrvd(i).lt.omin) omin=obsrvd(i)
               if(obsrvd(i).gt.omax) omax=obsrvd(i)
            end do
            write(*,*) 'calmin, calmax=',cmin,cmax
            write(*,*) 'obsmin, obsmax=',omin,omax
         endif

c  Calculate residuals
         call compute_residual(kconv,obsrvd,calcul,resids,nmp)

c Calculate RMS fit
         call vdot(resids,1,resids,1,sumr2,nmp)
         rmsocl=sqrt(sumr2/nmp)/cont_level
         if(abs(rmsocl).gt.1.E+15) rmsocl=sign(1.E+15,rmsocl) ! prevents rmsocl=Inf

         if (rmsocl.gt.9*rmswas) then  ! Fit got much worse,
            if(debug)write(*,*)'nit,rmswas,rms='
     &         ,nit,rmswas,rmsocl,thresh,kconv
            if(debug) write(6,*)'Retracing 90%',kconv,rmsocl,rmswas
            call vmul(dx,1,-0.9,0,dx,1,nfp)
         else
            if(debug) write(6,*)'Continuing',kconv,rmsocl,rmswas
            thresh=(64*kconv-63)*(rmsocl+0.005)/100000
c            thresh=(64*kconv-63)*(rmsocl+0.01)/100000
            if(debug)write(*,*)'nit,rmswas,rms=',
     &        nit,rmswas,rmsocl,thresh,kconv
            if(abs(rmsocl-rmswas).lt.thresh.or.nit+2*kconv-1.gt.mit)then
               kconv=kconv-1
               rmswas=big
               if(kconv.ge.1) then
c                  if(cf) write(*,*) 'call subtract_cf',nit,kconv
                  if(cf) call subtract_cf(obsrvd,calcul,resids,nmp,mmp)
                  call vsub(obsrvd,1,calcul,1,resids,1,nmp) ! residuals
               endif
            endif ! abs(rmsocl-rmswas).lt.thresh .or. nit+2*kconv-1.gt.mit

            ynoise=2.5*cont_level*corrld/sngl(0.1d0+snr)
c   The last few (i>NMP) elements of RESID & PD contains A PRIORI info
c   RESIDS(nmp+i) holds the values of (AX(i)-CX(i))*YNOISE/APU(i)
c   PD(nmp+i,i) holds YNOISE/APU(i), with off-diagonal elements of zero.
            call vdiv(ynoise,0,apu,1,pd(nmp+1),nmpfp+1,nfp) ! YNOISE/APU(i)
            call vsub(apx,1,cx,1,wk,1,nfp)              ! APX(i)-CX(i)
            call vmul(wk,1,pd(nmp+1),nmpfp+1,resids(nmp+1),1,nfp)
            if(debug) then
              write(*,122)'apx=',(apx(j),j=ntg+1,nfp),(apx(j),j=1,ntg)
              write(*,122)'cx =',(cx(j),j=ntg+1,nfp),(cx(j),j=1,ntg)
122           format(a4,20f10.5)
            endif
c  Solve matrix equation PD.dx=resids
            call vmov(zero,0,wk,1,nfp)
c         write(*,*)'calling shfti...',nmpfp,nfp
            call shfti(pd,nmpfp,nmpfp,nfp,resids,nmpfp,
     &      1,tau,krank,rdum,wk,ip)
c            write(*,*)'called shfti...',(ip(j),j=1,nfp)

            call vmov(resids,1,dx,1,nfp)
            if(debug) then
               write(*,122)'dx =',(dx(jtg),jtg=ntg+1,nfp),
     &        (dx(jtg),jtg=1,ntg)
               if(krank.lt.nfp) then
                  write(6,*)'Rank Deficient:',krank,'  /',nfp
               else
                  write(6,*)'Full Rank:',krank,'  /',nfp
               endif
            endif
            if(ipcl.gt.0) then
            if(dx(ipcl).gt.9.21) then
               fs=10000.
            elseif(dx(ipcl).lt.-9.21) then
               fs=.0001
            else
               fs=exp(dx(ipcl))
            endif
            dx(ipcl)=cx(ipcl)*(fs-1)
            endif
         endif  ! (rmsocl.gt.9*rmswas)

         if(kconv.eq.0 .or. mit.eq.0) go to 63   ! convergence
c  Limit maximum step size for fitted gases (non-linear).
         fr=1.0
         do jfp=1,nfp
            dxlimit=0.25+0.5*abs(cx(jfp))
            xfr=abs(dxlimit/(tiny+abs(dx(jfp))))
            if(fr .gt. xfr) then
               fr=xfr
               kfp=jfp
            endif
         end do
         if(nit.eq.0) fr=min(fr,0.4)
         if(nit.eq.1) fr=min(fr,0.8)
         if(debug .and. fr.lt.1.)
     &   write(*,*)'Limited step size: kfp,fr=',kfp,fr
         call vsma(dx,1,fr,cx,1,cx,1,nfp)
         rmswas=rmsocl
      end do   ! do nit=0,mit     ! Spectral fitting iteration loop

c  Compute upper-triangle of covariance matrix & move diagonal elements into EX
63    var=(tbar*cont_level*rmsocl*corrld)**2/
     & (1.-float(nfp)/(float(nmp)+0.1))
      var=var+(tbar*tau)**2  ! fudge to prevent EX=0 when rmsocl=0
      call scov2(pd,nmpfp,nfp,ip,var,ierr)
      if(debug) then
         write(*,*) 'SCOV2: ierr=',ierr
         write(*,*)' Upper triangle of correlation matrix:'
         write(*,'(8(a,i1))')('  Target_',j,j=1,ntg),
     &   '   Cntuum     Tilt      Shift     ZOff'
         do i=1,nfp
            write(*,'(11f10.6)')(pd((j-1)*nmpfp+i)
     &   /sqrt(pd((j-1)*nmpfp+j)*pd((i-1)*nmpfp+i)),j=i,nfp)
         end do
      endif
      if(ierr.eq.0 .and. krank.gt.0) then
         call vmov(pd,nmpfp+1,ex,1,nfp)
      else
         call vmov(1/(tbar*tau)**2,0,ex,1,nfp)
      endif
      if(ipcl.gt.0) ex(ipcl)=ex(ipcl)*cx(ipcl)**2
c      write(6,*)'cx=',(cx(j),j=1,nfp)
c      write(6,*)'dx=',(dx(j),j=1,nfp)
c      write(6,*)'ex=',(sqrt(ex(j)),j=1,nfp)
c=========================================================================

c  Compute transmittances (convolved with ILS) of individual target gases.
c  Place them in PD, which just so happens to be exactly the right size.
      jva=1
      jpd=1
      knn=1+ncp*(ntg+1) ! Start address of workspace

c  SH is shift expressed as fraction of primative GRID.
      if(ipfs.gt.0) then
         sh=rdec*(cx(ipfs)+sssss)
      else
         sh=rdec*sssss
      endif
      do jtg=0,ntg
         if(jtg.eq.0) then ! non-target gases
            call vexp(spxv(jva),1,spxv(knn),1,ncp)
         else        ! target gases
            call vmul(spxv(jva),1,cx(jtg),0,spxv(knn),1,ncp)
            call vexp(spxv(knn),1,spxv(knn),1,ncp)
         endif
         call newdec(spxv(knn),ncp,slit,nii,ldec,rdec,sh,pd(jpd),nmp)
c
c  Compute saturation error (ESAT)
         tt=0.0
         do i=1,nmp
            tt=tt+1/(1+pd(i+jpd)/eps)
         end do
         esat(jtg)=((tt/nmp-1/(1+1/eps))**2)/eps

         jva=jva+ncp
         jpd=jpd+nmpfp
      end do

c  Do the solar spectrum too.
      call newdec(spts,ncp,slit,nii,ldec,rdec,sh,ssnmp,nmp)

c  Determine the average optical depth of the first target gas

c  Add RSS error contributions from change in zero level (ZERR) and RMS fit
c      call vmul(dx,1,dx,1,dx,1,nfp)
c      call vadd (dx,1,ex,1,ex,1,nfp)
c      call vadd (ex,1,(3*rmsocl)**2,0,ex,1,nfp)   !  This is a fudge
c      call vsqrt(ex,1,ex,1,nfp)
      if(abs(cont_level).lt.tiny) cont_level=tiny
c      write(*,*)'ex=',sqrt(abs(ex(1))),(3*rmsocl),
c     & 5*mit*dx(1)/(mit+1),tbar
      do jfp=1,nfp
         ex(jfp)=dsqrt(
     &  dabs(dble(ex(jfp)))/(abs(tbar)+0.0001d0)/(tbar**2+1.0E-36)
c     &     +dx(jfp)**2      ! perturbation from adding ERROFF to residuals
     &     +(3*rmsocl)**2   ! fudge
c     &     +(100*(cx(jfp)-apx(jfp))*rmsocl)**2   ! fudge
c     &     +25*dxwas(jfp)**2   ! fudge
     &     +(5*mit*dx(jfp)/(mit+0.0001))**2   ! fudge
c     &     +0.5*esat(jfp)**2   ! fudge to increase error for saturated target lines
c     &     +0.04*(cx(ipcl)**2+1./cx(ipcl)**2)*(cx(jfp)-apx(jfp))**4  ! fudge
     &     )
      end do
      return
      end
