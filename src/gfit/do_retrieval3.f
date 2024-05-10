      subroutine do_retrieval3(obsrvd,nmp,apx,apu,slit,nhw,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,xzo,
     & cont_level,cont_tilt,cont_curv,
     & cfamp,cfperiod,cfphase,
     & nfov,z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     & ldec,rdec,spts,spxv,vac,splos,nlev,ncp,ntg,ncbf,nfp,snr,
     & corrld,debug,mit,nit,cont,calcul,rmsocl,cx,ex,
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
c     slit(1+2*ldec*nhw)     R*4  The ILS (slit function) needed by FM
c     nhw           I*4  Half-width of the SLIT vector
c     ldec          I*4  Decimation of SLIT vector
c     rdec          R*8  Ratio:  GINT/GRID (used by FM)
c     spxv(ncp,ntg) R*4  Extinction: VAC integrated along slant path
c     ncp           I*4  Number of primative spectral points (FM)
c     ntg           I*4  Number of target gases
c     nfp           I*4  Number of fitted parameters
c     snr           R*8  Nominal signal-to-noise ratio of OBSRVD
c     corrld        R*4  Noise amplification factor (non-independence of spectra points)
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


      implicit none

      real*4 zero, unity
      parameter(zero=0.0)
      parameter(unity=1.0)

      include "const_params.f"
      include "int_params.f"

      logical
     & debug,    !  activates debug write-statements when .true.
     & cf        !  Fits channel fringes when .true.

      integer*4
     & ncbf,idum,
     & nexpl,nexpr,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & nconv,kconv, ! Number of convergences
     & ncp,         ! Number of precomputed absorption frequencies in SPXV
     & krank,
     & ierr,
     & nlev,
     & nfov,
     & nmp,nhw,ldec,nmpfp,
     & ntg,jtg,
     & nfp,kfp,jfp,i,j,
     & mit,nit,
     & jva,jpd,knn

      real*4
     & cfamp,cffreq,cfperiod,cfphase,tcbar,spi,
     & solzen,roc,fbar,
     & slit(1+2*ldec*nhw),
     & cx(nfp),dx(mfp),ex(nfp),
     & obsrvd(nmp),calcul(nmp),cont(nmp),resid(mmp+nfp),
     & tcalc(mmp),
     & apx(nfp),apu(nfp),
     & wk(mfp),
     & rmsocl, rmswas,
     & thresh, cont_level, cont_tilt, cont_curv, xzo,
     & ti,cc,oc,sum_oc,sum_toc,sum_cc,sum_tcc,sum_ttcc,denom,
     & ynoise,
     & tau,
     & corrld,
     & cmin,cmax,omin,omax,
     & var,
     & slpd(nmp,nlev),
     & z(nlev),t(nlev),p(nlev),
     & vac(ncp,nlev,0:ntg),splos(nlev),ssnmp(nmp),
     & pd((mmp+mfp)*mfp),spts(ncp),
     & tpd((mmp+mfp)*mfp),
     & spxv(ncp*(ntg+3)), 
     & fs,fr,tt,eps,esat(0:mfp),
     & rdum,tbar,
     & dxlimit,xfr,
     & sumr2

      real*8
     & wavtkr,obalt,fovo,
     & hh,snr,rdec

      integer*4
     & ip(mfp)

c      character winfo*(*),ss(nfp)*4
      parameter (tau=7.e-05,eps=0.01,spi=3.14159265)


      rdum=ckm2cm
      idum=mtg

      cont_level=1.0
      cont_tilt=0.0
      cont_curv=0.0
      nmpfp=nmp+nfp
      cfamp=0.0
      cfperiod=0.0
      cfphase=0.0

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
c               write(6,*)' Warning: Limiting Frequency shift'
            endif
         endif
c  Limit TILT if it exceeds 1.0
c        if( abs(cx(ipct)) .ge. 1.0 ) then
c          if(debug) write(6,*)' Limiting TILT'
c          cx(ipct)=sign(1.0,cx(ipct))
c        endif

c  Calculate spectrum & PD's
         if(debug) write(*,*)'do_retrieval3 calling fm3: cx=',cx
         call fm3(0,slit,nhw,
     &   ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     &   ldec,spts,spxv,
     &   nfov,z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     &   vac,splos,nlev,ncp,rdec,
     &   cont_level, cont_tilt, cont_curv,xzo,
     &   cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,
     &   tcalc,tpd,nmp)

         if(debug) write(*,*)'calcul=',calcul(1),calcul(2),
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
               do i=1,nmp
                  ti=1.0-2*float(i-1)/(nmp-1)
                  cc=calcul(i)*calcul(i)
                  oc=obsrvd(i)*calcul(i)
                  sum_oc=sum_oc+oc
                  sum_toc=sum_toc+oc*ti
                  sum_cc=sum_cc+cc
                  sum_tcc=sum_tcc+cc*ti
                  sum_ttcc=sum_ttcc+cc*ti**2
               end do
               denom=sum_tcc**2-sum_cc*sum_ttcc
               if(abs(denom).gt.0.0 .and. ncbf.ge.2) then
                  cont_level=(sum_tcc*sum_toc-sum_oc*sum_ttcc)/denom
                  cont_tilt=(sum_tcc*sum_oc-sum_cc*sum_toc)/denom/
     &            cont_level
               else
                  cont_level=sum_oc/sum_cc
                  cont_tilt=0.0
               endif
               if(debug) write(*,*)
     &       'do_retrieval: cx(ipcl),cx(ipct) =',cont_level,cx(ipcl+1)
c
c  Apply the derived CL and CT to the calculated spectra and PD's
               do i=1,nmp
                  ti=1.0-2*float(i-1)/(nmp-1)
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
            write(*,*) 'cmin, cmax=',cmin,cmax
            write(*,*) 'omin, omax=',omin,omax
         endif

c  Calculate residuals
         call compute_residual(kconv,obsrvd,calcul,resid,nmp)

c Calculate RMS fit
         call vdot(resid,1,resid,1,sumr2,nmp)
         rmsocl=sqrt(sumr2/nmp)/cont_level
         if(abs(rmsocl).gt.1.E+15) rmsocl=sign(1.E+15,rmsocl) ! prevents rmsocl=Inf

         if (rmsocl.gt.9*rmswas) then  ! Fit got much worse,
            if(debug) write(6,*)'Retracing step',rmsocl,rmswas
            call vmul(dx,1,-0.9,0,dx,1,nfp)
         else
            if(debug) write(6,*)'Continuing',rmsocl,rmswas
            thresh=(64*kconv-63)*abs(cont_level)*(rmsocl+0.01)/100000
            if(abs(rmsocl-rmswas).lt.thresh.or.nit+2*kconv-1.gt.mit)then
               kconv=kconv-1
               rmswas=big
               if(kconv.ge.1) then
                  if(cf) then
c                     write(*,*) 'calling fit_channel_fringe',nit,kconv

                     call fit_channel_fringe(mmp,nmp,ncbf,resid,
     &                                           cfamp,cfperiod,cfphase)
c                     write(*,*)'cfamp,period,phase = ',cfamp,cfperiod,cfphase

c  If fringes are found, remove them from OBSRVD.
                     if(cfamp*cfperiod.gt.0.0) then   ! GCT 2016-05-12
                        call vdot(calcul,1,unity,0,tcbar,nmp)
                        tcbar=tiny+tcbar/nmp
                        cffreq=2*spi*float(nmp-1)/cfperiod
                        call vramp(resid,1,nmp)
                        call vsma(resid,1,cffreq,cfphase,0,resid,1,nmp)
                        call vcos(resid,1,resid,1,nmp)
                        call vmul(resid,1,calcul,1,resid,1,nmp)
                        call vsma(resid,1,-(cfamp/tcbar),obsrvd,1,
     &                  obsrvd,1,nmp)
                     endif
                     call vsub(obsrvd,1,calcul,1,resid,1,nmp) !  residuals
                  endif  !  if (cf)
               endif  !  if(kconv.ge.1) then
            endif ! abs(rmsocl-rmswas).lt.thresh .or. nit+2*kconv-1.gt.mit

            ynoise=2.5*cont_level*corrld/sngl(0.1d0+snr)
c   The last few (i>NMP) elements of RESID & PD contains A PRIORI info
c   RESIDS(nmp+i) holds the values of (AX(i)-CX(i))*YNOISE/APU(i)
c   PD(nmp+i,i) holds YNOISE/APU(i), with off-diagonal elements of zero.
            call vdiv(ynoise,0,apu,1,pd(nmp+1),nmpfp+1,nfp) ! YNOISE/APU(i)
            call vsub(apx,1,cx,1,wk,1,nfp)              ! APX(i)-CX(i)
            call vmul(wk,1,pd(nmp+1),nmpfp+1,resid(nmp+1),1,nfp)
            if(debug) then
               write(*,122)'apx=',(apx(j),j=ntg+1,nfp),(apx(j),j=1,ntg)
               write(*,122)'cx =',(cx(j),j=ntg+1,nfp),(cx(j),j=1,ntg)
122            format(a4,20f10.5)
            endif
c  Solve matrix equation PD.dx=resid
            call vmov(zero,0,wk,1,nfp)
c         write(*,*)'calling shfti...',nmpfp,nfp
            call shfti(pd,nmpfp,nmpfp,nfp,resid,nmpfp,
     &      1,tau,krank,rdum,wk,ip)
c            write(*,*)'called shfti...',(ip(j),j=1,nfp)

            call vmov(resid,1,dx,1,nfp)
            if(debug) then
               write(*,122)'dx =',(dx(jtg),jtg=ntg+1,nfp),
     &         (dx(jtg),jtg=1,ntg)
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
            dxlimit=0.2+abs(cx(jfp))
            xfr=abs(dxlimit/(tiny+abs(dx(jfp))))
            if(fr .gt. xfr) then
               fr=xfr
               kfp=jfp
            endif
         end do
         if(nit.eq.0) fr=min(fr,0.3)
         if(nit.eq.1) fr=min(fr,0.6)
         if(debug .and. fr.lt.1.)
     &   write(*,*)'Limited step size: kfp,fr=',kfp,fr
         call vsma(dx,1,fr,cx,1,cx,1,nfp)
         rmswas=rmsocl
      end do   ! do nit=0,mit     ! Spectral fitting iteration loop

c  Compute upper-triangle of covariance matrix & move diagonal elements into EX
63    var=(cont_level*rmsocl*corrld)**2/(1.-float(nfp)/(float(nmp)+0.1))
      var=var+tau**2  ! fudge to prevent EX=0 when rmsocl=0
      call scov2(pd,nmpfp,nfp,ip,var,ierr)
      if(debug) then
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
         call vmov(1/tau**2,0,ex,1,nfp)
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
      if(ipfs.gt.0) then
         hh=rdec*cx(ipfs)/(ifcsp+0.5*(ncp-1))
      else
         hh=0.0d0
      endif

      do jtg=0,ntg
         if(jtg.eq.0) then ! non-target gases
            call vexp(spxv(jva),1,spxv(knn),1,ncp)
         else        ! target gases
            call vmul(spxv(jva),1,cx(jtg),0,spxv(knn),1,ncp)
            call vexp(spxv(knn),1,spxv(knn),1,ncp)
         endif
c         write(*,*)'do_retrieval3: calling regrid2:'
         call regrid2(ifcsp,ncp,spxv(knn),nhw,slit,ldec,
     &    rdec*(1.0d0+hh),ifmsp,nmp,pd(jpd),nexpl,nexpr)
         if(nexpl.ne.0)write(*,*)'Warning: do_retrieval3: NEXPL=',nexpl
         if(nexpr.ne.0)write(*,*)'Warning: do_retrieval3: NEXPR=',nexpr
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
      call regrid2(ifcsp,ncp,spts,nhw,slit,ldec,
     & rdec*(1.0d0+hh),ifmsp,nmp,ssnmp,nexpl,nexpr)
      if(nexpl.ne.0) write(*,*)'Warning: do_retrieval3: NEXPL=',nexpl
      if(nexpr.ne.0) write(*,*)'Warning: do_retrieval3: NEXPR=',nexpr

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
         ex(jfp)=sqrt(
     &   abs(ex(jfp)/(tbar+0.0001))
c     &  +dx(jfp)**2      ! perturbation from adding ERROFF to residuals
     &  +(3*rmsocl)**2   ! fudge
c     &  +(100*(cx(jfp)-apx(jfp))*rmsocl)**2   ! fudge
c     &  +25*dxwas(jfp)**2   ! fudge
     &  +(5*mit*dx(jfp)/(mit+0.0001))**2   ! fudge
c     &  +0.5*esat(jfp)**2   ! fudge to increase error for saturated target lines
c     &  +0.04*(cx(ipcl)**2+1./cx(ipcl)**2)*(cx(jfp)-apx(jfp))**4  !  fudge
     &   )
      end do
      return
      end
