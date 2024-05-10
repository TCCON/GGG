      subroutine do_retrieval(obsrvd,nmp,apx,apu,slit,nhw,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,xzo,
     & cont_level,cont_tilt,cont_curv,
     & cfamp,cfperiod,cfphase,
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
c     debug         L    Activates diagnostic write-statements when .true.
c     mit           I*4  Maximum permitted number of iterations
c
c   Outputs:
c     nit           I*4  Number of iterations performed
c     calcul(nmp)   R*4  Calculated spectrum (f(x)) at final iteration
c                        CALCUL = CONT*[TRAN*(1-XZO)+XZO]
c     rmsocl        R*4  RMS difference between OBSRVD & CALCUL divided by CONT
c     cx(nfp)       R*4  State vector (x)
c     ex(nfp)       R*4  State vector uncertainties
c     pd(nmp,ntg)   R*4  Individual gas transmittance spectra
c
c   Notes/Explanations:
c    if nit>mit, then it didn't converge
c    if nit=mit, then it converged on the last iteration


      implicit none
      real*4  zero, unity, spi
      parameter(zero=0.0)
      parameter(unity=1.0)
      parameter(spi=3.14159265)

      include "const_params.f"
      include "int_params.f"

      logical
     & debug,    !  activates debug write-statements when .true.
     & cf        !  Fits channel fringes when .true.

      integer*4
     & ncbf, idum,
     & nexpl,nexpr,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & nconv,kconv, ! Number of convergences
     & ncp,         ! Number of precomputed absorption frequencies in SPXV
     & krank,
     & ierr,
     & nlev,
     & nmp,nhw,ldec,nmpfp,
     & ntg,jtg,
     & nfp,kfp,jfp,i,j,
     & mit,nit,
     & jva,jpd,knn

      real*4
c     & sflag,
     & slit(1+2*nhw*ldec),
     & cx(nfp),dx(mfp),ex(nfp),
     & obsrvd(nmp),calcul(nmp),cont(nmp),resid(mmp+nfp),
     & cfamp,cffreq,cfperiod,cfphase,tcbar,
     & apx(nfp),apu(nfp),
     & wk(mfp),
     & rmsocl, rmswas,
     & thresh, cont_level, cont_tilt, cont_curv, xzo, tilt_new,
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
     & fs,fr,eps,
     & rdum,tbar,
     & dxlimit,xfr,
     & sumr2

      real*8
     & xfs, snr, rdec

      integer*4
     & ip(mfp)

c      logical prf_flag(ntg)

      parameter (tau=4e-04,eps=0.01)

      rdum=ckm2cm
      idum=mtg

      cont_level=1.0
      cont_tilt=0.0
      cont_curv=0.0
      nmpfp=nmp+nfp
      cfamp=0.0
      cfperiod=0.0
      cfphase=0.0

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

      nconv=2
      if( ipcf .gt. 0 ) then
         cf=.true.
      else
         cf=.false.
      endif

      if (debug) write(*,*)'do_retrieval: ntg,ncbf,nfp=',ntg,ncbf,nfp
c      write(*,*)'inside do_retrieval: ntg,ncbf,nfp=',ntg,ncbf,nfp
      rmswas=big
      kconv=nconv
      do nit=0,mit     ! Spectral fitting iteration loop

c  Limit frequency shift to (gint+grid)/freq
         if(ipfs.gt.0) then
            if(abs(cx(ipfs)).gt.8E-06) cx(ipfs)=sign(8E-06,cx(ipfs))
         endif

c  Limit SG shift to 6 ppm
         if(ipsg.gt.0) then
            if(abs(cx(ipsg)).gt.6.E-06) cx(ipsg)=sign(6.E-06,cx(ipsg))
         endif

c  Limit TILT if it exceeds 1.0/Sqrt(3)
c  The TILT basis function has end values of +/- Sqrt(3).
         if( abs(cont_tilt) .gt. 1.0/1.732 ) then
            tilt_new=sign(1.0/1.732,cont_tilt)
c            write(6,*)' Limiting TILT: ',cont_tilt,tilt_new
            cont_tilt=tilt_new
            if(ncbf.gt.1) cx(ipcl+1)=cont_tilt
         endif

c  Calculate spectrum & PD's
         if(debug) write(*,*)'do_retrieval calling fm: cx=',cx
c         write(*,*)'calling fm: cx=',nit,(cx(j),j=1,nfp)

         call fm(0,slit,nhw,
     &   ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     &   ldec,spts,spxv,
     &   vac,splos,nlev,ncp,rdec,
     &   cont_level,cont_tilt,cont_curv,xzo,
     &   cx,ntg,ncbf,nfp,cont,calcul,slpd,pd,nmp)

c         write(*,*)'do_retrieval: calcul=',calcul(1),calcul(2),
c      &   calcul(nmp-1),calcul(nmp)
         if(debug)write(*,*)'do_retrieval: calcul=',calcul(1),calcul(2),
     &   calcul(nmp-1),calcul(nmp)

c    CALCUL = CONT*[TRAN*(1-XZO)+XZO]
c    TRANS = [CALCUL/CONT - XZO]/(1-XZO)
         call vdot(calcul,1,unity,0,tbar,nmp)
         tbar=tbar/nmp/(abs(cont_level)+tiny)
         tbar=(tbar-xzo)/(1-xzo)
c         write(*,*)'tbar,cx=',tbar ,(cx(j),j=1,nfp)
c  In the zeroth iteration, estimate CL & CT using the ratio of
c  the measured to calculated spectra.
         if(ncbf.gt.0) then  ! Estimate CL and CT
            if( nit.eq.0) then
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
c               write(*,*) 'do_retrieval5:',
c     &         sum_oc,sum_toc,sum_tcc,sum_ttcc,denom
               if(abs(denom).gt.0.0 .and. ncbf.ge.2) then
                  cont_level=(sum_tcc*sum_toc-sum_oc*sum_ttcc)/denom
                  cont_tilt=(sum_tcc*sum_oc-sum_cc*sum_toc)/denom/
     &            cont_level
               else
                  cont_level=sum_oc/sum_cc
                  cont_tilt=0.0
               endif
               if(cont_level .gt. 1.E+25) cont_level=1.E+25
               if(debug) write(*,*)
     &       'do_retrieval: cont_level,cont_tilt= ',cont_level,cont_tilt
c  Apply the derived CL and CT to the calculated spectra and PD's
               do i=1,nmp
                  ti=(1.0-2*float(i-1)/(nmp-1))*rnfac
                  cont(i)=cont_level*(1+cont_tilt*ti)
                  calcul(i)=calcul(i)*cont(i)
               end do
c               call vdot(pd,1,1.0,0,sumr2,nfp*(nmp+nfp))
c               write(*,*)'sum_pd1=',sumr2
c               write(*,*)'cont_level=',cont_level

               call vmul(pd,1,cont_level,0,pd,1,nfp*(nmp+nfp))

c               call vdot(pd,1,1.0,0,sumr2,nfp*(nmp+nfp))
c               write(*,*)'sum_pd2=',sumr2

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
         call compute_residual(kconv,obsrvd,calcul,resid,nmp)

c Calculate RMS fit
c         call vdot(obsrvd,1,obsrvd,1,sumr2,nmp)
c         write(*,*)'sumr2_obsrvd = ',sumr2
c         call vdot(calcul,1,calcul,1,sumr2,nmp)
c         write(*,*)'sumr2_calcul = ',sumr2
         call vdot(resid,1,resid,1,sumr2,nmp)
         rmsocl=sqrt(sumr2/nmp)/cont_level
         if(abs(rmsocl).gt.1.E+15) rmsocl=sign(1.E+15,rmsocl) ! prevents rmsocl=Inf
c         write(*,*)'sumr2_resid,cont_level = ',sumr2,cont_level,rmsocl

c         write(*,*) 'i1  obsrvd(i)  calcul(i)  resid(i) '
c         do i=1,nmp
c            write(*,*)i,obsrvd(i),calcul(i),resid(i)
c         end do

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
                  if(cf) then
c                     write(*,*) 'calling fit_channel_fringe',nit,kconv

                     call fit_channel_fringe(mmp,nmp,ncbf,resid,
     &                                           cfamp,cfperiod,cfphase)
c                     write(*,*)'cfamp,period,phase = ',cfamp,cfperiod,
c     &               cfphase,float(nmp-1)/cfperiod

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
                     call vsub(obsrvd,1,calcul,1,resid,1,nmp) ! residuals
                  endif  ! if (cf)
               endif  ! if (kconv.eq.1)
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
               write(*,*) 'resid=', (resid(j),j=nmp+1,nmp+ntg)
122            format(a4,20f10.5)
            endif
c  Solve matrix equation PD.dx=resid
            call vmov(zero,0,wk,1,nfp)
c            write(*,*)'calling shfti...',nmpfp,nfp
c            call vdot(pd,1,1.0,0,sumr2,nmpfp*nfp)
c            write(*,*)'Sum PD3 =',sumr2
c            write(*,*)'PD(),resid ='
c            do i=1,nmpfp
c              write(*,*) (pd(i+nmpfp*(j-1)),j=1,nfp),resid(i)
c            end do

            call shfti(pd,nmpfp,nmpfp,nfp,resid,nmpfp,
     &      1,tau,krank,rdum,wk,ip)
c            write(*,*)'called shfti...',krank,nfp

            call vmov(resid,1,dx,1,nfp)
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

         if(kconv.eq.0 .or. mit.eq.0) go to 66   ! convergence
c  Limit maximum step size for fitted gases (non-linear).
         fr=1.0
         do jfp=1,nfp
            dxlimit=0.2+0.5*abs(cx(jfp))
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
c         write(*,*)'nit, fr,cx(1),dx(1)=',nit,fr,cx(1),dx(1)
         call vsma(dx,1,fr,cx,1,cx,1,nfp)
         rmswas=rmsocl
      end do   ! do nit=0,mit     ! Spectral fitting iteration loop

c  Compute upper-triangle of covariance matrix & move diagonal elements into EX
66    var=((tbar+1.E-08)*cont_level*rmsocl*corrld)**2/
     & (1.-float(nfp)/(float(nmp)+0.1))
      var=var+((tbar+1.E-08)*tau)**2  ! fudge to prevent EX=0 when rmsocl=0
c      write(*,*) 'tbar,cont_level,rmsocl,var=',
c     & tbar,cont_level,rmsocl,corrld,var
      call scov2(pd,nmpfp,nfp,ip,var,ierr)
      if(debug) then
         write(*,*) 'SCOV2: ierr=',ierr
         write(*,*)' Correlation Coefficient Matrix:'
         write(*,'(8x,8(a,i1))')('  Target_',j,j=1,ntg),
     &   '   Cntuum     Tilt        FS       SG'
         do i=1,nfp
            write(*,'(a7,i1,32f10.6)') 'Target_',i,(pd((i-1)*nmpfp+j)
     &   /sqrt(pd((j-1)*nmpfp+j)*pd((i-1)*nmpfp+i)),j=1,i-1),
     &   (pd((j-1)*nmpfp+i)
     &   /sqrt(pd((j-1)*nmpfp+j)*pd((i-1)*nmpfp+i)),j=i,nfp)
         end do
      endif
      if(ierr.eq.0 .and. krank.gt.0) then
         call vmov(pd,nmpfp+1,ex,1,nfp)
      else
         call vmov(1/((tbar*tau)**2+1.E-36),0,ex,1,nfp)
      endif
      do j=1,nfp
         if(ex(j).gt.1.E+36) ex(j)=1.E+36
      end do
      if(ipcl.gt.0) ex(ipcl)=ex(ipcl)*cx(ipcl)**2
c      write(6,*)'cx=',(cx(j),j=1,nfp)
c      write(6,*)'dx=',(dx(j),j=1,nfp)
c      write(6,*)'ex=',(ex(j),j=1,nfp)
c      write(6,*)'sr(ex)=',(sqrt(ex(j)),j=1,nfp)
c=========================================================================

c  FS is shift expressed as fraction of primative GRID.
      if(ipfs.gt.0) then
         xfs=dble(cx(ipfs))
      else
         xfs=0.0d0
      endif

c  Compute transmittances (convolved with ILS) of individual target gases.
c  These will be written to the spt files
c  Place them in PD, which just so happens to be exactly the right size.
      jva=1
      jpd=1
      knn=1+ncp*(ntg+1) ! Start address of workspace
      do jtg=0,ntg
         if(jtg.eq.0) then ! non-target gases
            call vexp(spxv(jva),1,spxv(knn),1,ncp)
         else        ! target gases
            call vmul(spxv(jva),1,cx(jtg),0,spxv(knn),1,ncp)
            call vexp(spxv(knn),1,spxv(knn),1,ncp)
         endif
c       write(*,*) 'do_retrieval: regrid2: ',rdec,xfs
         call regrid2(ifcsp,ncp,spxv(knn),nhw,slit,ldec,
     &    rdec*(1.0d0+xfs),ifmsp,nmp,pd(jpd),nexpl,nexpr)
         if(nexpl.ne.0) write(*,*)'Warning: extrapolating NEXPL=',nexpl
         if(nexpr.ne.0) write(*,*)'Warning: extrapolating NEXPR=',nexpr

         jva=jva+ncp
         jpd=jpd+nmpfp
      end do  ! jtg=0,ntg

c  Do the solar spectrum too.
c       write(*,*) 'do_retrieval: solar: regrid2: ',rdec,xfs
      call regrid2(ifcsp,ncp,spts,nhw,slit,ldec,
     &rdec*(1.0d0+xfs),ifmsp,nmp,ssnmp,nexpl,nexpr)
      if(nexpl.ne.0) write(*,*)'Warning: extrapolating NEXPL=',nexpl
      if(nexpr.ne.0) write(*,*)'Warning: extrapolating NEXPR=',nexpr

c  Determine the average optical depth of the first target gas

      if(abs(cont_level).lt.tiny) cont_level=tiny
c      write(*,*)'tbar,ex=',tbar,ex(1),sqrt(abs(ex(1))),(3*rmsocl),
c     & 5*mit*fr*dx(1)/(mit+1),dx(1)
      do jfp=1,nfp
         ex(jfp)=sngl(dsqrt(dabs(dble(ex(jfp)))/(abs(tbar)+0.0001d0)/
     &   (tbar**2+1.D-36)
c     &     +dx(jfp)**2      ! perturbation from adding ERROFF to residuals
     &   +(3*rmsocl)**2   ! fudge
c     &     +(100*(cx(jfp)-apx(jfp))*rmsocl)**2   ! fudge
c     &     +25*dxwas(jfp)**2   ! fudge
     &   +(5*mit*fr*dx(jfp)/(mit+0.0001))**2   ! fudge
c     &     +0.04*(cx(ipcl)**2+1./cx(ipcl)**2)*(cx(jfp)-apx(jfp))**4  ! fudge
     &   ))
      end do
      return
      end
