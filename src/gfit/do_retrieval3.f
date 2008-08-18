      subroutine do_retrieval3(obsrvd,nmp,apx,apu,slit,nii,
     & z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     & ldec,rdec,spts,spxv,dspdzxv,vac,splos,nlev,ncp,ntg,snr,
     & corrld,sssss,winfo,debug,mit,nit,calcul,rms,cx,ex,pd,ssnmp)

c  Adjusts the state vector (cx) to get the best fit to measured spectrum (y).
c  Performs an iterative solution to the optimal estimation equation
c  (K'.S-1.K+Sa-1)dx = (K'.Sy-1.(y-f(x))+Sa-1.(x-xa))              
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
c     apu(nfp)      R*4  A priori state vector uncertainties (Sa)
c     slit(nii)     R*4  The ILS (slit function) needed by FM
c     nii           I*4  length of the SLIT vector
c     ldec          I*4  Decimation of SLIT vector
c     rdec          R*8  Ratio:  GINT/GRID (used by FM)
c     spxv(ncp,ntg) R*4  Extinction: VAC integrated along slant path
c  dspdzxv(ncp,ntg) R*4  Extinction: VAC integrated along slant path
c     ncp           I*4  Number of primative spectral points (FM)
c     ntg           I*4  Number of target gases
c     snr           R*8  Nominal signal-to-noise ratio of OBSRVD
c     corrld        R*4  
c     ssssss        R*8  Offset from start of primitive spectrum
c     winfo         C**  Command line from .ggg file
c     debug         L    Activates diagnostic write-statements when .true.
c     mit           I*4  Maximum permitted number of iterations
c
c   Outputs:
c     nit           I*4  Number of iterations performed
c     calcul(nmp)   R*4  Calculated spectrum (f(x)) at final iteration
c     rms           R*4  RMS difference between OBSRVD & CALCUL
c     cx(nfp)       R*4  State vector (x)
c     ex(nfp)       R*4  State vector uncertainties
c     pd(nmp,ntg)   R*4  Individual gas transmittance spectra


      implicit none

      logical
     & debug,    !  activates debug write-statements when .true.
     & cf        !  Fits chennel fringes when .true.

      integer*4
     & n1,n2,n3,n4,
     & nconv,kconv, ! Number of convergences
     & ncp,         ! Number of precomputed absorption frequencies in SPXV
     & krank,
     & ierr,
     & nlev,
     & mmp,nmp,imp,nii,ldec,nmpfp,
     & mfp,ntg,jtg,
     & nfp,kfp,jfp,i,j,
     & mit,nit,
     & jva,jpd,kn2

      parameter (mmp=360000,mfp=14)

      real*4
     & solzen,roc,fbar,
     & slit(nii),
     & cx(ntg+4),dx(mfp),ex(ntg+4),
     & obsrvd(nmp),calcul(nmp),resids(mmp+ntg+4),
     & tcalc(mmp),
     & apx(ntg+4),apu(ntg+4),
     & wk(mfp),
     & rms, rwas,
     & thresh,
     & ynoise,
     & tau,
     & corrld,
     & var,
     & zero,unity,
     & z(nlev),t(nlev),p(nlev),
     & vac(ncp,nlev,0:ntg),splos(nlev),ssnmp(nmp),
     & pd((mmp+mfp)*mfp),
     & tpd((mmp+mfp)*mfp),
     & spts(ncp),
     & spxv(ncp*(ntg+3)), dspdzxv(ncp*(ntg+3)),
     & tiny,big,
     & fs,fr,tt,eps,esat(0:mfp),
     & rdum,
     & dxlimit,xfr,
     & sumr2

      real*8
     & wavtkr,obalt,fovo,
     & sh,sssss,
     & snr,rdec

      integer*4
     & ip(mfp)

      character winfo*(*)
      parameter (zero=0.0,unity=1.0,tiny=1.0e-36,tau=6.e-06,big=1.0E+18,
     & eps=0.01)
c
      n1=ntg+1
      n2=ntg+2
      n3=ntg+3
      n4=ntg+4
      nfp=ntg+4
      nmpfp=nmp+nfp

c      call vdot(obsrvd,1,obsrvd,1,sumr2,nmp)
c      if(debug) write(*,*)'obsrvd=',sqrt(sumr2/nmp)

c  Check that static array dimensions are adequate
      if(nfp.gt.mfp) then
         write(*,*)'mfp,nfp=',mfp,nfp
         stop 'do_retrieval: Increase parameter MFP'
      endif
      if(nmp.gt.mmp) then
         write(*,*)'mmp,nmp=',mmp,nmp
         stop 'do_retrieval: Increase parameter MMP'
      endif

      if( index(winfo,' cf ') .gt. 0 ) then
         nconv=3
         cf=.true.
      else
         nconv=2
         cf=.false.
      endif

      rwas=big
      kconv=nconv
      do nit=0,mit     ! Spectral fitting iteration loop

c  Limit frequency shift to 0.8*GINT
         if(abs(cx(n3)) .gt. 0.8) then
            cx(n3)=sign(0.8,cx(n3))
c            write(6,*)' Warning: Limiting Frequency shift'
         endif
c  Limit TILT if it exceeds 1.0
c        if( abs(cx(n2)) .ge. 1.0 ) then
c          if(debug) write(6,*)' Limiting TILT'
c          cx(n2)=sign(1.0,cx(n2))
c        endif

c  Calculate spectrum & PD's
c         write(*,*)'do_retrieval calling fm: cx=',cx
         call fm3(0,winfo,slit,nii,ldec,spts,spxv,dspdzxv,
     &   z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     &   vac,splos,nlev,ncp,rdec,sssss,cx,ntg,calcul,pd,
     &   tcalc,tpd,nmp)
         call vdot(calcul,1,calcul,1,sumr2,nmp)
         if(debug) write(*,*)'calcul=',sqrt(sumr2/nmp)

c  Calculate residuals
c         if(kconv.eq.nconv) then ! use logarithmic residuals 1'st convergence
c            do imp=1,nmp
c               if(calcul(imp).le.tiny*obsrvd(imp).or.
c     &         obsrvd(imp).le.tiny*calcul(imp)) then
c                  resids(imp)=obsrvd(imp)-calcul(imp)
c               else
c                  resids(imp)=calcul(imp)*log(obsrvd(imp)/calcul(imp))
c               endif
c            end do   ! imp=1,nmp
c         else  !                 use linear residuals for subsequent convergences
c            call vsub(obsrvd,1,calcul,1,resids,1,nmp) ! residuals
c         endif   !  kconv.eq.nconv
         call compute_residual(kconv,obsrvd,calcul,resids,nmp)

c Calculate RMS fit
         call vdot(resids,1,resids,1,sumr2,nmp)
         rms=sqrt(sumr2/nmp)
         if(debug) write(*,*)'%rms=',100*rms
         if(abs(rms).gt.1.E+15) rms=sign(1.E+15,rms) ! prevents rms=Inf

         if (rms.gt.9*rwas) then  ! Fit got much worse,
            if(debug) write(6,*)'Retracing step',rms,rwas
            call vmul(dx,1,-0.9,0,dx,1,nfp)
         else
            if(debug) write(6,*)'Continuing',rms,rwas
            thresh=(64*kconv-63)*(rms+0.01*abs(cx(n1)))/100000
            if(abs(rms-rwas).lt.thresh .or. nit+2*kconv-1.gt.mit) then
               kconv=kconv-1
               rwas=big
               if(kconv.ge.1) then
                  if(cf) call subtract_cf(obsrvd,calcul,resids,nmp,mmp)
                  call vsub(obsrvd,1,calcul,1,resids,1,nmp) ! residuals
               endif
            endif ! abs(rms-rwas).lt.thresh .or. nit+2*kconv-1.gt.mit

            ynoise=2.5*cx(n1)*corrld/sngl(0.1d0+snr)
c   The last few (i > NMP) elements of RESID & PD contains A PRIORI information
c   RESIDS(nmp+i) holds the values of (AX(i)-CX(i))*RAE(i)*YNOISE
c   while  PD(nmp+i,i) holds RAE(i)*YNOISE, the other elements being zero.
            call vdiv(ynoise,0,apu,1,pd(nmp+1),nmpfp+1,nfp)
            call vsub(apx,1,cx,1,wk,1,nfp)
            call vmul(wk,1,pd(nmp+1),nmpfp+1,resids(nmp+1),1,nfp)
            if(debug) then
               write(*,122)'apx=',(apx(j),j=n1,n4),(apx(j),j=1,ntg)
               write(*,122)'cx=',(cx(j),j=n1,n4),(cx(j),j=1,ntg)
122            format(a3,7f11.5)
            endif
c  Solve matrix equation PD.dx=resids
            call vmov(zero,0,wk,1,nfp)
c         write(*,*)'calling shfti...',nmpfp,nfp
            call shfti(pd,nmpfp,nmpfp,nfp,resids,nmpfp,
     &      1,tau,krank,rdum,wk,ip)
c            write(*,*)'called shfti...',(ip(j),j=1,nfp)

            call vmov(resids,1,dx,1,nfp)
            if(debug) then
               write(*,122)'dx=',(dx(jtg),jtg=n1,n4),(dx(jtg),jtg=1,ntg)
               if(krank.lt.nfp) then
                  write(6,*)'Rank Deficient:',krank,'  /',nfp
               else
                  write(6,*)'Full Rank:',krank,'  /',nfp
               endif
            endif
            if(dx(n1).gt.9.21) then
               fs=10000.
            elseif(dx(n1).lt.-9.21) then
               fs=.0001
            else
               fs=exp(dx(n1))
            endif
            dx(n1)=cx(n1)*(fs-1)
         endif  ! (rms.gt.9*rwas)

         if(kconv.eq.0 .or. mit.eq.0) go to 63   ! convergence
c  Limit maximum step size for fitted gases (non-linear).
c       do jtg=1,ntg
c       dxlimit=0.2+abs(cx(jtg))
c       if(abs(dx(jtg)) .gt. dxlimit) dx(jtg)=sign(dxlimit,dx(jtg))
c       end do
c       write(*,*)'nit,dx=',nit,dx
c       call vadd(cx,1,dx,1,cx,1,nfp)
         fr=1.0
         do jfp=1,nfp
            dxlimit=0.2+abs(cx(jfp))
            xfr=abs(dxlimit/abs(dx(jfp)))
            if(xfr.lt.fr) then
               fr=xfr
               kfp=jfp
            endif
         end do
         if(debug .and. fr.lt.1.)
     &   write(*,*)'Limited step size: kfp,fr=',kfp,fr
         call vsma(dx,1,fr,cx,1,cx,1,nfp)
         rwas=rms
      end do

c  Compute upper-triangle of covariance matrix & move diagonal elements into EX
63    var=(rms*corrld)**2/(1.-float(nfp)/(float(nmp)+0.1))
      var=var+tau**2  ! fudge to prevent EX=0 when rms=0
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
      ex(n1)=ex(n1)*cx(n1)**2
c      write(6,*)(cx(j),j=1,nfp)
c      write(6,*)(dx(j),j=1,nfp)
c      write(6,*)(sqrt(ex(j)),j=1,nfp)
c=========================================================================

c  Compute transmittances (convolved with ILS) of individual target gases.
c  Place them in PD, which just so happens to be exactly the right size.
      jva=1
      jpd=1
      kn2=1+ncp*(n1) ! Start address of workspace
      sh=rdec*(cx(n3)+sssss)
      do jtg=0,ntg
         if(jtg.eq.0) then ! non-target gases
            call vexp(spxv(jva),1,spxv(kn2),1,ncp)
         else        ! target gases
            call vmul(spxv(jva),1,cx(jtg),0,spxv(kn2),1,ncp)
            call vexp(spxv(kn2),1,spxv(kn2),1,ncp)
         endif
         call newdec(spxv(kn2),ncp,slit,nii,ldec,rdec,sh,pd(jpd),nmp)
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
c      call vadd (ex,1,(3*rms/cx(n1))**2,0,ex,1,nfp)   !  This is a fudge
c      call vsqrt(ex,1,ex,1,nfp)
      if(abs(cx(n1)).lt.tiny) cx(n1)=tiny
      do jfp=1,nfp
         ex(jfp)=sqrt(abs(ex(jfp))
c     &    +dx(jfp)**2      ! perturbation from adding ERROFF to residuals
     &   +(3*rms/cx(n1))**2   ! fudge
c     &    +(100*(cx(jfp)-apx(jfp))*rms/cx(n1))**2   ! fudge
c     &    +25*dxwas(jfp)**2   ! fudge
     &   +25*dx(jfp)**2   ! fudge
c     &   +0.5*esat(jfp)**2   ! fudge to increase error for saturated target lines
c     &    +0.04*(cx(n1)**2+1./cx(n1)**2)*(cx(jfp)-apx(jfp))**4  ! fudge
     &   )
      end do
      return
      end
