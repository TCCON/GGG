c  texo.f
c  Program to derive an exo-atmospheric spectrum by extrapolating
c  a series of atmospheric spectra to zero airmass.
c
c Minimizes  SUMk SUMi [Y(i,k) - F(i,k)]^2
c
c where
c  Y(i,k) measured signal at the i'th spectral point of the k'th spectrum
c  F(i,k) is a multivariate function simulating the calculated radiances
c
c  F(i,k) = F(Dk,Ak,Ei,Ci) = ILS # (Ei(Dk)*exp(Ci*Ak))
c
c where
c  Ei   exo-atmospheric spectrum (function of the Doppler stretch)
c  Ci   -ve atmospheric absorption coefficient spectrum (per unit airmass)
c  Dk   Doppler stretch of the solar features in the k'th spectrum
c  Ak   airmass of the k'th spectrum.
c  #    convolution operator.
c
c Assumes:
c  Ak is known. The unknowns are therefore Ei, Ci, and Dk
c
c  Assumes that the shifted (solar) features are broader
c  than their shifts, so that the partial differentials
c  wrt shift can be computed numerically.
c
c  That each spectrum can be decomposed into the product of two components:
c  1) An airmass-independent solar/instrumental spectrum (Ei)
c  2) An airmass-dependent absorptance spectrum (Ci)
c  The program further assumes that component (1) is subject to a
c  spectrum-dependent doppler stretch whereas component (2) is fixed
c  in frequency.
c
c  The main limitation of this approach is when the spectra contain an
c  airmass-independent component which does NOT shift in frequency with
c  the solar lines. The can arise due to instrumental H2O in the case of
c  MkIV balloon spectra or cell-HCl in the case of TCCON measurements.
c  This causes the doppler stretches to be underestimated, which causes
c  the derived solar features to have an anti-symmetric airmass-dependent
c  component, which causes their position to appear to be airmass-dependent.

c  MATH
c  Ni is the number of points in each spectrum
c  Nk is the number of spectra in the occultation.
c
c  There are 2*Ni+Nk unknowns and Ni*Nk= observations.
c  So the problem is well conditioned for Nk>=2 provided
c  that we ensure that S1, D1, A1, and ILS are all known.
c
c  For MkIV spectra having Ni ~5x10^5  and Nk ~40, the matrix
c  of PD's is of dimension (Nk.Ni, 2*Ni) = 2x10^7 x 10^6 and
c  so the time taken to solve this problem (~M.N^2)
c  is ~2x10^19 FLOPS per iteration.  So it is completely
c  impractical to do a straight-forward simultaneous fit
c  of the million unknowns to the twenty-million observations.
c
c  Fortunately, there is a high degree of similarity
c  between different rows/columns of the PD matrix, i.e.
c  df/dDk looks very similar from one spectrum (k) to the next.
c  Aditionally, there is very little correlation between the
c  neighboring values of Ei and Ci, suggesting that it may be
c  feasible to solve for Ei and Ci one i at a time.
c
c  So the plan is to solve for Ei and Ci simultaneously,
c  but separately for each i, So at each step we'll have
c  2 unknowns and Nk measurements, so that will take 4.Nk
c  FLOPS per i-value. So the total will be 4.Nk.Ni FLOPS.
c
c  The we'll determine Dk independently for each spectrum.
c  So at each step we'll have 1 unknown and Ni measurements,
c  so this will take  9.Ni FLOPS per k-value,  or 9.Nk.Ni total.
c
c  So the grand total is about 13.Nk.Ni = 10^8 FLOPS per iteration,
c  as compared with 2x10^19 for the fully-simultaneous method.
c  Of course, the short-cut method will take longer to converge
c  than the full method, but as long as it takes < 10^11
c  iterations we'll still come out ahead.
c
c  The partial differentials are calculated analytically:
c     dF/dDk = SUMi [Sk * ILS # (dEi/dDk * exp(Ci*Ak))]
c
c     dF/dEi = SUMk [Sk * ILS # (exp(Ci*Ak))]
c     dF/dCi = SUMk [Sk * ILS # (Ak*Ei(Dk)*exp(Ci*Ak))]
c  where the # symbol represents the convolution operator.
c
c
c  Program Summary:
cc  Main iteration loop
c      do iter=1,mit
c
cc     Update exo-atmospheric and Absorption coeff spectra.
c         do i=1,nmp  ! One spectral point at a time
c             compute dFk/dEi & dFk/dCi k=1,nspe
c             compute Rik, k=1,nspe
c             solve PD.x=R  ! PD(NSPE,2)
c             Ei=Ei+dEi
c             Ci=Ci+dCi
c         end do
c
cc      Update Dk values
c         do kspe=1,nspe  ! One spectrum at a time
c            compute E(Dk), and exp(Ci*Ak), i=1,nmp
c            compute dFi/dDk,  i=1,nspe
c            compute Rik  i=1,nmp
c            solve PD.x=R  ! PD(NMP,1)
c            Dk=Dk+dDk
c         end do
c
c      end do  ! iter=1,mit

c  Program also has a self-test capability, in which it creates synthetic
c  spectra over a user-selected spectral region with which to test the texo program.
c  Then tries to fit them.
c
      implicit none
      integer*4 lunb,lunr,lunw,luno,lnbc,lloc,llfs,nii,iseed,
     & mmp,nmp,mspe,nspe,possp,bytepw,bytepw1,kspe,nfp,iabpw,
     & iyr,iset,i,j,k,mit,iter,kconv,mconv,iend,iline,nk
      parameter (lunr=14,lunw=15,luno=16, lunb=17, nk=4,
     &  mspe=99, mmp=995000,mit=15,nfp=2, mconv=7, nii=21)
      real*4 stretch(mspe),st,airmass(mspe),am,exo,
     & y(mmp,mspe),r4buf(mmp),rms,rtry,apwt,smin,var,gasdev,
     & twt,tj,tjj,ty,tjy,otee,otcc,de_err,dc_err,
     & coswt(nk),slit1(nii),slit2(nii),wt(mmp),tot,
     & y2,ty2,huge,rbest,rbestwas,lmp,rmsarr(mspe),f2,
     & ds,dpd(mmp),dd,ddlo,ddhi,amin,amax,denom,
     & tcc(mmp),tee(mmp),tec(mmp),tre(mmp),trc(mmp),
     & cbest(mmp),ebest(mmp),ctry(mmp),etry(mmp),dc,de,
     & f(mmp),se(mmp),res(mmp),cres(mmp,mspe),eres(mmp,mspe),
     & cpd(mmp,mspe),epd(mmp,mspe),absor,dpdc,dpd2

      real*8 resnog
      integer*4 kmin,kmax,ifirst,ilast,istat,lr,m1,m2,kfmin,kfmax,lrt,lo
      integer*2 i2buf(2*mmp)
      byte bbuf(4*mmp)
      equivalence(bbuf,i2buf,r4buf)

      real*8 pi,graw,freq,v(3),tdd,trd,fmin,fmax,fbar,fwid
      real*8
     & oblat,           ! observation latitude (deg).
     & oblon,           ! observation longitude (deg).
     & obalt,           ! observation altitude (km)
     & zpdtim,          ! Time of ZPD (UT hours)
     & asza,            ! astronomical solar zenith angle (unrefracted)
     & zenoff,          ! zenith pointing offset
     & fovi,            ! Internal angular diameter of FOV (radians)
     & fovo,            ! External angular diameter of FOV (radians)
     & amal,            ! angular misalignment of interferometer (radians)
     & zoff,            ! zero-level offset
     & snr,             ! Signal-to-Noise Ratio
     & opd,             ! Optical path difference (cm) of interferogram
     & tins,            ! Temperature INSide the instrument
     & pins,            ! Pressure INSide the instrument
     & hins,            ! Humidity INSide the instrument
     & tout,            ! Temperature OUTside the instrument
     & pout,            ! Pressure OUTside the instrument
     & hout,            ! Humidity OUTside the instrument
     & sia,             ! Solar Intensity (Average)
     & sis,             ! Solar Intensity (SD)
     & aipl,            ! Airmass-Independent Path Length (km)
     & lasf,            ! laser frequency (e.g. 15798.03 cm-1)
     & wavtkr           ! suntracker operating frequency (e.g. 9900 cm-1)

      character runlab*21, rlmin*21, rlmax*21, version*40,
     & path*128, dplist*40, cc*1, outpath*40,rlabs*21,rlext*21,rlexo*21,
     & runlog*80, col1*1, apf*2, rlheader*400,root*64, header*40
      logical self_test
      version='TEXO   Version 1.4.1   2008-02-21   GCT '
c
      runlab='                     '
      rlmin= '                     '
      rlmax= '                     '
      rlabs= '                     '
      rlext= '                     '
      rlexo= '                     '
      pi=4*datan(1.0d0)
      self_test=.true.
      self_test=.false.
      call getenv('GGGPATH',root)
      lrt=lnbc(root)
      dplist=root(:lrt)//'/config/data_part.lst'
      outpath=root(:lrt)//'/spectra/'
      lo=lnbc(outpath)
      call getendian(iend)
      ds=2.0e-07
      apwt=1.e+03
      huge=9.9E+29
      ddlo=0.01
      ddhi=1.2
      smin=0
      iseed=4444
      
      write(*,*)version
      write(*,*)'Enter path to input file/runlog:'
      read(*,'(a)') runlog
c      runlog='/home/gka/ggg/runlogs/gnd/paIn_040909s.grl'
      bytepw=4*iend
      open(lunr,file=runlog,status='old')
      read(lunr,'(a)') rlheader

      write(*,*)'Enter Starting & Ending frequencies:'
      write(*,*)'Enter 0 99999 to retain original spectral limits'
      read(*,*) fmin,fmax
      fbar=0.5*(fmin+fmax)
      fwid=fmax-fmin
      do iline=1,3
         v(iline)=fmin+(iline+3)*fwid/10
      end do

      llfs=lloc(runlog,'/')
      open(lunw,file=runlog(:llfs)//'texo_'//runlog(llfs+1:),
     & status='unknown')
      write(lunw,'(a)') rlheader(:lnbc(rlheader))
c
c  Loop over spectra 
      amax=0.0
      amin=huge
      do kspe=1,mspe
11       call read_runlog(lunr,col1,runlab,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     &   ilast,graw,possp,bytepw1,zoff,snr,apf,tins,pins,hins,
     &   tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
c         if(istat.ne.0) write(*,*) 'Error in read_runlog:  istat=',istat
         if(istat.ne.0) go to 99
         iabpw=iabs(bytepw1)

c Weed out the spectra that don't fill our requested interval.
c Redefine IFIRST, ILAST, POSSP to make it appear that the
c measured spectrum extends from FMIN to FMAX.
         m1=nint(fmin/graw)
         m2=nint(fmax/graw)
         write(*,'(4i10,2x,a)') ifirst,ilast,m1,m2,runlab
         if((m1-ifirst)*float(m2-ilast).gt.0.0) go to 11
         kfmin=m1
         if(m1.lt.ifirst) kfmin=ifirst
         kfmax=m2
         if(m2.gt.ilast) kfmax=ilast
         possp=possp+iabpw*(kfmin-ifirst)
         nmp=kfmax-kfmin+1
         if(nmp.le.0)  go to 11  ! wrong detector spectral region
         if(nmp.gt.mmp) then
            write(*,*)kfmin,kfmax,nmp,mmp
            stop 'nmp > mmp'
         endif

         airmass(kspe)=pout/1013.25/cos(asza*pi/180.)
         if(airmass(kspe).gt.amax) then
            amax=airmass(kspe)
            kmax=kspe
            rlmax=runlab
         endif
         if(airmass(kspe).lt.amin) then
            amin=airmass(kspe)
            kmin=kspe
            rlmin=runlab
         endif

      if(self_test .eqv. .true.) then      ! compute synthetic spectra
         stretch(kspe)=-(kspe-1)*1.0E-6
         st=-stretch(kspe)  !  +ve am; -ve pm
         do i=1,nmp
           absor=0.0
           freq=(i-1+kfmin)*graw
           exo=20000.0*exp(-35*((freq*(1-st)-fbar)/fwid)**2)
           absor=absor+0.2*exp(-((v(1)-freq+v(1)*st)/0.3)**2)
           absor=absor+0.2*exp(-((v(3)-freq+v(3)*st)/0.3)**2)
           absor=absor+0.02*airmass(kspe)/(0.03**2+(v(2)-freq)**2)
           y(i,kspe)=exo*exp(-absor)+ 10.*gasdev(iseed)
           iseed=iseed+1113
c            write(*,*)kspe,freq,y(i,kspe),absor
         end do  ! i=1,nmp
         call write_binary_file(19,outpath(:lo)//runlab,header,0,bytepw,
     &   y(1,kspe),nmp)
         possp=0
      else                                 !  Read real measured spectra
         call gindfile(dplist,runlab,path)
         write(*,*)runlab,path
         if(lnbc(path).eq.0) then
            write(*,*) runlab,' not found on disk'
         endif
         open(lunb,file=path,access='direct',status='old',
     &    form='unformatted',recl=possp+iabpw*nmp)
         read(lunb,rec=1) (cc,j=1,possp),(bbuf(i),i=1,iabpw*nmp)
         close(lunb)
         if(iend*bytepw1.lt.0) call rbyte(bbuf,iabpw,nmp)
         if(iabpw.eq.2) then
            do i=1,nmp
               y(i,kspe)=i2buf(i)
            end do    !  i=1,nmp
         elseif(iabpw.eq.4) then
            do i=1,nmp
               y(i,kspe)=r4buf(i)
            end do    !  i=1,nmp
         else
            write(*,*)'Unsupported format'
         endif
      endif                         !  self_test .eqv. .true.

c  Write possibly-truncated spectrum and runlog
      call write_runlog(lunw,col1,runlab,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,kfmin,
     & kfmax,graw,possp,bytepw1,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)

      stretch(kspe)=-(1-kspe)*0.0E-6     !  Initial guess
      write(*,'(i5,2f13.7)') kspe,stretch(kspe),airmass(kspe)

      end do       !  kspe=1,mspe
      call read_runlog(lunr,col1,runlab,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     & ilast,graw,possp,bytepw1,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
      if(istat.ne.0) go to 99
      write(*,*)'nspe,mspe=',kspe,mspe
      stop ' nspe.ge.mspe'
99    close(lunr)
      nspe=kspe-1
      write(*,*)'nspe,nmp,kmin,kmax',nspe,nmp,kmin,kmax
      if(nspe.le.0) stop 'nspe <=0'
      lr=lnbc(rlmin)

c Initialize arrays Ei and Ci using highest and lowest airmass spectra.
c This is simply to get a good starting guess for Ei and Ci.
      write(*,*) 'Initializing arrays Ei and Ci '
      do j=1,nmp
         if(y(j,kmin).le.0) then
            dd=ddlo
         else
            dd=y(j,kmax)/y(j,kmin)
            if(dd.lt.ddlo) dd=ddlo
            if(dd.gt.ddhi) dd=ddhi
         endif
         cbest(j)=log(dd)/(airmass(kmax)-airmass(kmin))
         ebest(j)=y(j,kmin)*exp(-(airmass(kmin)*cbest(j)))
      end do  ! j=1,nmp

c  Main Fitting loop starts here
      kconv=0
      rbestwas=huge
      lmp=1.0  ! Levenberg-Marquardt Parameter.
      do while (kconv.lt.mconv)
      rbest=huge
      call vmov(ebest,1,etry,1,nmp)
      call vmov(cbest,1,ctry,1,nmp)
      lmp=sqrt(lmp)
      do iter=1,mit  ! Iteration Loop
c  Compute spectrum from c() and e() using resample subroutine
         rtry=0.0
         ty2=0.0
         do kspe=1,nspe
            st=stretch(kspe)
            am=airmass(kspe)
            call resample(etry,st,kfmin,se,nmp)         ! Ei' = E(Dk)
            call vmul(ctry,1,am,0,f,1,nmp)              ! Ak.Ci
            call vexp(f,1,f,1,nmp)                      ! exp(Ak.Ci)
c            call vmul(f,1,yscale(kspe),0,f,1,nmp)      ! exp(Ak.Ci)
            call resample(f,-st,kfmin,epd(1,kspe),nmp)    ! exp(Ak.Ci')
            call vmul(se,1,f,1,f,1,nmp)                 ! f=Ei'.exp()
            call vmul(f,1,am,0,cpd(1,kspe),1,nmp)       ! Ak.Ei'.exp()
            call vdot(f,1,f,1,y2,nmp)
            call vsub(y(1,kspe),1,f,1,res,1,nmp)        ! Res = y - f
            call resample(res,-st,kfmin,eres(1,kspe),nmp)    
            call vdot(res,1,res,1,rms,nmp)
            call vmov(res,1,cres(1,kspe),1,nmp)
            rtry=rtry+rms
            ty2=ty2+y2
         end do  ! kspe=1,nspe
         rtry=sqrt(rtry/ty2)
         write(*,'(2i5,3f13.7)')kconv,iter,lmp,rtry,rbest
         if(rtry.ge.rbest) then
c    Fit didn't improve. Keep previous estimates.
            if(lmp.gt.500.0) go to 134
            lmp=lmp*4
         else
c    Fit got better.
            lmp=lmp/2
            rbest=rtry
            call vmov(etry,1,ebest,1,nmp)
            call vmov(ctry,1,cbest,1,nmp)

c           Update  coefficients tcc,tee,tec,trc,tre
            do i=1,nmp
               tee(i)=0.0
               tcc(i)=apwt
               tec(i)=0.0
               tre(i)=0.0
               trc(i)=-cbest(i)*apwt
               do kspe=1,nspe
                  tcc(i)=tcc(i)+cpd(i,kspe)**2
                  tee(i)=tee(i)+epd(i,kspe)**2
                  tec(i)=tec(i)+cpd(i,kspe)*epd(i,kspe)
                  tre(i)=tre(i)+epd(i,kspe)*eres(i,kspe)
                  trc(i)=trc(i)+cpd(i,kspe)*cres(i,kspe)
               end do
            end do  ! i=1,nmp
         endif !  if(rtry.ge.rbest) then

         do i=1,nmp
            otee=(1+lmp)*tee(i)
            otcc=(1+lmp)*tcc(i)
            denom=(otee*otcc-tec(i)**2)
            de=(tre(i)*otcc-trc(i)*tec(i))/denom
            dc=(trc(i)*otee-tre(i)*tec(i))/denom
            de_err=sqrt(tee(i)/denom)
            dc_err=sqrt(tcc(i)/denom)
            etry(i)=ebest(i)+de
            ctry(i)=cbest(i)+dc
         end do  ! i=1,nmp
      end do  ! iter=1,mit
134   continue
      write(*,*)'rbest,rwas=',rbest,rbestwas
      if(rbest.ge.rbestwas) then
         kconv=kconv+1
      else
         kconv=0  ! fit improved
         rbestwas=rbest
      endif

      write(*,*)'Calculating Sk and Dk (scale & shift)'
c  Iterate on Sk and Dk
Compute spectrum from c() and e() using resample subroutine
      do kspe=1,nspe
         st=stretch(kspe)
         am=airmass(kspe)
         call resample(ebest,st+ds/2,kfmin,dpd,nmp) ! Ei" = E(Dk)+ds/2
         call resample(ebest,st-ds/2,kfmin,se,nmp)  ! Ei" = E(Dk)-ds/2
         call vsub(dpd,1,se,1,dpd,1,nmp)
         call vdiv(dpd,1,ds,0,dpd,1,nmp)             ! dE/dDk
         call resample(ebest,st,kfmin,se,nmp)       ! Ei' = E(Dk)
         call vmul(cbest,1,am,0,f,1,nmp)             ! Ak.Ci
         call vexp(f,1,f,1,nmp)                      ! exp(Ak.Ci)
         call vmul(f,1,dpd,1,dpd,1,nmp)              ! dE/dDk.exp(Ak.Ci)
         call vmul(se,1,f,1,f,1,nmp)                 ! f=Ei'.exp()
         call vdot(f,1,f,1,y2,nmp)
         call vsub(y(1,kspe),1,f,1,res,1,nmp)        ! Res = y - f
         call vdot(res,1,res,1,rms,nmp)
         call vdot(f,1,f,1,f2,nmp)
         rmsarr(kspe)=sqrt(rms/f2)

c        Update  coefficients tcc,tee,tec,trc,tre
         tdd=0.0d0
         trd=0.0d0
         do i=1,nmp
            tdd=tdd+dpd(i)**2
            trd=trd+dpd(i)*res(i)
         end do  ! i=1,nmp
         if(kspe.eq.kmin) smin=trd/tdd
         stretch(kspe)=stretch(kspe)+trd/tdd
      end do  ! kspe=1,nspe

      call resample(ebest,+ds/2,kfmin,dpd,nmp) ! Ei" = +ds/2
      call resample(ebest,-ds/2,kfmin,se,nmp)  ! Ei" = -ds/2
      call vsub(dpd,1,se,1,dpd,1,nmp)
      call vdiv(dpd,1,ds,0,dpd,1,nmp)           ! dE/dDk
      call vdot(dpd,1,dpd,1,dpd2,nmp)
      call vmul(dpd,1,ebest,1,dpd,1,nmp)
      call vdot(dpd,1,cbest,1,dpdc,nmp)
c Reference D and Y to the lowest aimass spectrum.
      do  kspe=1,nspe
         stretch(kspe)=stretch(kspe)-smin
     &   +(airmass(kspe)-airmass(kmin))*dpdc/dpd2
      end do   ! kspe=1,nspe
      end do   ! while(kconv.lt.mconv) 
      close(lunr)

c  Write out results
      open(lunr,file=runlog,status='old')
      read(lunr,'(a)') rlheader
      do kspe=1,mspe
10       call read_runlog(lunr,col1,runlab,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     &   ilast,graw,possp,bytepw1,zoff,snr,apf,tins,pins,hins,
     &   tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         if(istat.ne.0) exit
         iabpw=iabs(bytepw1)

         m1=nint(fmin/graw)
         m2=nint(fmax/graw)
c         write(*,'(4i10,2x,a)') ifirst,ilast,m1,m2,runlab
         if((m1-ifirst)*float(m2-ilast).gt.0.0) go to 10
         kfmin=m1
         if(m1.lt.ifirst) kfmin=ifirst
         kfmax=m2
         if(m2.gt.ilast) kfmax=ilast
         possp=possp+iabpw*(kfmin-ifirst)
         nmp=kfmax-kfmin+1
         if(nmp.le.0)  go to 10  ! wrong detector spectral region
         write(*,'(i3,a22,3f12.6)')kspe,runlab,
     &   airmass(kspe),stretch(kspe),rmsarr(kspe)
      end do
c
c  Write absorption coefficient spectrum (for diagnostic purposes)
      col1=' '
      call write_binary_file(19,outpath(:lo)//rlmin(:lr-3)//'abs',
     &header,0,-4,cbest,nmp)
      rlabs=rlmin(:lr-3)//'abs'
      call write_runlog(lunw,col1,rlabs,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,kfmin,
     &   kfmax,graw,0,-4,zoff,snr,apf,tins,pins,hins,
     &   tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)

c-------------------------------------------------------------------
c  Perform shaving of smoothed exo-atmospheric spectrum
c  Precompute cosine weighting operator
      do k=1,nk
         coswt(k)=cos(k*pi/(2*nk+2))**2
      end do
c
c  Apodize.
      resnog=0.5d0/opd/graw
      call profzl(1,nii,resnog,0.0d0,0.0d0,slit1)
      call vdot(slit1,1,1.0,0,tot,nii)
      call vmul(slit1,1,1.0/tot,0,slit1,1,nii)

      call profzl(2,nii,resnog,0.0d0,0.0d0,slit2)
      call vdot(slit2,1,1.0,0,tot,nii)
      call vmul(slit2,1,1.0/tot,0,slit2,1,nii)

      call vmov(ebest,1,etry,1,nmp)
      call vmov(cbest,1,ctry,1,nmp)
      do i=1,nmp-nii
        call vdot(etry(i),1,slit1,1,ebest(i+nii/2),nii)
        call vdot(ctry(i),1,slit2,1,cbest(i+nii/2),nii)
      end do

c  Pre-calculate weights
      var=0.5
      do i=1,nmp
         wt(i)=var/(var+cbest(i)**2)
      end do
c
c  Perform weighted SLF to points from i-nk to i+nk
      do i=1+nk,nmp-nk
        twt=99*nk*wt(i)**2
        tj=0.0
        tjj=0.0
        ty=99*nk*ebest(i)*wt(i)**2
        tjy=0.0
        do j=1,nk
          twt=twt+(wt(i-j)+wt(i+j))*coswt(j)
          tj=tj+(-wt(i-j)+wt(i+j))*j*coswt(j)
          tjj=tjj+(wt(i-j)+wt(i+j))*j*j*coswt(j)
          ty=ty+(ebest(i-j)*wt(i-j)+ebest(i+j)*wt(i+j))*coswt(j)
          tjy=tjy+((-ebest(i-j))*wt(i-j)+ebest(i+j)*wt(i+j))*j*coswt(j)
        end do
        denom=twt*tjj-tj*tj
        if(denom.eq.0.0) then
           write(*,*) 'denom=0.0'
           denom=0.001
        endif
        etry(i)=(ty*tjj-tjy*tj)/denom
      end do

c  Write exo-atmospheric spectrum (for diagnostic purposes)
      call write_binary_file(19,outpath(:lo)//rlmin(:lr-3)//'ext',
     & header,0,bytepw,ebest,nmp)
      rlext=rlmin(:lr-3)//'ext'
      call write_runlog(lunw,col1,rlext,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,kfmin,
     & kfmax,graw,0,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)

c  Write shaved EXO file and its runlog entry
      call write_binary_file(19,outpath(:lo)//rlmin(:lr-3)//'exo',
     & header,0,bytepw,etry,nmp)
      rlexo=rlmin(:lr-3)//'exo'
      call write_runlog(lunw,col1,rlexo,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,kfmin,
     & kfmax,graw,0,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)

      close(lunw)
      close(lunr)
      close(luno)
      stop
      end
