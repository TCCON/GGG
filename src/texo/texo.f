c  texo.f
c  Program to derive an exo-atmospheric spectrum by extrapolating
c  a series of atmospheric spectra to zero airmass.
c
c  Minimizes  SUMk SUMi [ Y(i,k) - Sk*Ei(Dk)*exp(Ak*Ci) ]^2
c
c  where
c   Y(i,k) measured signal at i'th spectral point of k'th spectrum
c   Ei   exo-atmospheric spectrum (function of the Doppler stretch)
c   Dk   Doppler stretch of the solar features in the k'th spectrum
c   Ci   atmospheric absorption coefficient spectrum (per unit airmass)
c   Ak   airmass of the k'th spectrum.
c
c  Assumes:
c   Ak and Dk are known.
c   Ei, Ci  are unknown.
c   Each spectrum can be decomposed into the product of two components:
c   1) An airmass-independent solar/instrumental spectrum (Ei)
c   2) An airmass-dependent absorptance spectrum (Ci)
c  The program further assumes that component (1) is subject to a
c  spectrum-dependent doppler stretch whereas component (2) is fixed
c  in frequency.
c
c  A limitation of the airmass extrapolation is that lines may be
c  under-resolved and saturated and therefore growing non-linearly
c  with airmass (e.g. CO2 at 2330, 3600 and 3700 cm-1. Thus a linear
c  extrapolation to zero airmass won't completely fully remove these
c  types of lines. Hence the need for a subsequent "shaving" procedure. 
c  Even worse, this residual airmass-independent absorption arising
c  from saturated atmospheric lines, plus wavenumber-invariant
c  instrumental absorption features (e.g. H2O, cell HCl, channeling),
c  will not be Doppler stretched like the solar features and will
c  therefore cause the stretch to be under-estimated.
c
c  The Dk values are computed from the solar shifts, based on the
c  runlog OSDS (Observer-Sun Doppler Stretch) values. If the Dk values
c  were assumed to be zero, the solar features would exhibit anti-
c  symmetric residuals.  Since we are usually more interested in
c  spectral regions containing solar lines, rather than strongly
c  saturated CO2 lines, it is better to choose the Dk values that
c  minimize artifacts from the solar features.
c
c  Another problem is that if the balloon gets to float around noon,
c  the solar shifts and the airmass both increase monotonically in
c  the subsequent spectra. This means that the solar lines can be
c  fitted just as well by anti-symmetric absorption artifacts as by
c  a real frequency stretch.  This would not be true if there were
c  high sun spectra both before and after local noon, but this rarely
c  happens.  If frequency shifts are fitted, the derived absorption
c  spectrum (Ci) can have anti-symmetric artifacts at every solar
c  line, while the solar spectrum (Ei) itself looks fine.
c
c  A previous version of TEXO used to retrieve a Dk stretch value for
c  each spectrum based on mimimizing the fit to the measured spectra.
c  But the derived Dk values were found to be typically smaller than
c  the expected (OSDS) values and would vary considerably depending on
c  the number of solar versus instrumental features in the spectral
c  interval analyzed. So instead we use the computed OSDS values in
c  the runlog, although this assumes that the FOV is pointed at the
c  center of the sun for the whole flight.  If the pointing is offset,
c  this can cause a spurious doppler shift that varies between noon
c  and sunset as the solar spin axis rotates by 90 deg. Since the
c  sun is rotating at 1900 m/s, this has the potential to dwarf the
c  doppler stretch due to the Earth's rotation (380 m/s at 35N).
c
c  An inherent ambiguity of this approach occurs when a spectral
c  region has zero signal in all spectra, even in the lowest airmass
c  spectrum. This could arise because of a very strong absorption
c  (e.g. in the line centers of the nu3 CO2 band) or because the
c  instrument response is zero (at edges of the bandpass). By
c  looking at the adjacent spectal points this ambiguity can be
c  resolved, but since the airmass extra-polation analysis method
c  only sees one wavenumber at a time, it doesn't know how to handle
c  it.  Hence the need for "shaving".
c
c  We use a priori constraints to control this. In the center of the
c  bandpass we want to set the exo-spectrum to 20,000 +/- 10,000
c  whereas at the edges we want 0 +/- 10,000.
c
c  Jacobians are obtained by differentiating the cost function:
c    Chi_squared = SUMk SUMi [ Y(i,k) - Sk.Ei(Dk).Exp(Ak.Ci) ]^2
c
c    d/dCi     = SUMk [-Ak.Sk.Ei(Dk).Exp(Ak.Ci).RESi]
c    d/dEi(Dk) = SUMk [-Sk.Exp(Ak.Ci).RESi]
c    d/dSk     = SUMi [-Ei(Dk).Exp(Ak.Ci).RESi]
c    d/dDk     = SUMi [-Sk.dEi(Dk)/dDk.Exp(Ak.Ci).RESi]
c    d/dAk     = SUMi [-Sk.Ei(Dk).Ci.Exp(Ak.Ci).RESi]
c       where RESi = Y(i,k) - Sk.Ei(Dk).Exp(Ak.Ci)
c
c   A complication arises because the cost function contains Ei(Dk)
c   rather then Ei(0). The Jacobian wrt Ei can be expressed as the
c   Jacobian wrt Ei(Dk) stretched by -Dk such that
c       d/dEi = SUMk [-Sk.Exp(Ak.Ci(-Dk)).RESi(-Dk)]
c   A drawback with this formulation is that Ci, and hence RESi,
c   contains higher resolution structure than Ei, such that Ci(-Dk)
c   cannot be accurately approximated by a expansion of Ci(0).
c
c  At the minimum of cost function, the Jacobians are all zero and  
c     RESi' = Y(i,k) - Sk.(Ei+de)(Dk)*Exp(Ak.(Ci+dc))
c  where de and dc are the adjustments to Ei and Ci needed to
c reach the global minimum.
c   RESi' = Y(i,k) - Sk.(Ei+de)(Dk)*Exp(Ak.Ci+Ak.dc))
c   RESi' = Y(i,k) - Sk.(Ei+de)(Dk)*Exp(Ak.Ci).Exp(Ak.dc))
c   If Ak.dc <<1
c   RESi' = Y(i,k) - Sk.(Ei+de)(Dk)*Exp(Ak.Ci).(1-Ak.dc)
c   RESi' = RESi - de(Dk).Sk.Exp(Ak.Ci) + dc.Ak.Sk.Ei(Dk)*Exp(Ak.Ci))
c
c  MATH
c  Ni is the number of points in each spectrum
c  Nk is the number of spectra in the occultation.
c
c  There are 2*Ni+Nk unknowns and Ni*Nk= observations.
c  So the problem is well conditioned for Nk>=2 provided
c  that we ensure that S1, D1, A1 are all known.
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
c  Program Summary:
cc  Iteration loop
c   do iter=1,mit
c
cc  Update  i-dependent vectors (Ci and Ei)
c         do i=1,nmp  ! One spectral point at a time
c             compute dFk/dEi & dFk/dCi k=1,nspe
c             compute Rik, k=1,nspe
c             solve PD.x=R  ! PD(NSPE,2)
c             Ei=Ei+dEi
c             Ci=Ci+dCi
c         end do
c
cc   Update k-dependent vectors (Sk, Ak, Dk)
c         do kspe=1,nspe  ! One spectrum at a time
c            compute E(Dk), and exp(Ci*Ak), i=1,nmp
c            compute dFi/dDk,  i=1,nspe
c            compute Rik  i=1,nmp
c            solve PD.x=R  ! PD(NMP,1)
c            Dk=Dk+dDk
c         end do
c
c   end do  ! iter=1,mit

c  Program also has a self-test capability, in which it creates
c  synthetic spectra over a user-selected spectral region with which
c  to test the texo program. Then tries to fit them.
c
      implicit none
      include "../gfit/ggg_int_params.f"

      integer*4 lunr_bs,lunr_rlg,lunw_rlg,luno,lnbc,lloc,llfs,nii,iseed,
     & mmp,nmp,mspe,nspe,possp_rrl,possp_wrl,bytepw_rrl,bytepw_wrl,
     & kspe,iabpw,idum,imax,
     & iyr,iset,i,j,k,mit,iter,kconv,mconv,iend,iline,nkshav
      parameter (lunr_rlg=14,lunw_rlg=15,luno=16,lunr_bs=17,nkshav=6,
     &  mspe=124, mmp=695000,mit=12,mconv=2, nii=21)
      real*4 stretch(mspe),st,airmass(mspe),am,spi,
     & demax,emin,emax,fd,
     & e_apv,e_apu,
     & c_apv,c_apu,
c     & exo,f2,rms,
     & y(mmp,mspe),r4buf(mmp),rtry,apwt,eap,smin,var,gasdev,
     & twt,tj,tjj,ty,tjy,otee,otcc,de_err,dc_err,
     & coswt(nkshav),slit1(nii),slit2(nii),wt(mmp),tot,
     & y2,ty2,huge,rbest,rbestwas,lmp,rmsarr(mspe),
     & ds,dd,ddlo,ddhi,amin,amax,denom,cont(mmp),
     & tcc(mmp),tee(mmp),tec(mmp),tre(mmp),trc(mmp),
     & cbest(mmp),ebest(mmp),ctry(mmp),etry(mmp),dc,de,
     & f(mmp),se(mmp),res(mmp),cres(mmp,mspe),eres(mmp,mspe),
c     & dpd(mmp),dpdc,dpd2,
     & cpd(mmp,mspe),epd(mmp,mspe)
      parameter(spi=3.14159265359)

      real*8 resn,rect,resnog,w,absor
      integer*4 kmin,kmax,ifirst,ilast,istat,lr,m1,m2,kfmin,kfmax,lrt,lo
      integer*2 i2buf(2*mmp)
      byte bbuf(4*mmp)
      equivalence(bbuf,i2buf,r4buf)

      real*8 
c     & tdd,trd,
     & freq,v(3),fmin,fmax,fbar,fwid,effres,
     & graw,grawmin,    ! Spectral point spacing (cm-1)
     & oblat,oblatmin,  ! observation latitude (deg).
     & oblon,oblonmin,  ! observation longitude (deg).
     & obalt,obaltmin,  ! observation altitude (km)
     & zpdtim,zpdtimmin,! Time of ZPD (UT hours)
     & asza,            ! astronomical solar zenith angle (unrefracted)
     & zenoff,          ! zenith pointing offset
     & fovi,            ! Internal angular diameter of FOV (radians)
     & fovo,            ! External angular diameter of FOV (radians)
     & amal,            ! angular misalignment of interferometer (radians)
     & zoff,            ! zero-level offset
     & snr,             ! Signal-to-Noise Ratio
     & opd,opdmin,      ! Optical path difference (cm) of interferogram
     & tins,            ! Temperature INSide the instrument
     & pins,            ! Pressure INSide the instrument
     & hins,            ! Humidity INSide the instrument
     & tout,            ! Temperature OUTside the instrument
     & pout,            ! Pressure OUTside the instrument
     & hout,            ! Humidity OUTside the instrument
     & sia,             ! Solar Intensity (Average)
c    & sis,             ! Solar Intensity (SD)
     & fvsi,            ! Fractional Variation in Solar Intensity
     & wspd,            ! Wind Speed (m/s)
     & wdir,            ! Wind Direction (deg)
     & azim,azimmin,    ! Solar Azimuth Angle
     & osds,osdsmin,    ! Observer-Sun Doppler Stretch (ppm)
     & aipl,            ! Airmass-Independent Path Length (km)
     & lasf,            ! laser frequency (e.g. 15798.03 cm-1)
     & wavtkr           ! suntracker operating frequency (e.g. 9900 cm-1)

      character rlmin*57, rlmax*57,data_fmt_write_rl*256, 
     & path*128, dplist*90, cc*1, outpath*90,rlabs*57,rlext*57,rlexo*57,
     & col_labels_rl*320,
     & header*40

      character
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & sn(mspe)*(nchar),          !spectrum name vector
     & sptpath*128,
     & data_fmt_read_rl*256,      !
     & version*64,                !current program version
     & rlgfile*120,               !name of runlog file
     & fmin_str*80,fmax_str*80

      logical self_test


      data cont/mmp*20000.0/

      version=' TEXO   Version 1.62   2019-04-14   GCT '
c
      idum=mauxcol    ! avoid comlier warnings (Unused parameter)
      idum=mcolvav    ! avoid comlier warnings (Unused parameter)
      idum=mfilepath  ! avoid comlier warnings (Unused parameter)
      idum=mgas       ! avoid comlier warnings (Unused parameter)
      idum=mlev       ! avoid comlier warnings (Unused parameter)
      idum=mrow_qc    ! avoid comlier warnings (Unused parameter)
      idum=mspeci     ! avoid comlier warnings (Unused parameter)
      idum=mvmode     ! avoid comlier warnings (Unused parameter)
      idum=ncell      ! avoid comlier warnings (Unused parameter)

      specname='                     '
      rlmin= '                     '
      rlmax= '                     '
      rlabs= '                     '
      rlext= '                     '
      rlexo= '                     '
      self_test=.true.
      self_test=.false.
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)
      dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
      outpath=gggdir(:lrt)//'spectra'//dl
      lo=lnbc(outpath)
      call getendian(iend)
      bytepw_wrl=4*iend
      ds=2.0e-07
c      apwt=5.0
c      eap=10.0
      huge=9.9E+29
      ddlo=0.01
      ddhi=1.2
      smin=0
      iseed=4444
      possp_wrl=0
      
      write(*,*)version
      if (iargc() == 0) then
         write(*,*)'Enter path to input file/runlog:'
         read(*,'(a)') rlgfile
      elseif (iargc() == 3) then
         call getarg(1, rlgfile)
         call getarg(2, fmin_str)
         read(fmin_str,*)fmin
         call getarg(3, fmax_str)
         read(fmax_str,*)fmax
      else
         write(*,*) 'Usage: $gggpath/bin/texo path/runlog '//
     &  'start_freq end_freq (enter 0 99999 to retain spectral limits'
         stop
      endif

      open(lunr_rlg,file=rlgfile,status='old')
      call read_runlog_header(lunr_rlg,data_fmt_read_rl,col_labels_rl)
      write(*,*)'data_fmt_write_rl=',data_fmt_read_rl
c      read(lunr_rlg,*) nlhead,ncol
c      do i=2,nlhead
c      read(lunr_rlg,'(a)') col_labels_rl
c      end do

      if (iargc() == 0) then
         write(*,*)'Enter Starting & Ending frequencies:'
         write(*,*)'Enter 0 99999 to retain original spectral limits'
         read(*,*) fmin,fmax
      endif
      fbar=0.5*(fmin+fmax)
      fwid=fmax-fmin
      do iline=1,3
         v(iline)=fmin+iline*fwid/4
      end do

      llfs=lloc(rlgfile,'/')
      open(lunw_rlg,file=rlgfile(:llfs)//'texo_'//rlgfile(llfs+1:),
     & status='unknown')
c      write(lunw_rlg,*) 3,36
c      write(lunw_rlg,'(a)') 'Created by texo.f'
c      write(lunw_rlg,'(a)') col_labels_rl(:lnbc(col_labels_rl))
      call write_runlog_header(lunw_rlg,'Created by '//version,
     & data_fmt_write_rl)
c
      write(73,*)'   path      kspe       d_stretch       airmass(kspe)'

c  Loop over runlog entries, computing airmass & stretch, 
c  and reading in relevant spectra 
      amax=0.0
      amin=huge
      do kspe=1,mspe
11       call read_runlog_data_record(lunr_rlg,data_fmt_read_rl,col1,
     & specname,iyr,iset,zpdtim,oblat,oblon,obalt,asza,zenoff,azim,
     & osds,opd,fovi,fovo,amal,ifirst,ilast,
     & graw,possp_rrl,bytepw_rrl,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

c         if(istat.ne.0) write(*,*) 'Error in read_runlog:  istat=',istat
         if(istat.ne.0) go to 99
         if(col1.eq.':') go to 11
         iabpw=iabs(bytepw_rrl)
         sn(kspe)=specname

c Weed out the spectra that don't fill our requested interval.
c Redefine IFIRST, ILAST, POSSP to make it appear that the
c measured spectrum extends from FMIN to FMAX.
         m1=nint(fmin/graw)
         m2=nint(fmax/graw)
c         write(*,'(4i10,2x,a)') ifirst,ilast,m1,m2,specname
         if((m1-ifirst)*float(m2-ilast).gt.0.0) go to 11
         kfmin=m1
         if(m1.lt.ifirst) kfmin=ifirst
         kfmax=m2
         if(m2.gt.ilast) kfmax=ilast
         nmp=kfmax-kfmin+1
         if(nmp.le.0)  go to 11  ! wrong detector spectral region
         if(nmp.gt.mmp) then
            write(*,*)kfmin,kfmax,nmp,mmp
            stop 'nmp > mmp'
         endif

         stretch(kspe)=sngl(osds*1.E-06)
         stretch(kspe)=sngl(osds*0.E-06)
         airmass(kspe)=sngl(pout/1013.25/dcos(0.982*asza*spi/180.))
         if(airmass(kspe).gt.amax) then
            amax=airmass(kspe)
            kmax=kspe
            rlmax=specname
         endif
c  Find lowest airmass spectrum
         if(airmass(kspe).lt.amin) then
            amin=airmass(kspe)
            kmin=kspe
            rlmin=specname
            zpdtimmin=zpdtim
            oblatmin=oblat
            oblonmin=oblon
            azimmin=azim
            osdsmin=osds
            opdmin=opd
            grawmin=graw
         endif

         if(self_test .eqv. .true.) then   ! compute synthetic spectra
            st=stretch(kspe)  !  +ve am; -ve pm
            do i=1,nmp
               absor=0.0
               freq=(i-1+kfmin)*graw
c               exo=20000.0*exp(-35*((freq*(1-st)-fbar)/fwid)**2)
               absor=absor+.2*dexp(-((v(1)-freq-v(1)*st)/.1)**2) ! solar line
               absor=absor+.2*dexp(-((v(3)-freq-v(3)*st)/.1)**2) ! solar line
               absor=absor+.011*airmass(kspe)/(0.025**2+(v(2)-freq)**2)
               y(i,kspe)=20000.0*exp(-sngl(absor))+ 00.*gasdev(iseed)
               iseed=iseed+1113
c               write(*,*)kspe,freq,y(i,kspe),absor
            end do  ! i=1,nmp
         else                            !  Read real measured spectra
            call gindfile(dplist,specname,path)
            if(lnbc(path).eq.0) then
               write(*,*) specname,' not found on disk'
            endif
            open(lunr_bs,file=path,access='direct',status='old',
     &      form='unformatted',recl=possp_rrl+iabpw*(nmp+kfmin-ifirst))
            read(lunr_bs,rec=1) (cc,j=1,possp_rrl+iabpw*(kfmin-ifirst)),
     &      (bbuf(i),i=1,iabpw*nmp)
            close(lunr_bs)
            if(iend*bytepw_rrl.lt.0) call rbyte(bbuf,iabpw,nmp)
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

c  write binary spectrum
         call write_binary_file(19,outpath(:lo)//'texo_'//specname,
     &   header,0,bytepw_wrl,y(1,kspe),nmp)

c  Write runlog entry for possibly-truncated original spectrum
         call write_runlog_data_record(lunw_rlg,data_fmt_write_rl,
     &   col1,'texo_'//specname,iyr,iset,zpdtim,oblat,oblon,obalt,
     &   asza,zenoff,azim,osds,opd,fovi,fovo,amal,kfmin,kfmax,
     &   graw,possp_wrl,bytepw_wrl,zoff,snr,apf,tins,pins,hins,
     &   tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)


         write(73,'(a30,i5,2f13.7)') path,kspe,
     &   stretch(kspe)-stretch(kmin),airmass(kspe)

      end do       !  kspe=1,mspe
      write(*,*)'nspe,mspe=',kspe,mspe
      stop ' Increase parameter MSPE '
99    close(lunr_rlg)
      nspe=kspe-1
      write(*,*)'nspe,nmp,kmin,kmax',nspe,nmp,kmin,kmax
      if(nspe.le.0) stop 'nspe <=0'
      lr=lnbc(rlmin)

c Initialize arrays Ei and Ci using highest and lowest airmass spectra.
c This is simply to get a good starting guess for Ei and Ci.
      write(*,*) 'Initializing arrays Ei and Ci '
      emin=1e+36
      emax=-1e+36
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
         write(71,*)2,6
         write(71,*)'kconv iter kspe airmass y(10934) f(10934)'
      write(*,*)'                       iter          lmp        rtry
     &        rbest        ebest(10934)      cbest(10934)'
      kconv=0
      rbestwas=9999.
      lmp=1.0  ! Levenberg-Marquardt Parameter.
      do while (kconv.lt.mconv)
         rbest=9999.
         call vmov(ebest,1,etry,1,nmp)
         call vmov(cbest,1,ctry,1,nmp)
         lmp=sqrt(lmp)

         write(*,*)'kconv,mconv = ',kconv,mconv
         do iter=1,mit  ! Iteration Loop
c  Compute spectrum from c() and e() using resample subroutine
            rtry=0.0
            ty2=0.0
            do kspe=1,nspe
               st=stretch(kspe)-stretch(kmin)
               am=airmass(kspe)
               call resample(etry,st,kfmin,se,nmp)        ! Ei' = E(Dk)
               call vmul(ctry,1,am,0,f,1,nmp)             ! f=Ak.Ci
               call vexp(f,1,f,1,nmp)                     ! f=exp(Ak.Ci)
c               call vmul(f,1,yscale(kspe),0,f,1,nmp)     ! exp(Ak.Ci)
               call resample(f,-st,kfmin,epd(1,kspe),nmp) ! epd=exp(Ak.Ci')
               call vmul(se,1,f,1,f,1,nmp)                ! f=Ei'.exp()
               call vmul(f,1,am,0,cpd(1,kspe),1,nmp)      ! cpd=Ak.Ei'.exp()
               call vdot(f,1,f,1,y2,nmp)
               call vsub(y(1,kspe),1,f,1,res,1,nmp)       ! Res = y-f
               call resample(res,-st,kfmin,eres(1,kspe),nmp)    
               call vdot(res,1,res,1,var,nmp)
               call vmov(res,1,cres(1,kspe),1,nmp)
               rtry=rtry+var
               ty2=ty2+y2
               rmsarr(kspe)=sqrt(var/y2)
               if( kconv.eq.mconv-1 .and. lmp.gt.500.) then
c                  write(*,*) 'Writing SPT file ztexo_'//sn(kspe)
                  rect=0.5d0*(fmin+fmax)*(fovi**2+amal**2)/8
                  resn=0.6d0/opd
                  effres=dsqrt(rect**2+resn**2)
c                  write(*,'(a,3f8.4)')'effres,resn,rect = ',effres,resn,rect
                  sptpath='/home/toon/ddd/spt/ztexo_'//sn(kspe)
                  call write_spt(19,0,sptpath,y(1,kspe),f,cont,
     &            1.,1.,graw*kfmin,0.d0,effres,graw,1.,'texo test',0.d0,
     &            0.0d0,0.0d0,0.0,rmsarr(kspe),0.0,1.0,f,f,nmp,0,0)
               endif

            write(71,*)kconv,iter,kspe,airmass(kspe),
     &      y(10934,kspe),f(10934)
            end do  ! kspe=1,nspe

            rtry=sqrt(rtry/ty2)
c            write(*,'(2i5,3f13.7)')kconv,iter,lmp,rtry,rbest
            if(rtry.ge.rbest) then
c    Fit didn't improve. Keep previous estimates.
               write(*,*) 'Fit worsened',iter,lmp,rtry,rbest,
     &      ebest(10934),cbest(10934)
               if(lmp.gt.500.0) go to 134
               lmp=lmp*4
            else
c    Fit got better.
               write(*,*) 'Fit improved',iter,lmp,rtry,rbest,
     &      ebest(10934),cbest(10934)
               lmp=lmp/2
               rbest=rtry
               call vmov(etry,1,ebest,1,nmp)
               call vmov(ctry,1,cbest,1,nmp)
               
c      write(*,*)'2: emin,emax = ',vmin(ebest,1,nmp),vmax(ebest,1,nmp)

c  We want to choose the delta_Ci and delta_Ei at each i to minimize
c   SUMk [ RESi - delta_Ci.CPD(i,k) - delta_Ei.EPD(i,k) ]^2
c  Differentiating wrt delta_Ci and delta_Ei and setting to zero yields
c   SUMk [CPD(i,k).[ RESi - delta_Ci.CPD(i,k) - delta_Ei.EPD(i,k) ]] = 0
c   SUMk [EPD(i,k).[ RESi - delta_Ci.CPD(i,k) - delta_Ei.EPD(i,k) ]] = 0
c
c  So we have the following matrix equation to be solved at each i:
c  |SUMk CPD(i,k)^2  SUMk CPD(i,k).EPD(i,k)| |delta_Ci| = |SUMk RESi.CPD(i,k)|
c  |SUMk EPD(i,k).CPD(i,k)  SUMk EPD(i,k)^2| |delta_Ei| = |SUMk RESi.EPD(i,k)|
c
c  Let 
c    tcc = SUMk CPD(i,k)^2
c    tee = SUMk EPD(i,k)^2
c    tec = SUMk EPD(i,k).CPD(i,k)
c    trc = SUMk RESi.CPD(i,k)
c    tre = SUMk RESi.EPD(i,k)

c  Matrix equation becomes:
c  |tcc tec| |delta_Ci| = |trc|
c  |tec tee| |delta_Ei| = |tre|
c
c  Solution is
c    delta_Ci = (trc.tee-tre.tec)/(tee.tcc-tec^2)
c    delta_Ei = (tre.tcc-tre.tec)/(tee.tcc-tec^2)
c
c  if there were no correlation between EPD and CPD, tec=0
c    delta_Ci = trc/tcc
c    delta_Ei = tre/tee
c  this is the "steepest descent" solution in Levenberg-
c  Marquardt speak. By reducing tec or equivalently
c  increasing tee & tcc, the solution can be steered
c  away from least squares minimization and toward
c  steepest descent. Note that this also eliminates
c  the possibility of the denominator (tee.tcc-tec^2)
c  being zero.
c
c  In the event that the PDs are small, i.e. no information
c  from the measurements, then the a priori values dominate:
c  fd is 1.0 in the middle of the spectrum and 0.125 at the edges.
c  fd=0.125+3.5*(i-1)(nmp-i)/nmp^2
c   e = 20000*fd +/- 4000
c   c = 0 +/- 10
c
c   e_apv=20000*fd
c   e_apu=4000*fd
c   c_apv=0
c   c_apu=5
c
c   trc = (c_apv-cbest)/c_apu^2 + Sum[cpd(k)*cres(k)]
c   tcc = 1/c_apu^2 + Sum[cpd(k)^2]
c   tre = (e_apv-ebest)/e_apu^2 + Sum[epd(k)*eres(k)]
c   tee = 1/e_apu^2 + Sum[epd(k)^2]
c   tec = 0 +Sum[epd(k)*cpd(k)]
c  If the Sums are negligible
c    delta_Ci = trc/tcc = c_apv-cbest(i)
c    delta_Ei = tre/tee = e_apv-ebest(i)
c  
c  Compute coefficients tcc,tee,tec,trc,tre
               e_apv=20000.0
               e_apu=1000.0
               c_apv=0.0
               c_apu=1.0
               do i=1,nmp
c                  wte=4*(i-1)*(nmp-i)/eap**2/nmp**2
c                  tee(i)=wte
c                  tcc(i)=apwt
c                  tec(i)=0.0
c                  tre(i)=20000*wte
c                  trc(i)=-cbest(i)*apwt
                  fd=0.125+3.5*(i-1)*float(nmp-i)/nmp**2
                  tee(i)=1.0/(e_apu*fd)**2
                  tcc(i)=1.0/c_apu**2
                  tec(i)=0.0
                  tre(i)=(e_apv*fd-ebest(i))/(e_apu*fd)**2
                  trc(i)=(c_apv-cbest(i))/c_apu**2
                  do kspe=1,nspe
                     tcc(i)=tcc(i)+cpd(i,kspe)**2
                     tee(i)=tee(i)+epd(i,kspe)**2
                     tec(i)=tec(i)+cpd(i,kspe)*epd(i,kspe)
                     tre(i)=tre(i)+epd(i,kspe)*eres(i,kspe)
                     trc(i)=trc(i)+cpd(i,kspe)*cres(i,kspe)
                  end do
               end do  ! i=1,nmp
            endif !  if(rtry.ge.rbest) then

            imax=0
            demax=0.0
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
               if(abs(de).gt.abs(demax)) then
                 imax=i
                 demax=de
               endif
            end do  ! i=1,nmp
c            write(*,*) 'kconv,iter,ebest,cbest=',kconv,iter,
c     &      ebest(10934),cbest(10934)
         end do  ! iter=1,mit
134      continue
         write(*,*)'rbest,rwas,kconv=',rbest,rbestwas,kconv
         if(rbest.ge.rbestwas) then
            kconv=kconv+1
         else
            kconv=0  ! fit improved
            rbestwas=rbest
         endif

c         write(*,*)'Calculating Sk and Dk (scale & shift)'
cc  Iterate on Sk and Dk
cCompute spectrum from c() and e() using resample subroutine
c         do kspe=1,nspe
c            st=stretch(kspe)-stretch(kmin)
c            am=airmass(kspe)
c            call resample(ebest,st+ds/2,kfmin,dpd,nmp) ! Ei"=E(Dk)+ds/2
c            call resample(ebest,st-ds/2,kfmin,se,nmp)  ! Ei"=E(Dk)-ds/2
c            call vsub(dpd,1,se,1,dpd,1,nmp)
c            call vdiv(dpd,1,ds,0,dpd,1,nmp)            ! dE/dDk
c            call resample(ebest,st,kfmin,se,nmp)       ! Ei'=E(Dk)
c            call vmul(cbest,1,am,0,f,1,nmp)            ! Ak.Ci
c            call vexp(f,1,f,1,nmp)                     ! exp(Ak.Ci)
c            call vmul(f,1,dpd,1,dpd,1,nmp)             ! dE/dDk.exp(Ak.Ci)
c            call vmul(se,1,f,1,f,1,nmp)                ! f=Ei'.exp()
c            call vdot(f,1,f,1,y2,nmp)
c            call vsub(y(1,kspe),1,f,1,res,1,nmp)       ! Res = y-f
c            call vdot(res,1,res,1,rms,nmp)
c            call vdot(f,1,f,1,f2,nmp)
c            rmsarr(kspe)=sqrt(rms/f2)
cc
c            tdd=0.0d0
c            trd=0.0d0
c            do i=1,nmp
c               tdd=tdd+dpd(i)**2
c               trd=trd+dpd(i)*res(i)
c            end do  ! i=1,nmp
c            if(kspe.eq.kmin) smin=trd/tdd
cc            stretch(kspe)=stretch(kspe)+trd/tdd
c         end do  ! kspe=1,nspe
c
c         call resample(ebest,+ds/2,kfmin,dpd,nmp) ! Ei" = +ds/2
c         call resample(ebest,-ds/2,kfmin,se,nmp)  ! Ei" = -ds/2
c         call vsub(dpd,1,se,1,dpd,1,nmp)
c         call vdiv(dpd,1,ds,0,dpd,1,nmp)          ! dE/dDk
c         call vdot(dpd,1,dpd,1,dpd2,nmp)
c         call vmul(dpd,1,ebest,1,dpd,1,nmp)
c         call vdot(dpd,1,cbest,1,dpdc,nmp)
cc Reference D and Y to the lowest aimass spectrum.
c         do  kspe=1,nspe
cc            stretch(kspe)=stretch(kspe)-smin
cc     &      +(airmass(kspe)-airmass(kmin))*dpdc/dpd2
c         end do   ! kspe=1,nspe
      end do   ! while(kconv.lt.mconv) 
      close(lunr_rlg)

c  Write out results
      open(lunr_rlg,file=rlgfile,status='old')
      call read_runlog_header(lunr_rlg,data_fmt_read_rl,col_labels_rl)
c      read(lunr_rlg,*) nlhead, ncol
c      do i=2,nlhead
c         read(lunr_rlg,'(a)') col_labels_rl
c      end do
      write(*,*)
      write(72,*)' #     Spectrum       Airmass     Stretch     RMS fit'
      do kspe=1,nspe
10       call read_runlog_data_record(lunr_rlg,data_fmt_read_rl,
     &   col1,specname,iyr,iset,zpdtim,oblat,oblon,obalt,asza,zenoff,
     &   azim,osds,opd,fovi,fovo,amal,ifirst,ilast,
     &   graw,possp_rrl,bytepw_rrl,zoff,snr,apf,tins,pins,hins,
     &   tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
         if(istat.ne.0) exit
         if(col1.eq.':') go to 10
         iabpw=iabs(bytepw_rrl)

         m1=nint(fmin/graw)
         m2=nint(fmax/graw)
c         write(*,'(4i10,2x,a)') ifirst,ilast,m1,m2,specname
         if((m1-ifirst)*float(m2-ilast).gt.0.0) go to 10
         kfmin=m1
         if(m1.lt.ifirst) kfmin=ifirst
         kfmax=m2
         if(m2.gt.ilast) kfmax=ilast
c         possp_wrl=possp_rrl+iabpw*(kfmin-ifirst)
         nmp=kfmax-kfmin+1
         if(nmp.le.0)  go to 10  ! wrong detector spectral region
         write(72,'(i3,1x,a21,3f11.6)')kspe,specname(:21),airmass(kspe),
     &   1000000*(stretch(kspe)-stretch(kmin)),rmsarr(kspe)
      end do
c
c  Write absorption coefficient spectrum (for diagnostic purposes)
      col1=' '
      rlabs=rlmin(:lr-3)//'abs'
      call write_binary_file(19,outpath(:lo)//rlabs,
     & header,0,bytepw_wrl,cbest,nmp)
      call write_runlog_data_record(lunw_rlg,data_fmt_write_rl,
     & col1,rlabs,iyr,iset,zpdtim,oblat,oblon,obalt,asza,zenoff,azim,
     & osds,opd,fovi,fovo,amal,kfmin,kfmax,
     & graw,possp_wrl,bytepw_wrl,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

c  Write exo-atmospheric spectrum (for diagnostic purposes)
      rlext=rlmin(:lr-3)//'ext'
      call write_binary_file(19,outpath(:lo)//rlext,
     & header,possp_wrl,bytepw_wrl,ebest,nmp)
      call write_runlog_data_record(lunw_rlg,data_fmt_write_rl,
     & col1,rlext,iyr,iset,zpdtimmin,oblatmin,oblonmin,obaltmin,0.0,
     & zenoff,azimmin,osdsmin,opdmin,fovi,fovo,amal,kfmin,kfmax,
     & grawmin,0,bytepw_wrl,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

c-------------------------------------------------------------------
c  Perform shaving of smoothed exo-atmospheric spectrum
c  Precompute cosine weighting operator
      do k=1,nkshav
         coswt(k)=cos(k*spi/(2*nkshav+2))**2
      end do
c
c  Apodize.
      resnog=0.5d0/opdmin/grawmin
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

c  Pre-calculate weights.
c  Absorption lines reduce weights from 1.0.
      var=1.5
      do i=1,nmp
         wt(i)=var/(var+cbest(i)**2)
         if(cbest(i).gt.0.0) wt(i)=1.0  ! GCT 2015-03-30
      end do
c
c  Perform weighted SLF to points from i-nkshav to i+nkshav
      do i=1+nkshav,nmp-nkshav
         twt=99*nkshav*wt(i)**2
         tj=0.0
         tjj=0.0
         ty=99*nkshav*ebest(i)*wt(i)**2
         tjy=0.0
         do j=1,nkshav
            twt=twt+(wt(i-j)+wt(i+j))*coswt(j)
            tj=tj+(-wt(i-j)+wt(i+j))*j*coswt(j)
            tjj=tjj+(wt(i-j)+wt(i+j))*j*j*coswt(j)
            ty=ty+(ebest(i-j)*wt(i-j)+ebest(i+j)*wt(i+j))*coswt(j)
            tjy=tjy+((-ebest(i-j))*wt(i-j)+ebest(i+j)*wt(i+j))*
     &      j*coswt(j)
         end do
         denom=twt*tjj-tj*tj
         if(abs(denom).le.0.0) then
            write(*,*) 'denom=0.0'
            denom=0.001
         endif
         etry(i)=(ty*tjj-tjy*tj)/denom
         w=(kfmin+i)*grawmin
         if(abs(w-2330.0).lt.55.0) then
            if(etry(i).lt.12000.) etry(i)=12000.
            if(etry(i).gt.23000.) etry(i)=23000.
         endif
      end do

c  Write shaved EXO file and its runlog entry
      rlexo=rlmin(:lr-3)//'exo'
      call write_binary_file(19,outpath(:lo)//rlexo,
     & header,possp_wrl,bytepw_wrl,etry,nmp)
      call write_runlog_data_record(lunw_rlg,data_fmt_write_rl,
     & col1,rlexo,iyr,iset,zpdtimmin,oblatmin,oblonmin,obaltmin,0.0,
     & zenoff,azimmin,osdsmin,opdmin,fovi,fovo,amal,kfmin,kfmax,
     & grawmin,possp_wrl,bytepw_wrl,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

      close(lunw_rlg)
      close(lunr_rlg)
      close(luno)
      stop
      end
