c  PROGRAM: derive_airmass_correction.f
c
c  PURPOSE: To quantify any airmass-dependent artifacts in the XGAS data
c
c  USAGE: 
c       $GGGPATH/bin/derive_airmass_correction xx*.vsw
c   or
c       $GGGPATH/bin/derive_airmass_correction xx*.vav
c       
c  In the first case, airmass corrections for individual windows are
c  reported; in the second, airmass corrections are reported per gas.
c  NOTE: in the first case, there must be a corresponding .vav file
c  because the average H2O column is required to water-correct Xluft.
c  
c  METHODOLOGY:
c
c  Minimizes the function:
c      chi2 = Sum [ (yy(i)-f(i)) / uy(i) ]**2
c  where
c   yy(i) is the XCO2 value measured from the i'th spectrum
c   uy(i) is the associated uncertainty
c   f(i) = cx(1) + cx(2)*ABF(i) + cx(3)*SBF(i)
c   cx(1) is the solar noon XCO2 vmr
c   cx(2) is the coefficient of the anti-symmtric basis function ABF(i)
c   cx(3) is the coefficient of the symmtric basis function SBF(i)
c      ABF(i)=sin(2.Pi.(t(i)-t_solar_noon))
c      SBF(i)=((SZA(i)+13)/(90+13))**3 - ((45+13)/(90+13))**3
c
c  ABF is zero at solar noon, -1 at sunrise, and +1 at sunset.
c  SBF is zero at SZA=45 deg, is -0.18 at SZA=0, is 0.82 at SZA=90
c
c  The minimization is performed with respect to the three unknown
c  coefficients cx. This is done separately for each day.
c
c    Y_hat = CX(1)*1       Daily Mean solar-noon 32 deg XCO2 value
c  ASDC(i) = CX(2)*ABF(i)  Anti-Symmetric Diurnal Component
c   SDC(i) = CX(3)*SBF(i)  Symmetric Diurnal Component
c
c  For days with <3 observation, it is not possible to determine
c  all of the unknowns. In these situations an a priori constraint
c  stabilizes the (singular) matrix and sets the coefficients of
c  ABF and SBF to zero.
c
c  Since we are only really interested in SDC, why is it
c  necessary to also determine YBAR and ASDC for each day?
c  The seasonal variation of XCO2 is much larger than the
c  SDC, at least in the NH. So without first removing the
c  day-to-day variations, it is not possible to get a good
c  estimate of the SDC. Similarly, the variations within a
c  single day are comparable to the SDC, especially in the
c  summer when photo-synthesis draws down the XCO2 each day.
c  To prevent these real variations in XCO2 from confounding
c  the determination of the SDC (especially on days with
c  observations that are asymmetrical about noon), it is
c  essential to fit the real variations.
c
c  The underlying assumption is that the real variations
c  in XCO2 are anti-symmetric about noon. And that any
c  symmetric variations are artifacts. Since the diurnal
c  XCO2 variations driven by photosynthesis & respiration
c  result in a maximum at sunrise and a minimum at sunset,
c  this is probably a good approximation in unpolluted sites.
c  But in polluted locations (e.g. JPL) the pollution signal
c  must first be subtracted by using its correlation with CO.
c
c  Input Files:
c       runlog.vav 
c       airmass_dependence.dat  
c
c  Output Files:
c       dac_runlog.vav_gas.out
c       
      implicit none
      include "../gfit/ggg_int_params.f"
      include "../comn/postproc_params.f"

      integer*4 lunr,lunr_av,lunw,ncoml,ncol,mcol,kcol,icol,j,kgas,ko2,
     & lnbc,mrowpd,nopd,iopd,nfp,nday,doy,i,nlhead,kh2o,
     & iyear,iywas,idoy,idwas,li,naux,nrow,iunder,idum,
     & kfvsi,kyear,kdoy,klat,klon,ksza, kqcflag, mchar,
     & ncoml_av,ncol_av,nrow_av,naux_av,kcol_av, mchar_av
      parameter (lunr=14,lunw=16,lunr_av=18,mcol=200,mrowpd=1000,nfp=3)
      character header*2400,headarr(mcol)*12,gas*4,
     & header_av*800,headarr_av(mcol)*12,
     & inputfile*40,outputfile*40,version*62, specname*(nchar),
     & avfile*40, specname_av*(nchar), gas_window*10,gsza_cl*10,
     & p_cl*10,gas_cl*10
      real*4 
     & yy(mrowpd),uy(mrowpd),
     & bf(mrowpd,nfp),apx(nfp),apu(nfp),
     & rnorm, uscale,
     & qc_threshold
      real*8 
     & aa(nfp),ae(nfp), year, 
     & ybar,
     & tyy(nfp),ty(nfp),tee(nfp),
     & fugas,fuo2,
     & solar_noon,b,eot,diff,
     & yrow(mcol),chi2,dpi,
     & yrow_av(mcol),gsza,p
      logical isavg,onegas
      parameter(dpi=3.14159265359D0)


c  Do-nothing code to use parameters in include files
      idum=mauxcol    ! Avoid compiler warning (unused parameter)
      idum=mcolvsw    ! Avoid compiler warning (unused parameter)
      idum=mcolvav    ! Avoid compiler warning (unused parameter)
      idum=mfilepath  ! Avoid compiler warning (unused parameter)
      idum=mgas       ! Avoid compiler warning (unused parameter)
      idum=mlev       ! Avoid compiler warning (unused parameter)
      idum=mrow_qc    ! Avoid compiler warning (unused parameter)
      idum=mspeci     ! Avoid compiler warning (unused parameter)
      idum=mvmode     ! Avoid compiler warning (unused parameter)
      idum=ncell      ! Avoid compiler warning (unused parameter)
      idum=maddln     ! Avoid compiler warning (unused parameter)
      idum=mcharhead  ! Avoid compiler warning (unused parameter)
    

      version=
     &' derive_airmass_correction        1.27    2020-12-16   GCT/JLL'
      write(*,*) version

      qc_threshold=2.
      mchar=0
      mchar_av=0

      kgas=0
      ko2=0
      kh2o=0
      kfvsi=0
      kyear=0
      kdoy=0
      klat=0
      klon=0
      ksza=0
      kqcflag=0

      gsza = 13.0
      p = 3.0
      gas_cl = '          '
      if (iargc() == 0) then
         write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav):'
         read(*,'(a)') inputfile
      elseif (iargc() .ge. 1) then
         call getarg(1, inputfile)
      else
         stop 'Usage: $gggpath/bin/derive_airmass_correction vavfile'
      endif

      if (iargc() .ge. 2) then
         call getarg(2, gsza_cl)
         read(gsza_cl, '(f10.0)') gsza
      endif
      if (iargc() .ge. 3) then
         call getarg(3, p_cl)
         read(p_cl, '(f10.0)') p
      endif
      if (iargc() .ge. 4) then
         call getarg(4, gas_cl)
      endif
c      write(*,'(a,f10.2)') 'Will use gsza=', gsza
      li=lnbc(inputfile)
      onegas = lnbc(gas_cl) .gt. 0

c Determine is working from an averaged file or not
      if(inputfile(li-3:li) .eq. '.vav') then
        isavg = .true.
      elseif(inputfile(li-3:li) .eq. '.vsw') then
        isavg = .false.
      elseif(inputfile(li-1:li) .eq. 'av') then
        write(*,*) 'Warning: input file is one of the expected '//
     &  'types (.vav or .vsw), assuming averaged quantities'
        isavg = .true.
      else
        write(*,*) 'Warning: input file is one of the expected '//
     &  'types (.vav or .vsw), assuming unaveraged quantities'
        isavg = .false.
      end if

      if(.not. isavg) avfile=inputfile(:li-2)//'av'


      open(lunr,file=inputfile, status='old')
      read(lunr,countfmt)ncoml,ncol,nrow,naux
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml-1
         read(lunr,*)
      end do
      read(lunr,'(a)')header
      call substr(header,headarr,mcol,kcol)
      call lowercase(header)
      if (index(header,'spectrum') .gt. 0) mchar=1
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      close(lunr)

c  If deriving airmass corrections for specific windows then
c  we still need the average H2O to water-correct Luft, so 
c  we're going to need to read in the average file as well
      if (.not. isavg) then
        open(lunr_av, file=avfile, status='old')
        read(lunr_av,countfmt)ncoml_av,ncol_av,nrow_av,naux_av
        do j=2,ncoml_av-1
          read(lunr_av,*)
        end do
        read(lunr_av,'(a)')header_av
        call substr(header_av, headarr_av, mcol, kcol_av)
        call lowercase(header_av)
        write(*,*) 'header_av=',header_av(1:32)
        if (index(header_av, 'spectrum') .gt. 0) mchar_av=1
        if (kcol_av .ne. ncol_av) stop 'ncol_av/kcol_av mismatch'
        if (nrow_av .ne. nrow) stop 
     & '.vsw and .vav have different number of rows'
        close(lunr_av)
      end if
      
c      write(*,*)'Enter name of gas:'
c      read(*,'(a)') gas

c     write(*,*)'naux,ncol=',naux,ncol
      write(*,'(2(a,f10.6))') 'gsza = ', gsza, ', p = ', p
      write(*,'(a16,6a12)') 'gas','ybar','ybar_error','asdc',
     &'asdc_error','sdc','sdc_error'
      do kgas=naux+1,ncol,2

c  JLL: get which gas or window we are fitting in this loop
         if (isavg) then
            gas=headarr(kgas)(1:4)
            gas_window=gas
         else
            iunder=lnbc(headarr(kgas))
            gas_window=headarr(kgas)(1:iunder)
            iunder=index(headarr(kgas),'_')-1
            gas=headarr(kgas)(1:iunder)
         endif ! isavg

c  JLL: if we're only supposed to fit one gas, and this loop
c  isn't that gas, skip to the next one.
         if (onegas .and. gas_window(1:lnbc(gas_window)) .ne.
     &       gas_cl(1:lnbc(gas_cl))) cycle

c  Read the header of the primary input file 
c  (.vav or .vsw) and figure out which columns 
c  contain the values that we need.
         open(lunr,file=inputfile, status='old')
         read(lunr,countfmt)ncoml,ncol,nrow,naux
         if(ncol.gt.mcol) stop 'increase mcol'
         do j=2,ncoml
            read(lunr,*)
         end do
         do icol=1,ncol
            if(isavg) then
              if(headarr(icol) .eq.  'o2') ko2=icol
              if(headarr(icol) .eq. 'h2o') kh2o=icol
            else
              if(headarr(icol) .eq. 'o2_7885') ko2=icol
            end if
            if(headarr(icol) .eq.'fvsi') kfvsi=icol
            if(headarr(icol) .eq.'year') kyear=icol
            if(headarr(icol) .eq. 'day') kdoy=icol
            if(headarr(icol) .eq. 'lat') klat=icol
            if(headarr(icol) .eq.'long') klon=icol
            if(headarr(icol) .eq.'asza') ksza=icol
            if(headarr(icol) .eq.'solzen') ksza=icol
            if(headarr(icol) .eq.'qcflag') kqcflag=icol
         end do

c  If deriving window specific airmass corrections, then we
c  also need to find the average water column in the .vav
c  file to do the luft water correction
         if (.not. isavg) then
           open(lunr_av, file=avfile, status='old')
           read(lunr_av,countfmt)ncoml_av,ncol_av,nrow_av,naux_av
           if(ncol_av .gt. mcol) stop 'increase mcol (from vav)'
           do j=2,ncoml_av
             read(lunr_av,*)
           end do

           do icol=1,ncol_av
             if(headarr_av(icol) .eq. 'h2o') kh2o=icol
           end do ! icol=1,ncol_av
         endif ! (.not isavg)


c  Set loose A Priori constraints on ybar, ASDC, & SDC.
c  This prevents the solution going crazy whenever there are
c  fewer than 3 linearly-independent observations in a day.
         if(gas.eq.'co2'.or.gas.eq.'lco2'.or.gas.eq.'wco2') then
            apx(1)=400E-06
            apu(1)=300E-06
         elseif(gas.eq.'ch4') then
            apx(1)=1.8E-06
            apu(1)=1.4E-06
         elseif(gas.eq.'n2o') then
            apx(1)=0.3E-06
            apu(1)=0.3E-06
         elseif(gas.eq.'co') then
            apx(1)=0.1E-06
            apu(1)=0.2E-06
         elseif(gas.eq.'ocs') then
            apx(1)=0.1E-09
            apu(1)=0.2E-09
         elseif(gas.eq.'luft') then
            apx(1)=1.0E-00
            apu(1)=2.0E-00
         elseif(gas.eq.'hcl') then
            apx(1)=1E-10
            apu(1)=5E-11
         else
c           write(*,*)'Skipping column/gas: ',kgas, gas
            close(lunr)
            if (.not. isavg) close(lunr_av)
            cycle
         endif
c
c  Set large A Priori constraints on the coefficients of the ASDC and SDC
         do j=2,nfp
            apx(j)=0.0
            apu(j)=10*apu(1)
         end do

         outputfile='dac_'//inputfile(:li)//'_'//
     &      gas_window(:lnbc(gas_window))//'.out'
         open(lunw,file=outputfile,status='unknown')
         write(lunw,'(2i5)')2,11
         write(lunw,'(a)')'gas    year     doy  nopd    uscale      ybar
     &   ybar_error     asdc     asdc_error     sdc      sdc_error'

c  Read each day of data into memory.
         nday=0
         chi2=0.0d0
         idwas=99999
         iywas=99999
         iopd=1
         do while (iopd.lt.mrowpd-nfp)  ! Observations Per Day
            yrow(kyear)=0
            if (mchar .eq. 1) then
               read(lunr,*,end=88) specname, (yrow(j),j=1+mchar,ncol)
            else
               read(lunr,*,end=88) (yrow(j),j=1,ncol)
            endif

            if (mchar_av .eq. 1 .and. .not. isavg) then
              read(lunr_av,*) specname_av,
     &          (yrow_av(j),j=1+mchar_av,ncol_av)
              if(specname .ne. specname_av) stop 
     & '.vsw and .vav spectrum name mismatch'
            elseif(.not. isavg) then
              read(lunr_av,*) (yrow_av(j),j=1,ncol_av)
            endif

            if (kqcflag .ne. 0) then
               if( yrow(kqcflag) .lt. qc_threshold) cycle
            endif
            b=2*dpi*(yrow(kdoy)-81.0)/364             ! From Wikipedia (364 ???)
            eot=9.87*sin(2*b)-7.53*cos(b)-1.5*sin(b) ! Equation of Time (minutes)
            solar_noon=0.5-yrow(klon)/360-eot/60/24  ! solar noontime (days)
            diff=yrow(kdoy)-solar_noon               ! time wrt solar noon (days)
            idoy=nint(diff)                          ! Convert to Local time
            iyear=int(yrow(kyear))
            if(iyear.ne.iywas .or. idoy.ne.idwas) then ! New day
               if(iopd.gt.1) then
                  nopd=iopd-1
                  if(nopd.gt.nfp) then
c                     JLL: only compute and write the basis function
c                     coefficients if there are at least as many
c                     spectra as basis functions. While the prior
c                     constraints keep the fitting sane, it can
c                     result in a fit that has rnorm = 0 (i.e. 
c                     there is no residual between Ax and b)
c                     which makes uscale 0 and thus the errors 0,
c                     which results in a NaN for the airmass CF.

                      call wlsfit(mrowpd,nopd,yy,nfp,bf,apx,apu,rnorm)
                      uscale = rnorm*sqrt(1./nopd)            ! error scale factor 
                      chi2=chi2+rnorm**2
                      write(lunw,'(a,f11.5,2i5,7(1pe12.4))')gas_window,
     &                iywas+float(idwas)/365.25,
     &                idwas,nopd,uscale,
     &                (1.E+06*yy(j),1.E+06*sqrt(uscale*bf(j,j)),j=1,nfp)
                      iopd=1
                      nday=nday+1
                  endif ! (nopd.gt.nfp)
               endif
            endif

            if(yrow(kfvsi).lt.0.0) yrow(kfvsi)=0.0  ! Don't let missing data affect UY
            fuo2=yrow(ko2+1)/yrow(ko2)          ! fractional uncertainty in column O2
            fugas=yrow(kgas+1)/yrow(kgas)       ! fractional uncertainty in chosen gas
            yy(iopd)=0.2095*sngl(yrow(kgas)/yrow(ko2)) ! Xgas
c Xluft needs to have a water correction before its airmass dependence is calculated
            if(gas.eq.'luft') then 
               if(kh2o.gt.0 .and. isavg) then
                  yy(iopd)=sngl(0.2095*yrow(kgas)/yrow(ko2)-
     &                    0.2095*yrow(kh2o)/yrow(ko2)*18.02/28.964)
               elseif(kh2o.gt.0) then
                  yy(iopd)=sngl(0.2095*yrow(kgas)/yrow(ko2)-
     &                    0.2095*yrow_av(kh2o)/yrow(ko2)*18.02/28.964)
               else
                  write(*,*)'Warning: no H2O correction on Xluft!'
               endif
            endif
            uy(iopd)=yy(iopd)*sngl(dsqrt(0.00001d0+fuo2**2+fugas**2
     &     +0.1*yrow(kfvsi)**2))                ! contribution of solar intensity variations
c   JLL: removed de-weighting of high SZAs for GGG2020; testing showed
c   that allowing these points full weight did a better job capturing
c   the airmass dependence of individual windows.
c
c     &     +0.1*(yrow(ksza)/90)**8)             ! de-weight high zenith angles
c     &     +0.2*(yrow(ksza)/90)**12))             ! de-weight high zenith angles

c   Divide measurements (YY) and basis functions (BF) by uncertainties (UY)
            yy(iopd)=yy(iopd)/uy(iopd)
            bf(iopd,1)=sngl(1.0d0/uy(iopd))
            bf(iopd,2)=sngl(dsin(2*dpi*diff)/uy(iopd))
            bf(iopd,3)=sngl((((yrow(ksza)+gsza)/(90+gsza))**p
     &      -((45.+gsza)/(90+gsza))**p)/uy(iopd))
            iywas=iyear
            idwas=idoy
            iopd=iopd+1
         end do              ! while (iopd.lt.mrowpd-nfp)
88       close (lunr)
         if (.not. isavg) close (lunr_av)
         nopd=iopd-1
c  
c  Do the last day.
c  JLL: unless the last day had fewer spectra than the number
c       of basis functions.
         if (nopd .ge. nfp) then
             call wlsfit(mrowpd,nopd,yy,nfp,bf,apx,apu,rnorm)
             uscale = rnorm*sqrt(1./nopd)
             write(lunw,'(a,f11.5,2i5,7(1pe12.4))') gas_window,
     &        iywas+idwas/365.25,idwas,nopd,uscale,
     &        (1.E+06*yy(j),1.E+06*sqrt(uscale*bf(j,j)),j=1,nfp)
             chi2=chi2+rnorm**2
         endif ! (nopd .ge. nfp)

         close(lunw)
         nopd=iopd-1

         open(lunw,file=outputfile,status='old')
         read(lunw,*) nlhead,ncol
         do i=2,nlhead
            read(lunw,*)
         end do
         do j=1,nfp
            tee(j)=0.0d0
            tyy(j)=0.0d0
            ty(j)=0.0d0
         end do
         do i=1,nday
            read(lunw,*) gas_window,year,doy,nopd,uscale,
     & (aa(j),ae(j),j=1,nfp)
            do j=1,nfp
c               taa(j)=taa(j)+aa(j)/ae(j)**2
c               tae(j)=tae(j)+1/ae(j)**2
               ty(j)=ty(j)+aa(j)/ae(j)**2
               tyy(j)=tyy(j)+(aa(j)/ae(j))**2
               tee(j)=tee(j)+1/ae(j)**2
            end do
         end do
         close(lunw)
         write(*,*)
c        write(*,*) 'gas     ybar   ybar_error      asdc    asdc_error   
c    &    sdc     sdc_error'
         write(*,'(a,6(1x,f11.4))')'X'//gas_window//'(ppm)',
     &    (ty(j)/tee(j),
     &    sqrt(tyy(j)/tee(j)-(ty(j)/tee(j))**2),j=1,nfp)
         ybar=ty(1)/tee(1)
         write(*,'(a,6(1x,f11.5))') 'X'//gas_window//'/ybar',
     &    (ty(j)/tee(j)/ybar,
     &    sqrt(tyy(j)/tee(j)-(ty(j)/tee(j))**2)/ybar,j=1,nfp)

      end do        !  kgas=naux+1,ncol,2
      stop
      end
