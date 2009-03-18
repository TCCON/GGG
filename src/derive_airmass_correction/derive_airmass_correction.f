c  Program: derive_airmass_dependence.f
c
c  Purpose: To quantify any airmass-dependent artifacts in the XCO2 data
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
c  SBF is zero at SZA=32 deg, is -0.122 at SZA=0, is 0.875 at SZA=77
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
c  stabilizes the (singular) matrix and sets the ceofficients of
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
c       runlog.vav_ddac.out
c       
      implicit none
      integer*4 lunr,luns,lunw,ncoml,ncol,mcol,kcol,icol,j,kgas,ko2,
     & lnbc,mrow,nmp,imp,ntot,mfp,nfp,
     & iyear,iywas,idoy,idwas,li,
     & kfvsi,kyear,kdoy,klat,klon,ksza, kqcflag, nchar
      parameter (lunr=14,luns=15,lunw=16,mcol=50,mrow=1000,mfp=3)
      character header*800,headarr(mcol)*20,gas*20,
     & inputfile*40,outputfile*40,version*62, specname*35
      real*4 yy(mrow),uy(mrow),bf(mrow,mfp),apx(mfp),apu(mfp),
     & rnorm, uscale,
     & fugas,fuo2,solar_noon,b,eot,diff, qc_threshold
      real*8 yrow(mcol),pi,d2r,chi2

      version=
     &' derive_airmass_correction        1.1.2     2008-03-04     GCT'
      write(*,*) version

      pi=4*datan(1.0d0)
      d2r=pi/180
      nfp=3

      qc_threshold=2.
      nchar=0

      kgas=0
      ko2=0
      kfvsi=0
      kyear=0
      kdoy=0
      klat=0
      klon=0
      ksza=0
      kqcflag=0
      write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav):'
      read(*,'(a)') inputfile
      li=lnbc(inputfile)
      if(inputfile(li-3:li) .ne. '.vav') write(*,*)
     & ' Warning: input file is not of expected type (.vav)'
      open(lunr,file=inputfile, status='old')

      write(*,*)'Enter name of gas:'
      read(*,'(a)') gas

c  Set loose A Priori constraints on ybar, ASDC, & SDC.
c  This prevent the solution going crazy whenever there are
c  fewer than 3 linearly-independent observations in a day.
      if(gas.eq.'co2') then
         apx(1)=380E-06
         apu(1)=200E-06
      elseif(gas.eq.'ch4') then
         apx(1)=1.8E-06
         apu(1)=1.0E-06
      elseif(gas.eq.'n2o') then
         apx(1)=0.3E-06
         apu(1)=0.2E-06
      elseif(gas.eq.'co') then
         apx(1)=0.1E-06
         apu(1)=0.1E-06
      else
         write(*,*)'Unrecognized gas', gas
         stop
      endif
c
c  Set large A Priori constraints on the coefficients of the ASDC and SDC
      do j=2,nfp
        apx(j)=0.0
        apu(j)=10*apu(1)
      end do

c  Read the header of the .vav file and figure out
c  which columns contain the values that we need.
      read(lunr,'(i2,i4)')ncoml,ncol
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml-1
         read(lunr,*)
      end do
      read(lunr,'(a)')header
      call substr(header,headarr,mcol,kcol)
      if (index(header,'Spectrum') .gt. 0) nchar=1
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      do icol=1,ncol
         if(headarr(icol) .eq. gas) kgas=icol
         if(headarr(icol) .eq.  'o2') ko2=icol
         if(headarr(icol) .eq.'fvsi') kfvsi=icol
         if(headarr(icol) .eq.'year') kyear=icol
         if(headarr(icol) .eq. 'day') kdoy=icol
         if(headarr(icol) .eq. 'lat') klat=icol
         if(headarr(icol) .eq.'long') klon=icol
         if(headarr(icol) .eq.'asza') ksza=icol
         if(headarr(icol) .eq.'qcflag') kqcflag=icol
      end do

      outputfile='dac_'//inputfile(:li)//'_'//gas(:lnbc(gas))//'.out'
      open(lunw,file=outputfile,status='unknown')
      write(lunw,*)2,10
      write(lunw,*)'   year     doy  nmp     uscale      ybar    '// 
     &' ybar_error    asdc     asdc_error    sdc    sdc_error'

c  Read each day of data into memory.
      ntot=0
      chi2=0.0d0
      idwas=99999
      iywas=99999
      imp=1
      do while (imp.lt.mrow-nfp)
         yrow(kyear)=0
         if (nchar .eq. 1) then
             read(lunr,*,end=88) specname, (yrow(j),j=1+nchar,ncol)
         else
             read(lunr,*,end=88) (yrow(j),j=1,ncol)
         endif
         if (kqcflag .ne. 0) then
            if( yrow(kqcflag) .lt. qc_threshold) cycle
         endif
         b=2*pi*(yrow(kdoy)-81.0)/364             ! From Wikipedia (364 ???)
         eot=9.87*sin(2*b)-7.53*cos(b)-1.5*sin(b) ! Equation of Time (minutes)
         solar_noon=0.5-yrow(klon)/360-eot/60/24  ! solar noontime (days)
         diff=yrow(kdoy)-solar_noon               ! time wrt solar noon (days)
         idoy=nint(diff)                          ! Convert to Local time
         iyear=int(yrow(kyear))
         if(iyear.ne.iywas .or. idoy.ne.idwas) then ! New day
         if(ntot.gt.0) then
            nmp=imp-1
            call wlsfit(mrow,nmp,yy,nfp,bf,apx,apu,rnorm)
            uscale = rnorm*sqrt(1./nmp)            ! error scale factor 
            chi2=chi2+rnorm**2
            write(lunw,'(f11.5,2i5,7(1pe12.4))') iywas+idwas/365.25,
     &      idwas,nmp,uscale,
     &      (1.E+06*yy(j),1.E+06*sqrt(uscale*bf(j,j)),j=1,nfp)
            imp=1
         endif
         endif

         if(yrow(kfvsi).lt.0.0) yrow(kfvsi)=0.0  ! Don't let missing data affect UY
         fuo2=yrow(ko2+1)/yrow(ko2)          ! fractional uncertainty in column O2
         fugas=yrow(kgas+1)/yrow(kgas)       ! fractional uncertainty in chosen gas
         yy(imp)=0.2095*yrow(kgas)/yrow(ko2) ! Xgas
         uy(imp)=yy(imp)*sqrt(0.00001+fuo2**2+fugas**2
     &  +0.1*yrow(kfvsi)**2                  ! contribution of solar intensity variations
     &  +0.1*(yrow(ksza)/90)**8)             ! de-weight high zenith angles

c   Divide measurements (YY) and basis functions (BF) by uncertainties (UY)
         yy(imp)=yy(imp)/uy(imp)
         bf(imp,1)=1.0/uy(imp)
         bf(imp,2)=sin(2*pi*diff)/uy(imp)
         bf(imp,3)=(((yrow(ksza)+13)/(90+13))**3-((45.0+13)/(90+13))**3)
     &   /uy(imp)
         iywas=iyear
         idwas=idoy
         imp=imp+1
         ntot=ntot+1
      end do              ! while (imp.lt.mrow-nfp)
88    close (lunr)
      nmp=imp-1
      call wlsfit(mrow,nmp,yy,nfp,bf,apx,apu,rnorm)
      uscale = rnorm*sqrt(1./nmp)
      write(lunw,'(f11.5,2i5,7(1pe12.4))') iywas+idwas/365.25,
     & idwas,nmp,uscale,
     & (1.E+06*yy(j),1.E+06*sqrt(uscale*bf(j,j)),j=1,nfp)
      close(lunw)
      nmp=imp-1
      chi2=chi2+rnorm**2
      write(*,'(a5,i6,f9.6)')'ntot=',ntot,sqrt(chi2/ntot)
      stop
      end
