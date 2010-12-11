c  Program derive_insitu_correction.f
c  
c  Input Files:
c       runlog.vav.ada 
c       aircraft_overpasses.dat  
c
c  Output Files:
c       Writes the O2 and CO2 calibration factors to the screen.
c       
c  Program compares the O2 columns with the calculated air columns
c  to determine an average XO2 value for the entire dataset.
c
c  Program also compares the FTS XCO2 with values derived from
c  various aircraft over-flights.  This is evaluated only for
c  spectra that overlap the duration of the aircraft profile.
c
      implicit none
      integer*4 lunr,luns,lunw,ncoml,ncol,mcol,kcol,icol,i,j,kco2,ko2,
     & iy1,id1,ih1,im1,is1,iy2,id2,ih2,im2,is2,nap,lnbc,lg,kap,
     & nchar,
     & kqcflag,
     & nrow,irow,kair,kh2o,kfvsi,kyear,klat,klon,ksza
      parameter (lunr=14,luns=15,lunw=16,mcol=50)
      character header*800, headarr(mcol)*20, gggdir*80, inputfile*40,
     & specname*35, version*62
      real*8 yrow(mcol),alat,alon,zz,tk,ratio,
     & tt,ty,ty2,tsza, 
     & qc_threshold,
     & s0,s1,s2,t0,t1,t2,tn,tt0,tt1,tt2,ttn,tstart,tstop,
     & xco2_ac,xco2_ac_err,yy,wt,fco2,fxco2_ac,fvsi,fair,fh2o

      version=
     &' derive_insitu_correction         1.1.3     2010-11-23     GCT'
      write(*,*) version

      kair=0
      kh2o=0
      kco2=0
      ko2=0
      kfvsi=0
      kyear=0
      klat=0
      klon=0
      ksza=0
      kqcflag=0

      nchar=0
      qc_threshold=2.

      call getenv("GGGPATH", gggdir)
      write(*,*)gggdir
      lg=lnbc(gggdir)
      write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav.ada):'
      read(*,'(a)') inputfile
      open(lunr,file=inputfile, status='old')
      read(lunr,*)ncoml,ncol
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml-1
         read(lunr,*)
      end do
      read(lunr,'(a)')header
      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      if (index(header,'Spectrum') .gt. 0) nchar=1
      open(lunw,file='dic_'//
     & inputfile(:lnbc(inputfile))//'.out',status='unknown')
      do icol=1,ncol
c         write(*,*)icol,headarr(icol)
c  Find out which columns of the .vav file contain the required parameters
c  This allows flexibility in the order of the columns.
         if(headarr(icol) .eq. 'xair') kair=icol
         if(headarr(icol) .eq. 'xh2o') kh2o=icol
         if(headarr(icol) .eq. 'xco2') kco2=icol
         if(headarr(icol) .eq.  'xo2') ko2=icol
         if(headarr(icol) .eq. 'fvsi') kfvsi=icol
         if(headarr(icol) .eq.'year') kyear=icol
         if(headarr(icol) .eq. 'lat') klat=icol
         if(headarr(icol) .eq.'long') klon=icol
         if(headarr(icol) .eq.'asza') ksza=icol
         if(headarr(icol) .eq.'qcflag') kqcflag=icol
      end do

      if(kco2.eq.0) stop 'xco2 missing from input file  '
      if(ko2.eq.0)  stop ' xo2 missing from input file  '
      if(kh2o.eq.0) stop 'xh2o missing from input file  '
      if(kair.eq.0) stop 'xair missing from input file  '
c      write(*,*)kair,kh2o,kco2,ko2,kfvsi,kyear,klat,klon,ksza

      tk=0.0d0
      tt=0.0d0
      ty=0.0d0
      ty2=0.0d0
      do irow=1,999999
         if (nchar .eq. 1) then
             read(lunr,*,end=88) specname, (yrow(j),j=1+nchar,ncol)
         else
             read(lunr,*,end=88) (yrow(j),j=1,ncol)
         endif
         if (kqcflag .ne. 0) then
            if( yrow(kqcflag) .lt. qc_threshold) cycle
         endif
         fvsi=yrow(kfvsi)                   ! fractional uncertainty
         if(fvsi.lt.0.) fvsi=0.             ! Don't let missing data affect WT
         fair=yrow(kair+1)/yrow(kair)       ! fractional uncertainty
         fh2o=yrow(kh2o+1)/yrow(kh2o)       ! fractional uncertainty
         wt=1/(0.00001+fair**2+(fh2o*18.02/28.964)**2+fvsi**2)
         zz=yrow(ko2)/(yrow(kair)-yrow(kh2o)*18.02/28.964)
         tk=tk+1
         tt=tt+wt
         ty=ty+wt*zz
         ty2=ty2+wt*zz**2
      end do   ! irow=1,999999
      close (lunr)
      nrow=irow-1
      write(*,*)
      write(lunw,*)
      write(*,*)' Total Number of FTS observations=',nrow
      write(lunw,*)'Total Number of FTS observations=',nrow
      write(*,*)' Mean O2/(Air-H2O)/0.2095 =', ty/tt/0.2095,
     & ' +/-', sqrt(ty2*tt-ty**2)/tt/0.2095
      write(lunw,*)' Mean O2/(Air-H2O)/0.2095 =', ty/tt/0.2095,
     & ' +/-', sqrt(ty2*tt-ty**2)/tt/0.2095
      write(*,*)
      write(*,*)
      write(lunw,'(a)')
     & 'Year  DOY  NOBS   SZA    XCO2_AC Error  XCO2_FTS Error
     &  Ratio   Error'
      write(*,*)'Year  DOY  NOBS   SZA    XCO2_AC Error  XCO2_FTS Error
     &  Ratio   Error'

      ttn=0.0
      tt0=0.0
      tt1=0.0
      tt2=0.0
c  Read in line of aircraft profile XCO2
      open(luns,file=gggdir(:lg)//'/tccon/aircraft_overflights.dat',
     &  status='old')
      read(luns,*)
      do kap=1,999
         tn=0.0
         t0=0.0
         t1=0.0
         t2=0.0
         tsza=0.0
         s0=0.0
         s1=0.0
         s2=0.0
         read(luns,*,end=88)alat,alon,iy1,id1,ih1,im1,is1,
     &   iy2,id2,ih2,im2,is2,xco2_ac,xco2_ac_err
         fxco2_ac=xco2_ac_err/xco2_ac
         tstart=iy1+(id1+(ih1+(im1+float(is1)/60)/60)/24)/365.25
         tstop =iy2+(id2+(ih2+(im2+float(is2)/60)/60)/24)/365.25

c Read the FTS data file, again
         open(lunr,file=inputfile, status='old')
         do j=1,ncoml
           read(lunr,*) 
         end do
         do i=1,nrow
            if (nchar .eq. 1) then
               read(lunr,*,end=88) specname, (yrow(j),j=1+nchar,ncol)
            else
               read(lunr,*,end=88) (yrow(j),j=1,ncol)
            endif
            if (kqcflag .ne. 0) then
               if( yrow(kqcflag) .lt. qc_threshold) cycle
            endif
c  If a measurement coincides with the current aircraft overpass      
            if(yrow(kyear).ge.tstop) exit
            if(yrow(kyear).ge.tstart) then
               if(abs(yrow(klat)-alat) .lt. 0.2) then
                  if(abs(yrow(klon)-alon) .lt. 0.5) then
                     fco2=yrow(kco2+1)/yrow(kco2)
                     fvsi=yrow(kfvsi)   ! fractional uncertainty
                     if(fvsi.lt.0.0) fvsi=0.0  ! Don't let missing data affect WT
                     wt=1/(0.00001+fco2**2+fxco2_ac**2+fvsi**2)
                     yy=yrow(kco2)/xco2_ac
                     tn=tn+1
                     t0=t0+wt
                     t1=t1+wt*yy
                     t2=t2+wt*yy**2
                     tsza=tsza+wt*yrow(ksza)
                     zz=yrow(kco2)
                     wt=1/(0.00001+fco2**2+fvsi**2)
                     s0=s0+wt
                     s1=s1+wt*zz
                     s2=s2+wt*zz**2
                  endif
               endif
            endif
         end do   ! irow=1,nrow
         close(lunr)
         ratio=t1/t0  !  Ratio: FTS/Aircraft XCO2
         if(nint(tn).gt.0) write(*,'(3i5,f7.1,2(f10.2,f6.2),2f9.4)')
     &  iy1,id1,nint(tn), tsza/t0,
     &  1.0e+06*xco2_ac, 1.0e+06*xco2_ac_err,
     &  1.0e+06*s1/s0,1.0e+06*sqrt(s0*s2-s1**2)/s0,
     &  ratio,ratio*sqrt(t0*t2-t1**2)/t1
         if(nint(tn).gt.0) write(lunw,'(3i5,f7.1,2(f10.2,f6.2),2f9.4)')
     &  iy1,id1,nint(tn), tsza/t0,
     &  1.0e+06*xco2_ac, 1.0e+06*xco2_ac_err,
     &  1.0e+06*s1/s0,1.0e+06*sqrt(s0*s2-s1**2)/s0,
     &  ratio,ratio*sqrt(t0*t2-t1**2)/t1
        ttn=ttn+tn
        tt0=tt0+t0
        tt1=tt1+t1
        tt2=tt2+t2
      end do   ! kap=1,999
88    close(luns)
      nap=kap-1  ! Number of aircraft profiles
      ratio=tt1/tt0  !  Ratio: FTS/Aircraft XCO2
      write(*,*) 'Average over all aircraft profiles ',
     & ratio,' +/-',ratio*sqrt(tt0*tt2-tt1**2)/tt1
      write(lunw,*) 'Average over all aircraft profiles ',
     & ratio,' +/-',ratio*sqrt(tt0*tt2-tt1**2)/tt1
      close(lunw)
      stop
      end
