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
      include "../gfit/ggg_int_params.f"

      integer*4 lunr,luns,lunw,ncoml,ncol,mcol,kcol,icol,i,j,kco2,ko2,
     & iy1,id1,ih1,im1,is1,iy2,id2,ih2,im2,is2,nap,lnbc,lg,lgs,kap,
     & mchar,idum,
     & kqcflag,
     & nrow,irow,kluft,kh2o,kfvsi,kyear,klat,klon,ksza
      parameter (lunr=14,luns=15,lunw=16,mcol=50)
      character header*800, headarr(mcol)*20, gggdir*(mpath),
     &  inputfile*40,
     & specname*(nchar), version*62, gas*4,dl*1
      real*8 yrow(mcol),alat,alon,zz,tk,ratio,
     & tt,ty,ty2,tsza, 
     & qc_threshold,
     & s0,s1,s2,t0,t1,t2,tn,tt0,tt1,tt2,ttn,tstart,tstop,
     & xco2_ac,xco2_ac_err,yy,wt,fco2,fxco2_ac,fvsi,fluft,fh2o


      idum=mauxcol   ! Avoid compiler warning (unused parameter)
      idum=mcolvav   ! Avoid compiler warning (unused parameter)
      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mgas      ! Avoid compiler warning (unused parameter)
      idum=mlev      ! Avoid compiler warning (unused parameter)
      idum=mrow_qc   ! Avoid compiler warning (unused parameter)
      idum=mspeci    ! Avoid compiler warning (unused parameter)
      idum=mvmode    ! Avoid compiler warning (unused parameter)
      idum=ncell     ! Avoid compiler warning (unused parameter)

      version=
     &' derive_insitu_correction         1.15      2017-02-16     GCT'
      write(*,*) version

      kluft=0
      kh2o=0
      kco2=0
      ko2=0
      kfvsi=0
      kyear=0
      klat=0
      klon=0
      ksza=0
      kqcflag=0

      mchar=0
      qc_threshold=2.
      gas = '    '

      call get_ggg_environment(gggdir, dl)
      write(*,*)gggdir
      lg=lnbc(gggdir)
      if (iargc() == 0) then
         write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav.ada):'
         read(*,'(a)') inputfile
      elseif (iargc() == 2) then
         call getarg(1, inputfile)
         call getarg(2, gas(1:1))
      else
         write(*,*)'Use: $gggpath/bin/derive_insitu_correction '//
     &   'adafile gas# (1=xco2,2=xco,3=xch4,4=xn2o)'
         stop
      endif
      open(lunr,file=inputfile, status='old')
      read(lunr,*)ncoml,ncol
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml-1
         read(lunr,*)
      end do
      read(lunr,'(a)')header

c Determine which gas to process
      if (iargc() == 0) then
         write(6,9913)
9913     format(' Gas (1=xco2,2=xco,3=xch4,4=xn2o) ? ',$)
         read(5,'(a)') gas(1:1)
      endif
      if(gas(1:1).eq.'1') gas(1:3)='co2'
      if(gas(1:1).eq.'2') gas(1:3)='co '
      if(gas(1:1).eq.'3') gas(1:3)='ch4'
      if(gas(1:1).eq.'4') gas(1:3)='n2o'
      if(gas(2:3).eq.'  ') stop 'Unknown gas'
      write(*,*)'Processing ',gas
      lgs=lnbc(gas)

      call substr(header,headarr,mcol,kcol)
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'
      call lowercase(header)
      if (index(header,'spectrum').gt.0) mchar=1
      open(lunw,file='dic_'//
     & inputfile(:lnbc(inputfile))//'_'//gas(:lgs)//'.out',
     & status='unknown')
      do icol=1,ncol
c         write(*,*)icol,headarr(icol)
c  Find out which columns of the .vav file contain the required parameters
c  This allows flexibility in the order of the columns.
         if(headarr(icol) .eq. 'xluft') kluft=icol
         if(headarr(icol) .eq. 'xh2o') kh2o=icol
         if(headarr(icol) .eq. 'x'//gas(:lgs)) kco2=icol
         if(headarr(icol) .eq.  'xo2') ko2=icol
         if(headarr(icol) .eq. 'fvsi') kfvsi=icol
         if(headarr(icol) .eq.'year') kyear=icol
         if(headarr(icol) .eq. 'lat') klat=icol
         if(headarr(icol) .eq.'long') klon=icol
         if(headarr(icol) .eq.'asza') ksza=icol
         if(headarr(icol) .eq.'qcflag') kqcflag=icol
      end do

      if(kco2.eq.0) stop ' xgas missing from input file  '
      if(ko2.eq.0)  stop '  xo2 missing from input file  '
      if(kh2o.eq.0) stop ' xh2o missing from input file  '
      if(kluft.eq.0)stop 'xluft missing from input file  '
c      write(*,*)kluft,kh2o,kco2,ko2,kfvsi,kyear,klat,klon,ksza

      tk=0.0d0
      tt=0.0d0
      ty=0.0d0
      ty2=0.0d0
c     write(*,*)'mchar=',mchar
      do irow=1,999999
         if (mchar .eq. 1) then
            read(lunr,*,end=77) specname, (yrow(j),j=1+mchar,ncol)
         else
            read(lunr,*,end=77) (yrow(j),j=1,ncol)
         endif
         if (kqcflag .ne. 0) then
            if( yrow(kqcflag) .lt. qc_threshold) cycle
         endif
         fvsi=yrow(kfvsi)                   ! fractional uncertainty
         if(fvsi.lt.0.) fvsi=0.             ! Don't let missing data affect WT
         fluft=yrow(kluft+1)/yrow(kluft)    ! fractional uncertainty
         fh2o=yrow(kh2o+1)/yrow(kh2o)       ! fractional uncertainty
         wt=1/(0.00001+fluft**2+(fh2o*18.02/28.964)**2+fvsi**2)
         zz=yrow(ko2)/(yrow(kluft)-yrow(kh2o)*18.02/28.964)
         tk=tk+1
         tt=tt+wt
         ty=ty+wt*zz
         ty2=ty2+wt*zz**2
c        write(*,*)'fluft,fvsi,fh2o,fyear=',fluft,fvsi,fh2o,yrow(kyear)
c         write(*,*)'yrow(kyear)',yrow(kyear)
      end do   ! irow=1,999999
77    close (lunr)
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
     & 'Year  DOY  NOBS   SZA    X'//gas(:lgs)
     &  //'_AC Error  X'//gas(:lgs)//'_FTS Error  Ratio    Error'
      write(*,*)'Year  DOY  NOBS   SZA    X'//gas(:lgs)//
     & '_AC Error  X'//gas(:lgs)//'_FTS Error  Ratio    Error'

      ttn=0.0
      tt0=0.0
      tt1=0.0
      tt2=0.0
c  Read in line of aircraft profile XCO2

      open(luns,file=gggdir(:lg)//'tccon'//dl//'aircraft_overflights_'
     & //gas(:lgs)//'.dat',status='old')
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
c        tstart=iy1+(id1+(ih1+(im1+float(is1)/60)/60)/24)/365.25
c        tstop =iy2+(id2+(ih2+(im2+float(is2)/60)/60)/24)/365.25
         tstart=iy1+(id1+(ih1+(im1+float(is1)/60)/60)/24)/366 ! this is to match collate_results.f
         tstop =iy2+(id2+(ih2+(im2+float(is2)/60)/60)/24)/366
c        write(*,*)'tstart,tstop,yrow(kyear)',tstart,tstop,yrow(kyear)

c Read the FTS data file, again
         open(lunr,file=inputfile, status='old')
         do j=1,ncoml
            read(lunr,*) 
         end do
         do i=1,nrow
            if (mchar .eq. 1) then
               read(lunr,*,end=66) specname, (yrow(j),j=1+mchar,ncol)
            else
               read(lunr,*,end=66) (yrow(j),j=1,ncol)
            endif
            if (kqcflag .ne. 0) then
               if( yrow(kqcflag) .lt. qc_threshold) cycle
            endif
c  If a measurement coincides with the current aircraft overpass      
c           write(*,*)'tstart,tstop,fyear',tstart,tstop,yrow(kyear)
            if(yrow(kyear).ge.tstop) exit
            if(yrow(kyear).ge.tstart) then
c              write(*,*)'yrow(klat),alat,yrow(klon),alon,yrow(kyear)=',
c    & yrow(klat),alat,yrow(klon),alon,yrow(kyear)
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
66       close(lunr)
         ratio=t1/t0  !  Ratio: FTS/Aircraft XCO2
c        write(*,*)'tn=',tn
         if(nint(tn).gt.0) write(*,'(3i5,f7.1,2(f10.2,f6.2),2f9.4)')
     &   iy1,id1,nint(tn), tsza/t0,
     &   1.0e+06*xco2_ac, 1.0e+06*xco2_ac_err,
     &   1.0e+06*s1/s0,1.0e+06*sqrt(s0*s2-s1**2)/s0,
     &   ratio,ratio*sqrt(t0*t2-t1**2)/t1
         if(nint(tn).gt.0) write(lunw,'(3i5,f7.1,2(f10.2,f6.2),2f9.4)')
     &   iy1,id1,nint(tn), tsza/t0,
     &   1.0e+06*xco2_ac, 1.0e+06*xco2_ac_err,
     &   1.0e+06*s1/s0,1.0e+06*sqrt(s0*s2-s1**2)/s0,
     &   ratio,ratio*sqrt(t0*t2-t1**2)/t1
         ttn=ttn+tn
         tt0=tt0+t0
         tt1=tt1+t1
         tt2=tt2+t2
      end do   ! kap=1,999
88    close(luns)
      nap=kap-1  ! Number of aircraft profiles
c     write(*,*)'nap=',nap
      ratio=tt1/tt0  !  Ratio: FTS/Aircraft XCO2
c     write(*,*)'tt1,tt0=',tt1,tt0
      write(*,*) 'Average over all aircraft profiles ',
     & ratio,' +/-',ratio*sqrt(tt0*tt2-tt1**2)/tt1
      write(lunw,*) 'Average over all aircraft profiles ',
     & ratio,' +/-',ratio*sqrt(tt0*tt2-tt1**2)/tt1
      close(lunw)
      stop
      end
