c  Program dayav
c  Reads the output files (e.g. fts93avg.vav) produced  by GGGAVG,
c  then averages these results over a user-selected period (AVPER days)
c  and writes the averaged results.

      implicit none
      integer k,lunr,lunw,maux,naux,ncol,mcol,nrow,ngas,
     & nlhead,nlabel,lnbc,lr,
     & ios,navg,ntotavg,nspe
      parameter (lunr=12)
      parameter (lunw=14)
      parameter (maux=25) ! number of auxiliary variables (e.g. Year, Day, Lat, Long)
      parameter (mcol=500)  ! maximum allowed number of primary variables
      character pabel*1000,outfile*40,infile*40,version*44,keyword*8,
     & spectrum*35,specwas*35,
     & clabel(maux+2*mcol)*32,gtext*1200
      real*8 yaux(maux),yobs(mcol),yerr(mcol),
     & tt(mcol),ty(mcol),t2(mcol),
     & ta(maux),avper,time,twas,valmiss
c====================================================================
c  Prompt user for input
      version=' DAYAV    Version 7.6.0    17-Nov-2009   GCT'
      write(6,*)version
      lr=0
      do while (lr .le. 0)
         write(6,'(a)') ' File to be averaged (e.g. gndallav.vav): '
         read(5,85)infile
         lr=lnbc(infile)
      end do   !  while (lr .le. 0)
      open(lunr,file=infile,status='old')
 85   format(a)
c
      write(6,'(a)')' Averaging period in days (e.g. 30.5 = 1 month):'
      read(5,*) avper
c====================================================================
c  Open input & output files and read/write header information.
      outfile=infile(:lr-1)//'d'
      open(lunw,file=outfile,status='unknown')
      read(lunr,*)nlhead,ncol,nrow,naux
      if(ncol.gt.mcol) stop 'ncol > mcol'
      if(naux.gt.maux) stop 'naux > maux'
      write(lunw,*)nlhead+1,ncol
      write(lunw,'(a,f11.6)')version//'  Averaging period(days)=',avper
      do k=2,nlhead-2
        read(lunr,'(a)')gtext
        write(lunw,'(a)')gtext(:lnbc(gtext))
      end do
      read(lunr,'(a8,1pe12.4)') keyword,valmiss
      write(lunw,'(a8,1pe12.4)') keyword,valmiss
      read(lunr,85)pabel
      write(lunw,85)pabel(:lnbc(pabel)+1)
      call substr(pabel,clabel,naux+2*mcol,nlabel)
      ngas=(nlabel-naux)/2
c=================================================================
c  Main loop
      twas=0.0
      nspe=0
      ntotavg=0
      ios=0
      write(*,*)'naux, ngas=',naux,ngas
      do while ( ios .eq. 0 )
         read(lunr,75,iostat=ios) spectrum, (yaux(k),k=1,naux),
     &   (yobs(k),yerr(k),k=1,ngas)
 75      format(a35,f14.8,23f13.5,1000(1pe12.4))
         time=365.25d0*yaux(1)  ! convert from years to days
         if ( dabs(time-twas) .le. avper .and. ios .eq. 0 ) then
            navg=navg+1
            do k=1,naux
               ta(k)=ta(k)+yaux(k)
            end do
            do k=1,ngas
               if(yobs(k).ne.valmiss.and.yerr(k).gt.0.0)then
                  tt(k)=tt(k)+1.d0/yerr(k)**2
                  ty(k)=ty(k)+yobs(k)/yerr(k)**2
                  t2(k)=t2(k)+(yobs(k)/yerr(k))**2
               endif
            end do
            specwas=spectrum
         else  !  Write averaged values to output file & reset TA, TT, TY.
            if(twas.gt.0.0) then
            do k=1,ngas
               if( tt(k).gt.0 ) then
                  t2(k)=t2(k)/tt(k)            
                  ty(k)=ty(k)/tt(k)            ! mean
                  tt(k)=1.d0/dsqrt(tt(k))      ! SE
                  if(t2(k).ge.ty(k)**2) then
                    t2(k)=sqrt(t2(k)-ty(k)**2)   ! SD
                  else
                    t2(k)=0.0
                  endif
               else
                  tt(k)=valmiss
                  ty(k)=valmiss
                  t2(k)=valmiss
               endif
            end do
            ta(4)=navg**2
            write(lunw,75)specwas,(ta(k)/navg,k=1,naux),
     &        (ty(k),t2(k),k=1,ngas)
            nspe=nspe+navg
            ntotavg=ntotavg+1
c            write(*,'(a35,3i6,3f8.5)')specwas,
c     &      int(ta(1)/navg),nint(ta(2)/navg),navg,
c     &      t2(4)/ty(4),t2(5)/ty(5),ty(7)
            write(34,'(a25,3i6,4f8.5)')specwas(:25),
     &      int(ta(1)/navg),nint(ta(2)/navg),navg,
     &      t2(10)/ty(10), ! HF
     &      t2(4)/ty(4),   ! CH4
     &      t2(5)/ty(5),   ! CO2
     &      ty(7)          ! H2O
            endif
c  Set TA, TT, TY to zero.
            navg=1
            do k=1,naux
               ta(k)=yaux(k)
            end do
            do k=1,ngas
               if(yobs(k).ne.valmiss.and.yerr(k).gt.0.0)then
                  tt(k)=1.d0/yerr(k)**2
                  ty(k)=yobs(k)/yerr(k)**2
                  t2(k)=t2(k)+(yobs(k)/yerr(k))**2
               else
                  tt(k)=0.0
                  ty(k)=0.0
                  t2(k)=0.0
               endif
            end do
            twas=time
         endif
      end do   !  while ( ios .eq. 0 )
      write(*,*)' Original number of input data records = ',nspe
      write(*,*)' Resulting number of averaged records  = ',ntotavg
      close(lunr)
      close(lunw)
      stop
      end
