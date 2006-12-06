c  Program reads the output files (e.g. fts93avg.vav) produced  by GGGAVG,
c  then averages these results over a user-selected period (AVPER days)
c  and writes the averaged results.
      integer j,k,lunr,lunw,naux,ncol,mcol,nlhead,nfmt,nlabel,lnbc,lr,
     & ios,navg,ntotavg,nspe
      parameter (lunr=12)
      parameter (lunw=14)
      parameter (naux=15)   ! number of auxiliary variables (e.g. Year, Day, Lat, Long)
      parameter (mcol=500)  ! maximum allowed number of primary variables
      character pabel*1000,outfile*40,infile*40,version*44,keyword*8,
     & clabel(naux+2*mcol)*32,gtext*1200
      real*8 yaux(naux),yobs(mcol),yerr(mcol),tt(mcol),ty(mcol),
     & ta(naux),avper,time,twas,valmiss(naux+2*mcol)
c====================================================================
c  Prompt user for input
      version=' DAYAV    Version 7.5.1    21-Jun-2004   GCT'
      write(6,*)version
      lr=0
      do while (lr .le. 0)
         write(6,'(a,$)') ' File to be averaged (e.g. gndallav.vav): '
         read(5,85)infile
         lr=lnbc(infile)
      end do   !  while (lr .le. 0)
      open(lunr,file=infile,status='old')
 85   format(a)
c
      write(6,'(a,$)')' Averaging period in days (e.g. 30.5 = 1 month):'
      read(5,*) avper
c====================================================================
c  Open input & output files and read/write header information.
      outfile=infile(:lr-1)//'d'
      open(lunw,file=outfile,status='unknown')
      read(lunr,*)nlhead,nfmt
      write(lunw,*)nlhead+1,nfmt
      write(lunw,'(a,f11.6)')version//'  Averaging period(days)=',avper
      do k=2,nlhead-1
        read(lunr,'(a)')gtext
        write(lunw,'(a)')gtext(:lnbc(gtext))
      end do
      if ( gtext(:8) .eq. 'MISSING:' ) then
        read(gtext,9000) keyword,(valmiss(k),k=1,nfmt)
9000    format(a8,1200(1pe12.4))
      else
        do j=1,nfmt
           valmiss(j)=9.9999E+29
        end do
      endif
      read(lunr,85)pabel
      write(lunw,85)pabel(:lnbc(pabel)+1)
      call substr(pabel,clabel,naux+2*mcol,nlabel)
      ncol=(nlabel-naux)/2
      if(ncol.gt.mcol) then
        write(*,*)' Increase parameter MCOL to ',ncol
        stop
      endif
c=================================================================
c  Main loop
      twas=0.0
      nspe=0
      ntotavg=0
      ios=0
      do while ( ios .eq. 0 )
         read(lunr,75,iostat=ios) (yaux(k),k=1,naux),
     &   (yobs(k),yerr(k),k=1,ncol)
 75      format(f12.6,14f12.5,1000(1pe12.4))
         time=365.25d0*yaux(1)
c         time=365.25d0*yaux(1)+yaux(2)+yaux(4)/24  ! old format files
         if ( dabs(time-twas) .le. avper .and. ios .eq. 0 ) then
            navg=navg+1
            do k=1,naux
               ta(k)=ta(k)+yaux(k)
            end do
            do k=1,ncol
               if(yobs(k).ne.valmiss(naux+2*k-1).and.yerr(k).gt.0.0)then
                  tt(k)=tt(k)+1.d0/yerr(k)**2
                  ty(k)=ty(k)+yobs(k)/yerr(k)**2
               endif
            end do
         else  !  Write averaged values to output file & reset TA, TT, TY.
            if(twas.gt.0.0) then
            do k=1,ncol
               if( tt(k).gt.0 ) then
                  ty(k)=ty(k)/tt(k)
                  tt(k)=1.d0/dsqrt(tt(k))
               else
                  tt(k)=valmiss(naux+2*k-1)
                  ty(k)=valmiss(naux+2*k-1)
               endif
            end do
            write(lunw,75)(ta(k)/navg,k=1,naux),(ty(k),tt(k),k=1,ncol)
            nspe=nspe+navg
            ntotavg=ntotavg+1
            write(*,'(2i6,2f7.0)')ntotavg,navg,ta(1)/navg,ta(2)/navg
            endif
c  Set TA, TT, TY to zero.
            navg=1
            do k=1,naux
               ta(k)=yaux(k)
            end do
            do k=1,ncol
               if(yobs(k).ne.valmiss(naux+2*k-1).and.yerr(k).gt.0.0)then
                  tt(k)=1.d0/yerr(k)**2
                  ty(k)=yobs(k)/yerr(k)**2
               else
                  tt(k)=0.0
                  ty(k)=0.0
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
