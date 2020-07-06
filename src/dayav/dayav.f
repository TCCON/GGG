c  Program dayav
c  Reads the output files (e.g. fts93avg.vav) produced  by GGGAVG,
c  then averages these results over each solar day.
c  Writes the averaged results.

      implicit none
      integer k,lunr,lunw,maux,naux,ncol,mcol,nrow,ngas,
     & nlhead,nlabel,lnbc,lr,
     & iyr,iywas,ldoy,ldwas,
     & ios,navg,ntotavg,nspe
      parameter (lunr=12)
      parameter (lunw=14)
      parameter (maux=25) ! number of auxiliary variables (e.g. Year, Day, Lat, Long)
      parameter (mcol=500)  ! maximum allowed number of primary variables
      character pabel*1250,outfile*40,infile*40,version*72,
c     & spectrum*35,specwas*35,
     & clabel(maux+2*mcol)*32,gtext*1200,rwformat*72,
     & chmiss*8,col1*1
      real*8 yaux(maux),yobs(mcol),yerr(mcol),
     & tt(mcol),ty(mcol),
     & ta(maux),time,twas,valmiss

      rwformat='format=(f14.8,nnf13.5,1000(1pe12.4))'
c====================================================================
c  Prompt user for input
      version=
     & ' DAYAV                    Version 8.01        2015-03-29   GCT'
      write(6,*)version
      lr=0
      do while (lr .le. 0)
         if (iargc() == 0) then
            write(6,'(a)') 'File to be averaged (e.g. gndallav.vav): '
            read(5,'(a)')infile
         elseif (iargc() == 1) then
            call getarg(1, infile)
         else
            stop 'Use: $gggpath/bin/dayav vavfile '
         endif

         lr=lnbc(infile)
      end do   !  while (lr .le. 0)
c====================================================================
c  Open input/output files and read/write header information.
      outfile=infile(:lr-1)//'d'
      open(lunr,file=infile,status='old')
      open(lunw,file=outfile,status='unknown')
      read(lunr,*)nlhead,ncol,nrow,naux
      if(ncol.gt.mcol) stop 'ncol > mcol'
      if(naux.gt.maux) stop 'naux > maux'
      write(lunw,*)nlhead+1,ncol
      write(lunw,'(a)') version
      do k=2,nlhead-3
         read(lunr,'(a)')gtext
         write(lunw,'(a)')gtext(:lnbc(gtext))
      end do
      read(lunr,'(a,1pe12.4)') chmiss,valmiss
      write(lunw,'(a,1pe12.4)') chmiss,valmiss
      read(lunr,'(a)')rwformat
      write(lunw,'(a)')rwformat
      read(lunr,'(a)')pabel
      write(lunw,'(a)')pabel(:lnbc(pabel)+1)
      call substr(pabel,clabel,naux+2*mcol,nlabel)
      if(nlabel.ne.ncol) then
         write(*,*)'nlabel,ncol=',nlabel,ncol
         stop 'nlabel.ne.ncol'
      endif
      ngas=(nlabel-naux)/2
      write(*,*) rwformat(8:)
c=================================================================
c  Read first data record and initialize
      read(lunr,rwformat(8:)) col1,(yaux(k),k=1,naux),
     & (yobs(k),yerr(k),k=1,ngas)
      iywas=int(yaux(1))
      ldwas=int(yaux(2)+yaux(6)/360.)  ! DOY + longit/360
      twas=365.25d0*yaux(1)  ! convert from years to days
      navg=1
      do k=1,naux
         ta(k)=yaux(k)
      end do
      do k=1,ngas
         if(yobs(k).ne.valmiss.and.yerr(k).gt.0.0)then
            tt(k)=(1.d0/yerr(k))**2
            ty(k)=(yobs(k)/yerr(k))/yerr(k)
         else
            tt(k)=0.0
            ty(k)=0.0
         endif
      end do
c=================================================================
c  Main loop
      nspe=0
      ntotavg=0
      ios=0
      write(*,*)'naux, ngas=',naux,ngas
      do while ( ios .eq. 0 )
         read(lunr,rwformat(8:),iostat=ios) col1,(yaux(k),k=1,naux),
     &   (yobs(k),yerr(k),k=1,ngas)
         time=365.25d0*yaux(1)  ! convert from years to days
         iyr=int(yaux(1))
         ldoy=int(yaux(2)+yaux(6)/360)
         write(61,*) yaux(1),yaux(2),iyr,iywas,ldoy,ldwas
         if( iyr.eq.iywas .and. ldoy.eq.ldwas .and. ios.eq.0 ) then  ! accumulate averages
            write(61,*)'same day'
            navg=navg+1
            do k=1,naux
               ta(k)=ta(k)+yaux(k)
            end do
            do k=1,ngas
               if(yobs(k).ne.valmiss.and.yerr(k).gt.0.0)then
                  tt(k)=tt(k)+(1.d0/yerr(k))**2
                  ty(k)=ty(k)+(yobs(k)/yerr(k))/yerr(k)
               endif
            end do
c            specwas=spectrum
         else  !  Compute averaged values, write then to output file, & reset TA, TT, TY.
c            if(twas.gt.0.0) then
            write(61,*)'different day'
            do k=1,ngas
               if( tt(k).gt.0 ) then
                  ty(k)=ty(k)/tt(k)            ! mean
                  tt(k)=1.d0/dsqrt(tt(k))      ! SE
               else
                  tt(k)=valmiss
                  ty(k)=valmiss
               endif
            end do
            write(lunw,rwformat(8:)) ' ',(ta(k)/navg,k=1,naux),
     &        (ty(k),tt(k),k=1,ngas)
            nspe=nspe+navg
            ntotavg=ntotavg+1
c            endif     !  if(twas.gt.0.0) then
c  Re-set TA, TT, TY to latest value.
            navg=1
            twas=time
            ldwas=ldoy
            iywas=iyr
            do k=1,naux
               ta(k)=yaux(k)
            end do
            do k=1,ngas
               if(yobs(k).ne.valmiss.and.yerr(k).gt.0.0)then
                  tt(k)=(1.d0/yerr(k))**2
                  ty(k)=(yobs(k)/yerr(k))/yerr(k)
               else
                  tt(k)=0.0
                  ty(k)=0.0
               endif
            end do
         endif   ! if(iyr.eq.iywas .and. ldoy.eq.ldwas .and. ios.eq.0)
      end do   !  while ( ios .eq. 0 )
      write(*,*)' Original number of input data records = ',nspe
      write(*,*)' Resulting number of averaged records  = ',ntotavg
      close(lunr)
      close(lunw)
      stop
      end
