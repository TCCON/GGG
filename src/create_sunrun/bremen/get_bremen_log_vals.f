      subroutine get_bremen_log_vals(specpath,iend,iy,
     & im,id,hh,mm,ss,ms,dur,tout,pout,hout)
!     NMD 101101
!     subroutine to return the matching log values for
!     temp,press and humidity for each spectrum

!     Inputs:
!     Path of spectrum (specpath)
!     Endian-ness (iend)
!     Date/Time info
!     Spectrum duration

!     Outputs:
!     Temp, Press, Hum, averaged over the spectrum

      implicit none

!     Declare variables
      character*38 specpath
      integer*4 iend
      integer*4 iy,im,id,hh,mm,ss,ms
      real*8 dur
      real*8 tout,pout,hout
    
      integer(4) one,lulog,fhlog
      character*128 logpath, fhlogpath
      character*160 header

      parameter(one=1,lulog=21,fhlog=22)

      character logfile*80,logdate*10,logtime*8,root*64,dl*1
      character fhlogdate*18
      integer*4  j,jj,jd,doy,nrpts,nlogvals,lnbc,lr,lrt
      integer*4 logmon(1440),logday(1440),logyear(1440),
     & loghour(1440),logmin(1440),logsec(1440),ierr,last,first
      integer*4 fhlogyear(999999),fhlogmon(999999),
     & fhlogday(999999),fhloghour(999999),fhlogmin(999999)
      real*8 logpress(1440),logtemp(1440),loghum(1440),
     & r8spectime,r8logtime(1440),r8specfin,avg_p,avg_t,avg_h,
     & tot_time,difftime,fhlogpress(999999),fhlogtemp(999999),
     & r8fhlogtime(999999)

!     Begin
      dl='/'
      call getenv('GGGPATH',root)
      root=root(:lnbc(root))//dl
      lrt=lnbc(root)
    
      lr=0

!      write(*,*) specpath(1:2)

!     We need to create the start and finish times in the log to check against
      if(specpath(1:2).eq.'br') then
        logfile=specpath(3:10)//'.log'
      elseif(iy.gt.1999) then
c        write(*,*) specpath(1:6)
        logfile='20'//specpath(1:6)//'.log'
      else
        logfile='19'//specpath(1:6)//'.log'
      endif
      
      call julian(iy,im,id,jd)
      call julian(iy,one,one,jj)
      doy=jd-jj+1

!      logpath='/home/home/n_deutscher/ggg/bremen_ptu/'//logfile
      logpath=root(:lrt)//'src/create_sunrun/bremen/bremen_ptu/'//
     & logfile
!      fhlogpath='/home/home/n_deutscher/ggg/bremen_ptu/'
      fhlogpath=root(:lrt)//'src/create_sunrun/bremen/bremen_ptu/'
     &  //'br_flughafen_2003-9.txt'
      r8spectime=iy+(doy+hh/24.0d0+mm/(24.0d0*60.0d0)+
     & ss/(24.0d0*3600.0d0))/366.0d0
      r8specfin=r8spectime+dur/(24.0d0*3600.0d0*366.0d0)
        
      nrpts=0
      nlogvals=0
      loghour=0
      logyear=0

!     Initialise Pout, Hout, Tout to 9999.0, 0.0, 0.0
      Pout=99.0
      Hout=0.0
      Tout=0.0

      open(LUlog,file=logpath,status='old',iostat=ierr)

      if(ierr.ne.0)then
        write(6,*) 'file not found: ', trim(logpath)
        write(6,*) 'attempting to use airport P/T data'
        open(fhlog,file=fhlogpath,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(6,*) 'error opening: ',trim(fhlogpath)
          return
        endif
        goto 95
      endif

      read(Lulog,*,err=97) header
      do j=1,2440
        read(Lulog,*,end=98) logdate,logtime,logpress(j-nrpts),
     &   logtemp(j-nrpts),loghum(j-nrpts)
        read(logtime(1:2),'(i2)') loghour(j-nrpts)
        read(logtime(4:5),'(i2)') logmin(j-nrpts)
        read(logtime(7:8),'(i2)') logsec(j-nrpts)
        read(logdate(1:2),'(i2)') logmon(j-nrpts)
        read(logdate(4:5),'(i2)') logday(j-nrpts)
        read(logdate(7:10),'(i4)') logyear(j-nrpts)

        logsec(j-nrpts)=logsec(j-nrpts)+30

        call julian(logyear(j-nrpts),logmon(j-nrpts),
     &   logday(j-nrpts),jd)
        call julian(logyear(j-nrpts),one,one,jj)
        doy = jd-jj + 1
        r8logtime(j-nrpts)=logyear(j-nrpts)+(doy+
     &   (loghour(j-nrpts)+(logmin(j-nrpts)+
     &    logsec(j-nrpts)/60.0d0)/60.0d0)/24.0d0)/366.0d0

!       Check for duplicate entries. Crude at the moment.
        if(j.gt.1)then
!          if(logmin(j-nrpts).le.logmin(j-nrpts-1)) then
!            if(loghour(j-nrpts).le.loghour(j-nrpts-1)) then
!              if(logday(j-nrpts).le.logday(j-nrpts-1)) then
!                if(logmon(j-nrpts).le.logmon(j-nrpts-1)) then
!                  if(logyear(j-nrpts).le.logyear(j-nrpts-1)) then
           if(r8logtime(j-nrpts).le.r8logtime(j-nrpts-1)) then
!                    write(6,*) loghour(j-nrpts),logmin(j-nrpts)
                    nrpts=nrpts+1
                  else
                    nlogvals=nlogvals+1
                  endif
!                else  
!                  nlogvals=nlogvals+1
!                endif
!              else
!                nlogvals=nlogvals+1
!              endif
!            else
!              nlogvals=nlogvals+1
!            endif
!          else
!            nlogvals=nlogvals+1
!          endif
        else
          nlogvals=nlogvals+1
        endif

!        write(16,*) logdate,' ',logtime,' ',logyear(j),
!     &   logmon(j),logday(j),loghour(j),logmin(j),
!     &   logsec(j),logpress(j),logtemp(j),loghum(j)

      enddo

98    close(Lulog)

      tot_time=0
      avg_p=0
      avg_h=0
      avg_t=0

!     Now check how the time matches with the spectrum time.
!        write(6,*) nlogvals, 'nlogvals'
      do j=1,nlogvals

!       because log times are instantaneous we add 30secs
!       to counter the inherent assumption that the data
!       are representative of the previous minute
!        write(6,*) logyear(j),r8logtime(j),r8spectime,r8specfin
        if(r8logtime(j).gt.r8spectime .and. 
     &   r8logtime(j).lt.r8specfin) then
          if(j.gt.1)then
            if(r8logtime(j-1).lt.r8spectime) then
!             spectrum doesn't cover the whole minute of the log entry
              difftime=(r8logtime(j)-r8spectime)*365.25*24*3600
              avg_p=avg_p+logpress(j)*difftime
              avg_t=avg_t+logtemp(j)*difftime
              avg_h=avg_h+loghum(j)*difftime
              tot_time=tot_time+difftime
            elseif(r8logtime(j-1).ge.r8spectime) then
              difftime=(r8logtime(j)-r8logtime(j-1))*365.25
     &          *24*3600
              avg_p=avg_p+logpress(j)*difftime
              avg_t=avg_t+logtemp(j)*difftime
              avg_h=avg_h+loghum(j)*difftime
              tot_time=tot_time+difftime
            else
              write(6,*) 'why the hell are we here?'
            endif
          endif
        elseif(r8logtime(j).gt.r8specfin)then
          if(j.gt.1)then
            if(r8logtime(j-1).lt.r8specfin)then
              difftime=(r8specfin-r8logtime(j-1))*365.25
     &          *24*3600
              avg_p=avg_p+logpress(j)*difftime
              avg_t=avg_t+logtemp(j)*difftime
              avg_h=avg_h+loghum(j)*difftime
              tot_time=tot_time+difftime
            endif 
          else
            write(6,*) 'no entries in log file for this spectrum'
            write(6,*) 'attempting to use airport P/T data'
            open(fhlog,file=fhlogpath,status='old',iostat=ierr)
            if(ierr.ne.0)then
              write(6,*) 'error opening: ',trim(fhlogpath)
              return
            endif
            goto 95
          endif
        endif
      enddo

      if(tot_time.gt.0)then
        Pout=avg_p/tot_time
        Tout=avg_t/tot_time
        Hout=avg_h/tot_time
      endif

      return
        
95    do j=1,999999
        read(FHlog,*,end=96) fhlogdate,fhlogtemp(j),fhlogpress(j)

        read(fhlogdate(6:9),'(i4)') fhlogyear(j)
        read(fhlogdate(10:11),'(i2)') fhlogmon(j)
        read(fhlogdate(12:13),'(i2)') fhlogday(j)
        read(fhlogdate(14:15),'(i2)') fhloghour(j)
        read(fhlogdate(16:17),'(i2)') fhlogmin(j)

        call julian(fhlogyear(j),fhlogmon(j),fhlogday(j),jd)
        call julian(fhlogyear(j),one,one,jj)
        doy= jd-jj+1
        r8fhlogtime(j)=fhlogyear(j)+(doy+(fhloghour(j)-1+
     &   (fhlogmin(j))/60.0d0)/24.0d0)/366.0d0

c        write(6,*) r8fhlogtime(j),fhlogdate,fhlogtemp(j),
c     &   fhlogpress(j)

c   this line is here to try to speed things up a little. 
        if(r8fhlogtime(j).gt.(r8spectime+1/(60.0d0*24.0d0*
     &   366.0d0))) then
          goto 96
        endif
      enddo

96    close(FHlog)  

c    Now match the airport temp/press to the correct spectrum time.
c    To do this, find the MET readings spanning the measurement time.
c    If it turns out that the spectrum is wholly bracketed by two MET
c    readings, then linearly interpolate to the average spectrum time.
c    Otherwise, weight the spectrum to how much of it is contained 
c    within each MET interval.

c    Before we start that, do some simple checks:
      if(r8spectime+10/(60.0d0*24.0d0*366.0d0).lt.r8fhlogtime(1))
     &  then
        write(6,*) 'no entry for this spectrum'
        return
      endif

      do j=1,999999
        if(r8spectime.gt.r8fhlogtime(j)) then
          first=j
        endif
        if(r8specfin.lt.r8fhlogtime(j)) then
          last=j
        goto 94
        endif
      enddo
      
94    if(last.eq.first+1) then
        Pout=(fhlogpress(last)-fhlogpress(first))*
     &   (r8spectime-r8fhlogtime(first))/(r8fhlogtime(last)
     &   -r8fhlogtime(first))+fhlogpress(first)
        Pout=Pout/10.0d0
        Tout=(fhlogtemp(last)-fhlogtemp(first))*
     &   (r8spectime-r8fhlogtime(first))/(r8fhlogtime(last)
     &   -r8fhlogtime(first))+fhlogtemp(first)
        Tout=Tout/10.0d0
        Hout=0
      elseif(last.eq.first+2) then
        Pout=(fhlogpress(last)-fhlogpress(first))*
     &   (r8spectime-r8fhlogtime(first))/(r8fhlogtime(last)
     &   -r8fhlogtime(first))+fhlogpress(first)
        Pout=Pout/10.0d0
        Tout=(fhlogtemp(last)-fhlogtemp(first))*
     &   (r8spectime-r8fhlogtime(first))/(r8fhlogtime(last)
     &   -r8fhlogtime(first))+fhlogtemp(first)
        Tout=Tout/10.0d0
        Hout=0
c        Pout=
c        Tout=
c        Hout=0
      else
        write(6,*) 'what on earth!'
c       we've already initialised these values, no need to set them again
c        Pout=99.0
c        Hout=0.0
c        Tout=0.0
      endif

c       conversion because of altitude difference
        Pout=Pout*exp(-9.81*(27-4)/(287*(Tout+273.15)))

c      write(6,*) Pout,Tout,Hout
      return

97    write(6,*) 'log file ',logpath, ' not found'
      return

      end
