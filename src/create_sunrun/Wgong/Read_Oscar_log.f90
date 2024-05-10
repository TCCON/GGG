    subroutine Read_Oscar_log (specname,specpath,logname,logpath,Pout,Tout,SIA,SIS,foundinlog) 

!   Read OSCAR log to retrieve pressure, temperature sis and sia.
!   DG Nov 07, updated and tidied April 2009


    implicit none
    character   specname*(*), logname*(*), specpath*(*), logpath*(*)
    integer*4   i,j,luns, Tcomma, Pcomma, SIAcomma, SIScomma
    integer*4   idate, ierror, lr, idum, lnblnk
    logical*1   foundinlog, foundinheader
    parameter   (luns=19)

    real(8)::   Pout,Tout,SIA, SIS

    character*256   string
    character*256   searchname
    character*1     dummy

    idum=lnblnk(specpath)  ! Avoid compiler warning (unused)
    idum=lnblnk(logname)  ! Avoid compiler warning (unused)
    open(luns,file=logpath, status='old', iostat=ierror)
    if(ierror.ne.0) then
        print *, 'Read_Oscar_log: Log file not found ', trim(logpath)
        return
    endif

!   read the first 10 lines
!   The 10th line contains the column headings
    do i=1,10
        read(luns,'(a)')string
    end do

!   Find the columns headed Troom, Pbarametric, SIA and SIS
      Tcomma=0
      Pcomma=0
      SIAcomma=0
      SIScomma=0
    do i=1,len(trim(string))
        if(string(i:i).eq.',')Tcomma=Tcomma+1
        if(string(i+1:i+5).eq.'Troom')exit
    end do
    do i=1,len(trim(string))
        if(string(i:i).eq.',')Pcomma=Pcomma+1
        if(string(i+1:i+4).eq.'Pbar')exit
    end do
    
    foundinheader=.false.
    do i=1,len(trim(string))
        if(string(i:i).eq.',')SIAcomma=SIAcomma+1
        if(string(i+1:i+3).eq.'SIA') then
            foundinheader=.true.
            exit
        endif
     end do
     if(.not.foundinheader)SIAcomma=0

    foundinheader=.false.
    do i=1,len(trim(string))
        if(string(i:i).eq.',')SIScomma=SIScomma+1
            if(string(i+1:i+3).eq.'SIS')then
            foundinheader=.true.
            exit
        endif
     end do
     if(.not.foundinheader)SIScomma=0

    lr = len(trim(specname))
    if(specname(lr:lr).eq.'1') then
      searchname=specname(:(lr-1))//'0'
    else
      searchname=trim(specname)
    endif

!   Read through log file and find line for specname/searchname
    do
        ierror=0
        read(luns,'(a)',iostat=ierror) string
        if(ierror<0)then        !EOF
            Pout=0.0
            Tout=0.0
            SIS=0.0
            SIA=0.0
            print *, 'Read_Oscar_log: Spectrum ', trim(searchname), ' not found in log'
            foundinlog=.false.
            return
        endif
        if(index(string,trim(searchname)).ne.0)exit
!        j=index(searchname,"\",.true.)
!        if(index(string,trim(searchname(j+1:))).ne.0)exit
       enddo

!   read the required values from this line
    foundinlog=.true.
    j=0
    do i=1,len(trim(string))
        read(string(i:i),'(a1)')dummy
        if(dummy.eq.',')j=j+1
        if(j.eq.Tcomma)exit
    end do
    read(string(i+1:),*)Tout

    j=0
    do i=1,len(trim(string))
        read(string(i:i),'(a1)')dummy
        if(dummy.eq.',')j=j+1
        if(j.eq.Pcomma)exit
    end do
    read(string(i+1:),*)Pout

!   Pressure calibration: 
!   calibrated 01 July 2008 => correction = -0.8mb
!    j=index(searchname,"\",.true.)
!    read(specpath(j+5:j+10),'(I6)')idate
    read(specname(5:10),'(I6)')idate
    if(idate < 080702)Pout=Pout-0.8
!   End pressure calibration

    if(SIAcomma/=0)then
        j=0
        do i=1,len(trim(string))
            read(string(i:i),'(a1)')dummy
            if(dummy.eq.',')j=j+1
            if(j.eq.SIAcomma)exit
        end do
        read(string(i+1:),*)SIA
    else
        SIA=0
    endif

    if(SIScomma/=0)then
        j=0
        do i=1,len(trim(string))
            read(string(i:i),'(a1)')dummy
            if(dummy.eq.',')j=j+1
            if(j.eq.SIScomma)exit
        end do
        read(string(i+1:),*)SIS
    else
    SIS=0
    endif

!   Convert SIS (now in %) to absolute units to be compatible with ipp and ggg
    sis=sis*sia/100.

    return
    end
