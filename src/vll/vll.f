c  Program VLL (View LineList)
c  For rapid viewing (and subsetting) of a HITRAN-format linelist
c  Assumes a linelist naming convention  xxxxx.nnn  where nnn is the number of characters per line
c     nnn = 101 (HITRAN_2000, Unix)
c     nnn = 102 (HITRAN_2000, DOS)
c     nnn = 161 (HITRAN_2004, Unix)
c     nnn = 162 (HITRAN_2004, DOS)
c  There must be only one period inthe file name.

      implicit none
      include "../ggg_int_params.f"

      integer*4 isot,molno,reclen,mgas,kgas,jgas,kiso,jiso,
     & lunr,nlps,jline,kline,lnbc,lloc,mchar,nline,posnall,
     $ ldot,ierr,fsib,file_size_in_bytes,lr
      parameter (lunr=14,mgas=75,nlps=18)
      real pbhw,pshift,sbhw,tdpbhw
      real*8 eprime,fpos,freq,stren
      character  ans*1,ccc*8,quantum*100,llformat*47,linfil*(mfilepath),
     $ molnam(mgas)*8,gggdir*(mpath), version*44,dl*1
c============================================================
      data linfil/' '/

      version=' VLL    Version 2.2.0    10-Aug-2011    GCT '
      write(*,*) version
      write(*,*)

      llformat=
     $'(i2,i1,f12.6,e10.3,10x,2f5.0,f10.4,f4.2,f8.6,a)'
      call get_ggg_environment(gggdir, dl)
      lr=lnbc(gggdir)
      ierr=0
      jgas=0
c
c  Read names of gases.
      open(lunr,file=gggdir(:lr)//'isotopologs'//dl//'isotopologs.dat')
      do while (jgas.lt.mgas .and. ierr.eq.0)
         read(lunr,'(1x,2i2,1x,a8)',iostat=ierr)jgas,jiso,molnam(jgas)
      end do
      close(lunr)
      write(*,*)' Reading isotopologs.dat which contains ',jgas,' gases'
      write(*,*)
      if(jgas.ge.mgas) write(*,*) 'Warning: Increase parameter MGAS'
c
c Prompt user for name of linelist
      do while (lnbc(linfil) .eq. 0) 
        write(6,*)'atm [a], gct [g], '//
     $  'hitran04 [4], hitran08 [8] other [o]'
        read(5,'(a1)') ans
        if(ans.eq.'a') linfil='atm.101'
        if(ans.eq.'g') linfil='gct.101'
        if(ans.eq.'f') linfil='fcia.101'
        if(ans.eq.'s') linfil='scia.101'
        if(ans.eq.'c') linfil='solar_dc.101'
        if(ans.eq.'i') linfil='solar_di.101'
        if(ans.eq.'0') linfil='hitran00.101'
        if(ans.eq.'4') linfil='hitran04.162'
        if(ans.eq.'8') linfil='hitran08.162'
        if(ans.eq.'o') then
           do while (lnbc(linfil).eq.0)
              write(6,'(a16,$)') 'other linelist: '
              read(5,'(a)')linfil
           end do
        endif
      end do  ! while (lnbc(linfil).eq.0) 
      linfil=gggdir(:lr)//'linelist'//dl//linfil(:lnbc(linfil))
      write(*,*)linfil

c Open the linelist and determine:
c   "reclen" from the filename extension.
c   "fsib" from subroutine file_size_in_bytes
c   "nline" by dividing "fsib" by "reclen".
      ldot=lloc(linfil,'.')   ! look for the period
      read(linfil(ldot+1:),'(i3)') reclen
c      write(*,*)ldot,linfil(ldot+1:),reclen
      fsib=file_size_in_bytes(lunr,linfil)
      open(lunr,file=linfil,access='direct',form='formatted',
     & status='old',recl=reclen)
      nline=fsib/reclen
      if(mod(fsib,nline) .ne. 0) then
         write(*,*)' Linelist size not divisible by reclen'
         write(*,*)linfil,fsib
         stop
      endif
      write(*,*)' Number of spectral lines       = ',nline
      write(*,*)' Number of characters per line  = ',reclen
      write(*,*)' Total number of characters     = ',fsib
      write(*,*)
c============================================================
      kgas=0
      kiso=0
      kline=1
      do                             ! G77
c      do while ( ccc(1:1) .ne. 'q')  ! sun
         write(6,*)' Enter frequency <cm-1>, continue <cr>,'//
     &    ' select gas <g#>, select isotope <i#>, or quit <q>'
         read(5,'(a)') ccc
         if(ccc(1:1).eq.'q') then
            close(lunr)
            stop
         elseif(ccc(1:1).eq.'g') then
            read(ccc(2:),*)kgas
         elseif(ccc(1:1).eq.'i') then
            read(ccc(2:),*)kiso
         else
            mchar=lnbc(ccc)
            if(mchar.gt.0) then   
               read(ccc,*)fpos  ! read real frequency entered
               ccc='        '
               if(fpos.ge.0)    ! if +ve position to frequency
     $         kline=posnall(lunr,fpos,nline)
            endif
         endif
c
         write(6,'(a)')
     $   'Gas     ID Iso  Frequency  Strength     E"     ABHW  SBHW'//
     $   '  TD  Shift'
         write(6,'(a)')'--------------------------------------------'//
     $   '----------------------------'
         jline=1
         do while (jline.le.nlps .and. kline.lt.nline)
            kline=kline+1
            read(lunr,llformat,rec=kline)molno,isot,freq,stren,
     &      pbhw,sbhw,eprime,tdpbhw,pshift,quantum
            if(index(linfil,'hitran').gt.0) then
               call hitran_to_atmos_gas_numbering(molno,isot)
            endif
            if(kgas.eq.0 .or. kgas.eq.molno) then 
               if(kiso.eq.0 .or. kiso.eq.isot) then 
                  jline=jline+1
                  write(6,'(a8,i2,i2,f13.6,1pe10.3,0pf10.4,2(1x,f5.4),
     &            1x,f3.2,2x,f5.3,2x,a)') molnam(molno),molno,
     &            isot,freq,stren,eprime,pbhw,sbhw,tdpbhw,pshift,
     &            quantum(:reclen-68)
               endif
            endif
         enddo
         if(kline.eq.nline) write(6,'(a)')' End of File'
         write(6,'(a)')'--------------------------------------------'//
     $   '----------------------------'
      end do
      end
