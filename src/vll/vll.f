c  Program VLL (View LineList)
c  For rapid viewing of a HITRAN-format linelist
c  Assumes a linelist naming convention  xxxxx.nnn
c  where nnn is the number of characters per line (typically 101)
      implicit none
      integer*4 isot,molno,reclen,mgas,kgas,jgas,kiso,jiso,
     & lunr,nlps,jline,kline,lnbc,lloc,nchar,nline,posnall,
     $ ldot,ierr,file_size,file_size_in_bytes,lr
      parameter (lunr=14,mgas=64,nlps=18)
      real pbhw,pshift,sbhw,tdpbhw
      real*8 eprime,fpos,freq,stren
      character  ans*1,ccc*8,fsw*3,llformat*55,linfil*80,
     $ molnam(mgas)*8,ref*6,root*48
c============================================================
      data linfil/' '/
      llformat=
     $'(i2,i1,f12.6,e10.3,10x,2f5.0,f10.4,f4.2,f8.6,24x,a3,a6)'
      call getenv('GGGPATH',root)
      lr=lnbc(root)
c       root='/home/toon/ggg'
      ierr=0
      jgas=0
      write(6,*)'View LineList   Version 2.0.1    23-Sep-06'
c
c  Read names of gases.
      open(lunr,file=root(:lr)//'/isotopologs/isotopologs.dat')
      do while (jgas.lt.mgas .and. ierr.eq.0)
         read(lunr,'(1x,2i2,1x,a8)',iostat=ierr)jgas,jiso,molnam(jgas)
c         write(*,*)jgas,molnam(jgas),llformat
      end do
      close(lunr)
      if(jgas.ge.mgas) write(*,*) 'Warning: Increase parameter MGAS'
c
c Prompt user for name of linelist
      do while (lnbc(linfil) .eq. 0) 
        write(6,*)'atm [a], gct [g], '//
     $  'hitran00 [0], hitran04 [4] other [o]'
        read(5,'(a1)') ans
        if(ans.eq.'a') linfil='atm.101'
        if(ans.eq.'g') linfil='gct.101'
        if(ans.eq.'s') linfil='sup.101'
        if(ans.eq.'x') linfil='solar.101'
        if(ans.eq.'0') linfil='hitran00.101'
        if(ans.eq.'4') linfil='hitran04.162'
        if(ans.eq.'o') then
           do while (lnbc(linfil).eq.0)
              write(6,'(a16,$)') 'other linelist: '
              read(5,'(a)')linfil
           end do
        endif
      end do  ! while (lnbc(linfil).eq.0) 
      linfil=root(:lr)//'/linelist/'//linfil(:lnbc(linfil))
      write(*,*)linfil

c Open the linelist.
c Determine "reclen" from the filename extension.
c Determin file_size from subroutine file_size_in_bytes
c Determine NLINE by dividing file_size by "reclen".
      ldot=lloc(linfil,'.')
      read(linfil(ldot+1:),'(i3)') reclen
c      write(*,*)ldot,linfil(ldot+1:),reclen
      file_size=file_size_in_bytes(lunr,linfil)
      open(lunr,file=linfil,access='direct',form='formatted',
     & status='old',recl=reclen)
      nline=file_size/reclen
      if(mod(file_size,nline) .ne. 0) then
         write(*,*)' Linelist size not divisible by reclen'
         write(*,*)linfil,file_size
         stop
      endif
      write(*,*)nline
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
            nchar=lnbc(ccc)
            if(nchar.gt.0) then   
               read(ccc,*)fpos  ! read real frequency entered
               ccc='        '
               if(fpos.ge.0)    ! if +ve position to frequency
     $         kline=posnall(lunr,fpos,nline)
            endif
         endif
c
         write(6,'(a)')
     $   'Gas     ID Iso  Frequency  Strength     E"     ABHW  SBHW'//
     $   '  TD  Shift  FSW  Refs.'
         write(6,'(a)')'--------------------------------------------'//
     $   '------------------------------------'
         jline=1
         do while (jline.le.nlps .and. kline.lt.nline)
            kline=kline+1
            read(lunr,llformat,rec=kline)molno,isot,freq,stren,
     &      pbhw,sbhw,eprime,tdpbhw,pshift,fsw,ref
            if(index(linfil,'hitran').gt.0) then
               call hitran_to_atmos_gas_numbering(molno,isot)
            endif
            if(kgas.eq.0 .or. kgas.eq.molno) then 
               if(kiso.eq.0 .or. kiso.eq.isot) then 
                  jline=jline+1
                  write(6,'(a8,i2,i2,f13.6,1pe10.3,0pf10.4,2(1x,f5.4),
     &            1x,f3.2,2x,f5.3,2x,a3,1x,a6)') molnam(molno),molno,
     &            isot,freq,stren,eprime,pbhw,sbhw,tdpbhw,pshift,fsw,ref
               endif
            endif
         enddo
         if(kline.eq.nline) write(6,'(a)')' End of File'
         write(6,'(a)')'--------------------------------------------'//
     $   '------------------------------------'
      end do
      end
