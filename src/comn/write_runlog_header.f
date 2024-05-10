      subroutine write_runlog_header(lunw_rlg,version,data_fmt_wrlg)
      integer*4 lunw_rlg,mcol,ncol,ldfr,lnblnk
      parameter (mcol=40)
      character
     &  col_labels*320,outarr(mcol)*40,data_fmt_wrlg*(*),version*(*)

c  Inputs:
c    lunw_rlg  I*4
c    version   C*(*)
c    
c Outputs:
c    data_fmt_wrlg   C*(*)
c
c  Usage:
c    open(lunw_rlg,file='runlog_name')
c    call write_runlog_header(lunw_rlg,header,data_fmt_wrlg)
c    do ispec=1,nspec
c       call write)runlog_data_record(lunw_rlg,data_fmt_wrlg, ....)
c    end do
c    close(lunw_rlg)

      col_labels=
     &  '    Spectrum_File_Name                                  '//
     &  '   Year  Day  Hour   oblat    oblon   obalt    ASZA'//
     &  '   POFF    AZIM   OSDS    OPD   FOVI  FOVO'//
     &  '  AMAL   IFIRST    ILAST     DELTA_NU    POINTER BPW ZOFF'//
     &  '  SNR APF  tins  pins  hins   tout   pout  hout'//
     &  '  sia    fvsi   wspd  wdir  lasf    wavtkr  aipl'

      call substr(col_labels, outarr, mcol, ncol)
      if(ncol.gt.mcol) stop 'write_runlog_header: Increase MCOL'

      data_fmt_wrlg='(a1,a57,1x,2i4,f8.4,f8.3,f9.3,2f8.3,1x,f6.4,'//
     & 'f8.3,f7.3,f7.2,3(1x,f5.4),2i9,1x,f14.11,i9,i3,1x,f5.3,i5,1x,'//
     & 'a2,2(f6.1,f8.2,f5.1),f7.1,f7.4,f6.1,f6.0,f10.3,f7.0,f7.3)'
      ldfr=lnblnk(data_fmt_wrlg)

      write(lunw_rlg,*)4,ncol
      write(lunw_rlg,'(a)') version
      write(lunw_rlg,'(a7,a)')'format=',data_fmt_wrlg(:ldfr)
      write(lunw_rlg,'(a)') col_labels(:lnblnk(col_labels))
      return
      end
