      subroutine write_runlog_header(lunw_rl,version,data_fmt_rl)
      integer*4 lunw_rl,mcol,ncol
      parameter (mcol=40)
      character col_labels*320,outarr(mcol)*40,
     & data_fmt_rl*(*),version*(*)
c
c  Inputs:
c    lunw_rl  I*4
c    version   C*(*)
c    
c Outputs:
c    data_fmt_rl   C*(*)
c
c  Usage:
c    open(lunw_rl,file='runlog_name')
c    call write_runlog_header(lunw_rl,header,data_fmt_rl)
c    do ispec=1,nspec
c       call write)runlog_data_record(lunw_rl,data_fmt_rl, ....)
c    end do
c    close(lunw_rl)

      col_labels=
     &  '    Spectrum_File_Name                                  '//
     &  '   Year  Day  Hour   oblat    oblon   obalt    ASZA'//
     &  '   POFF    AZIM   OSDS    OPD   FOVI  FOVO'//
     &  '  AMAL   IFIRST    ILAST     DELTA_NU    POINTER BPW ZOFF'//
     &  '  SNR APF  tins  pins  hins   tout   pout  hout'//
     &  '  sia    fvsi   wspd  wdir  lasf    wavtkr  aipl'

      call substr(col_labels, outarr, mcol, ncol)
      if(ncol.gt.mcol) stop 'write_runlog_header: Increase MCOL'

      data_fmt_rl='(a1,a57,1x,2i4,f8.4,f8.3,f9.3,2f8.3,1x,f6.4,f8.3,'//
     & 'f7.3,f7.2,3(1x,f5.4),2i9,1x,f14.11,i9,i3,1x,f5.3,i5,1x,'//
     & 'a2,2(f6.1,f8.2,f5.1),f7.1,f7.4,f6.1,f6.0,f10.3,f7.0,f7.3)'

      write(lunw_rl,*)4,ncol
      write(lunw_rl,'(a)') version
      write(lunw_rl,'(a7,a)')'format=',data_fmt_rl(:lnbc(data_fmt_rl))
      write(lunw_rl,'(a)') col_labels(:lnbc(col_labels))
      return
      end
