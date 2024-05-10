c  Program spectrum_average.f

c  Averages a series of binary spectra.

c  Works for any binary spectra having a GGG-format runlog,
c  which is adapted for use as the input file (including header).
c  [New extra lines are inserted containing the string "average"
c  to delimit the groups of spectra to be averaged.]
c
c  Creates separate averages for interleaved HgCd/InSb or InGaAs/Si spectra
c
c  Works transparently for I*2 and R*4 binary data types.
c
c  Byte-reversal handled automatically, if the computer
c  that you are working on has a different endian-ness
c  from the computer that wrote the binary spectra.
c
c  Headerless binary average output spectra are created
c  in the local directory.
c
      implicit none

      integer*4 lst,lr,
     & lunr_ss_rl,     ! LUN to read input runlogs from
     & lun_rpt,     ! LUN to write rms deviation
     & lunw_avg_rl,     ! LUN to write average spectrum to
     & lnbc,        ! function Last Non-Black Character
     & lloc,        ! function Last Location Of Character
     & ntype,jtype, ! Periodicity of spectra to be averaged
     & nspe         ! number of single spectra being averaged

      parameter (lun_rpt=30,lunr_ss_rl=15,lunw_avg_rl=16)

      character 
     & version*56,   ! version
     & rptfile*40,
     & data_fmt_read_rl*256,
     & data_fmt_write_rl*256,
     & stuff*256,cdum*12,
c     & menuinputfile*40,
     & runlog*120,col_labels_rl*512 

      integer*4
     & krec,
     & ldot,
     & istat         ! status flag (0=success, 1=EOF)

      real*8
     & nus, nue      ! selected frequency range of interest

      version=' spectrum_averaging    Version 0.35    2018-04-05   GCT '
      write(6,*) version

      ldot=0 ! avoid compiler warning (may be used uninitialized)

      if (iargc() == 0) then
         write(*,*)'Enter path to input_file/runlog:'
         read(*,'(a)') runlog

         write(*,*)'Enter Starting & Ending frequencies (cm-1):'
         write(*,*)'Enter 0 99999 to retain original spectral limits'
         read(*,*) nus,nue

         write(*,*)
     &    'Enter number of spectrum types interleaved in runlog:'
         write(*,*)'(e.g., HgCd/InSb=2; InGaAs/Si=2; InGaAs-only=1)'
         read(*,*) ntype

      elseif (iargc() == 1) then
         call getarg(1,runlog)
         nus=0.0d0
         nue=999999.9d0
         ntype=2
      elseif (iargc() == 4) then
         call getarg(1,runlog)
         call getarg(2,cdum)
         read(cdum,*) nus
         call getarg(3,cdum)
         read(cdum,*) nue
         call getarg(4,cdum)
         read(cdum,*) ntype
      else
         write(*,*)'Usage: $gggpath/bin/spectrum_averaging runlog'//
     & '   nus  nue   ntype'
         write(*,*)'(ntype= HgCd/InSb=2; InGaAs/Si=2; InGaAs-only=1)'
         stop
      endif
      write(*,*) 'Assuming nus= ',sngl(nus),' nue= ',sngl(nue),
     & 'ntype= ',ntype

      open(lunr_ss_rl,file=runlog,status='old')
c      read(lunr_ss_rl,*) nhl,ncol
c      call skiprec(lunr_ss_rl,nhl-2)
c      read(lunr_ss_rl,'(a)') col_labels_rl
c      call substr(runlog_header, outarr, mcol, kcol)
c      if (kcol.ne.ncol) stop 'SUBSTR mismatch: in runlog header'
      call read_runlog_header(lunr_ss_rl,data_fmt_read_rl,col_labels_rl)
      write(*,*)' read runlog header'

      ldot=lloc(runlog,'.')
      open(lunw_avg_rl,file='sa_'//runlog,status='unknown')
      call write_runlog_header(lunw_avg_rl,version,data_fmt_write_rl)
      write(*,*)' wrote runlog header'
c      write(lunw_avg_rl,*)3,ncol
c      write(lunw_avg_rl,'(a)') version
c      write(lunw_avg_rl,'(a)') col_labels_rl(:lnbc(col_labels_rl))

c  Open a seperate report file for each spectrum type
      rptfile='spectrum_averaging.rpt#'
      lr=index(rptfile,'#')
      do jtype=1,ntype
         write(rptfile(lr:lr),'(i1)') jtype
         open(lun_rpt+jtype,file=rptfile,status='unknown')
      end do
      istat=0
      do while (istat.eq.0)  ! Loop over average spectra
         do jtype=1,ntype  ! Loop over spectrum types (e.g. HgCd/InSb)

c        Compute average
c           write(*,*)'calling rravgcom: mode=1'
            write(*,*)'Reading single spectra...'
            call rravgcom(1,nus,nue,
     &      lunr_ss_rl,lunw_avg_rl,lun_rpt,
     &      data_fmt_read_rl,data_fmt_write_rl,
     &      ntype,jtype,nspe,krec,istat)
            if(jtype.eq.1 .and. istat.ge.1) exit  ! hit EOF
c            write(*,*)'called rravgcom 1',istat,krec

            call skiprec(lunr_ss_rl,-krec)

c        Compare with average
c            write(*,*)'calling rravgcom: mode=2'
            call rravgcom(2,nus,nue,
     &      lunr_ss_rl,lunw_avg_rl,lun_rpt,
     &      data_fmt_read_rl,data_fmt_write_rl,
     &      ntype,jtype,nspe,krec,istat)
c            write(*,*)'called rravgcom 2',istat

            write(*,*)
            if(jtype.lt.ntype) call skiprec(lunr_ss_rl,-krec)
         end do  !  jtype=1,ntype
      end do  !  do while loop over average spectra

      close(lunr_ss_rl)
      close(lunw_avg_rl)
      do jtype=1,ntype
         close(lun_rpt+jtype)
      end do
 
c  Merge multiple .rpt files into one file
      write(*,*)' Merging .rpt files into single report file'
      open(lun_rpt,file=runlog(:ldot)//'rpt',status='unknown')
      rptfile='spectrum_averaging.rpt#'
      lr=index(rptfile,'#')
      do jtype=1,ntype
         write(rptfile(lr:lr),'(i1)') jtype
         open(lun_rpt+jtype,file=rptfile,status='unknown')
      end do
      do   ! loop over records
         lst=1
         do jtype=1,ntype
            read(lun_rpt+jtype,'(a)',end=89) stuff(lst:)
            lst=lnbc(stuff)+5
         end do
         write(lun_rpt,'(a)') stuff(:lst)
      end do
89    continue
      do jtype=1,ntype
         close(lun_rpt+jtype,status='delete')
      end do
      close(lun_rpt)
      stop
      end
