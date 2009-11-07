c  Program spectrum_average.f

c  Averages a series of binary spectra.

c  Works for any binary spectra having a GGG-format runlog,
c  which is adapted for use as the input file.
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

      integer*4 lst,
     & lun_rlg,     ! LUN to read input runlogs from
     & lun_rpt,     ! LUN to write rms deviation
     & lun_wbs,     ! LUN to Write Binary Spectra (averages)
     & lun_wav,     ! LUN to write acerage spectrum to
     & lr,nhl,      ! number of header lines
     & mcol, ncol, kcol, ! Number of columns in runlog
     & lnbc,        ! function Last Non-Black Character
     & ntype,jtype, ! Periodicity of spectra to be averaged
     & irec,nrec,   ! number of records
     & nspe         ! number of single spectra being averaged

      parameter (lun_rpt=30,lun_rlg=15,lun_wbs=16,lun_wav=17)
      parameter (mcol=40)

      character 
     & version*56,   ! version
     & rptfile*40,
     & stuff*256,
     & runlog*120,outarr(mcol)*20,runlog_header*512

      integer*4
     & krec,
     & istat         ! status flag (0=success, 1=EOF)

      real*8
     & nus, nue      ! selected frequency range of interest

      version=' spectrum_averaging   Version 0.1.1   30-Oct-2009   GCT '
      write(6,*) version

      write(*,*)'Enter Starting & Ending frequencies (cm-1):'
      write(*,*)'Enter 0 99999 to retain original spectral limits'
      read(*,*) nus,nue
c
      write(*,*)'Enter path to input_file/runlog:'
      read(*,'(a)') runlog

      write(*,*)' Enter number of spectrum types interleaved in runlog:'
      write(*,*)' (e.g., HgCd/InSb=2; InGaAs/Si=2; InGaAs-only=1)'
      read(*,*) ntype

      open(lun_rlg,file=runlog,status='old')
      read(lun_rlg,*) nhl,ncol
      call skiprec(lun_rlg,nhl-2)
      read(lun_rlg,'(a)') runlog_header
      call substr(runlog_header, outarr, mcol, kcol)
      if (kcol.ne.ncol) stop 'mismatch: runlog header'

      open(lun_wav,file='spectrum_averaging.grl',status='unknown')
      write(lun_wav,*)3,ncol
      write(lun_wav,'(a)') version
      write(lun_wav,'(a)') runlog_header(:lnbc(runlog_header))

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
         call rravgcom(1,nus,nue,
     &   lun_rlg,lun_wav,lun_rpt,
     &   ntype,jtype,nspe,krec,istat)

         call skiprec(lun_rlg,-krec)

c        Compare with average
         call rravgcom(2,nus,nue,
     &   lun_rlg,lun_wav,lun_rpt,
     &   ntype,jtype,nspe,krec,istat)

         write(*,*)
         if(jtype.lt.ntype) call skiprec(lun_rlg,-krec)
      end do  !  jtype=1,ntype
      end do  !  Loop over average spectra

99    close(lun_rlg)
      close(lun_wav)
      do jtype=1,ntype
        close(lun_rpt+jtype)
      end do
 
c  Merge multiple .rpt files into one file
      open(lun_rpt,file='spectrum_averaging.rpt',status='unknown')
      rptfile='spectrum_averaging.rpt#'
      lr=index(rptfile,'#')
      do jtype=1,ntype
         write(rptfile(lr:lr),'(i1)') jtype
         open(lun_rpt+jtype,file=rptfile,status='unknown')
      end do
      do irec=1,999999
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
