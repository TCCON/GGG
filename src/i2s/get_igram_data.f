      subroutine get_igram_data(inpath,fpfn,slice,slice1,slicen,
     & verbose,msl,mip,mch,nptvec,bpdata,ichan,counter,buf_igram)
c
c  Input:
c    inpath         C*(*)  Directory path to input OPUS filea/slicess
c    fpfn           C*(*)  Name of input OPUS file
c    slice           I*4   Slice number from catalog, first slice of scan set
c    slice1          I*4   Starting slice number for this run
c    slicen          I*4   Ending slice number for this run
c    verbose         I*4   Level of verbosity for displayed messages
c    msl             I*4   Maximum number of interferogram slices per scan set
c    mip             I*4   Maximum number of input points
c    mch             I*4   Maximum number of data channels
c    nptvec(msl)     I*4   Number of PoinTs, from OPUS header
c    bpdata(msl,mch) I*4  Byte pointers into the data blocks of slices
c    ichan           I*4   Channel number (1=InGaAs=slave, 2=Si=master)
c
c  Output:
c    counter         I*4   Point count, data destination into 'buf_igram'
c    buf_igram(mip)  R*4   Data storage and FFT buf_igram
c
      implicit none

      integer*4
     & slice,        ! Subroutine input argument (see above)
     & slice1,       ! Subroutine input argument (see above)
     & slicen,       ! Subroutine input argument (see above)
     & verbose,      ! Subroutine input argument (see above)
     & msl,          ! Subroutine input argument (see above)
     & mip,          ! Subroutine input argument (see above)
     & mch,          ! Subroutine input argument (see above)
     & nptvec(msl),  ! Subroutine input argument (see above)
     & bpdata(msl,mch),! Subroutine input argument (see above)
     & ichan,        ! Subroutine input argument (see above)
     & counter,      ! Subroutine output argument (see above)
     & iendian,      ! Endianness of computer
     & bigendian,    ! Named constant for big endian detection 
     & luns,         ! Logical Unit Number for OPUS input file/slice
     & reclen,       ! Record length in bytes for single values
     & blklen,       ! Record length in bytes for data blocks
     & blkpnt,       ! Number of data points in a block
     & slicind,      ! Slice Index
     & lnbc,         ! Integer function Last Non-Blank Character in string
     & lf,
     & npts,         ! Number of points in OPUS input file/slice
     & inpoint,      ! Byte pointer to data values in OPUS file/slice
     & rs1,          ! Starting record number single value reads at start
     & re1,          ! Ending record number single value reads at start
     & rs2,          ! Starting record number block reads
     & re2,          ! Ending record number block reads
     & rs3,          ! Starting record number single value reads at end
     & re3,          ! Ending record number single value reads at end
     & recnum,       ! Number of current data record from OPUS file
     & ijk        ! Index into buf_igram during block read and byte swapping

      parameter (bigendian=1)
      parameter (luns=21)
      parameter (reclen=4)
      parameter (blklen=2**9)
      parameter (blkpnt=blklen/reclen)

      real*4 buf_igram(mip) ! Subroutine output argument (see above)

      character
     & inpath*(*),  ! Subroutine input argument (see above)
     & fpfn*(*),    ! Subroutine input argument (see above)
     & filename*256 ! Full file name for OPUS input fila/slicese
c
c  Initialize variables.
      call getendian(iendian)
      counter=0

c  Loop over the slices that make up this run.
      do slicind=slice1,slicen
         inpoint=bpdata(slicind-slice+1,ichan)
         npts=nptvec(slicind-slice+1)

         write(filename,'(2a)')inpath(:lnbc(inpath)),fpfn(:lnbc(fpfn))
         lf=lnbc(filename)
         if(slice1.gt.0) write(filename(lf+1:),'(i0,a2)') slicind,'.0'

         if(verbose.ge.4) write(*,'(a,a7,i0,a7,i0)') filename(:lf),
     &  '  npts=',npts,'  byte=',inpoint

c  We pre-compute the limits of three zones for reading the data.  Zone 1
c  is read in one-point increments until we get to a block boundary.
c  Zone 2 is read in full blocks.  Zone 3 goes back to one-point reads
c  to the end of the data section.
c
c  First compute rsX and reX in units of bytes, numbered from 0.
         rs1=inpoint
         rs2=(inpoint/blklen)*blklen
         if(rs2.lt.inpoint) rs2=rs2+blklen
         rs3=((inpoint+(npts*reclen))/blklen)*blklen

         re1=rs2-1
         re2=rs3-1
         re3=(inpoint+(npts*reclen))-1 
c
c  Now convert rsX and reX to blocks, numbered from 1.
         rs1=1+(rs1/reclen)
         re1=1+(re1/reclen)
         rs2=1+(rs2/blklen)
         re2=1+(re2/blklen)
         rs3=1+(rs3/reclen)
         re3=1+(re3/reclen)
c
c  Convert sample count to a pointer into an array numbered from 1.
         counter=counter+1
c
c  Perform file read in three zones (points/blocks/points).
         open(luns,file=filename,form='unformatted',status='old',
     &   access='direct',recl=reclen) 
         do recnum=rs1,re1
            read(luns,rec=recnum) buf_igram(counter)
            counter=counter+1
         enddo    ! recnum=rs1,re1
         close(luns)

         open(luns,file=filename,form='unformatted',status='old',
     &   access='direct',recl=blklen) 
         do recnum=rs2,re2
            read(luns,rec=recnum)
     &      (buf_igram(ijk),ijk=counter,counter+blkpnt-1)
            counter=counter+blkpnt
         enddo    ! recnum=rs2,re2
         close(luns)

         open(luns,file=filename,form='unformatted',status='old',
     &   access='direct',recl=reclen) 
         do recnum=rs3,re3
            read(luns,rec=recnum) buf_igram(counter)
            counter=counter+1
         enddo    ! recnum=rs3,re3
         close(luns)
c
c  Convert pointer back to sample count.
         counter=counter-1
c
c  On big-endian computer (e.g. SPARC and G4) swap data bytes.
         if(iendian.eq.bigendian) then
c            do ijk=counter-npts+1,counter
c               call r4rev(buf_igram(ijk))
c            enddo
            call rbyte(buf_igram(counter-npts+1),4,npts)
         endif

      end do ! do slicind=slice1,slicen

      return
      end
