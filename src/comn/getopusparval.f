      subroutine getopusparval(luns,ipoint,param,iend,mrs,i2val,rs)

c  Reads an OPUS file and extracts the value of the parameter PARAM
c  from the parameter block starting at IPOINT.
c
c  Inputs:
c          LUNS   I*4   Logical Unit Number of file (already open)
c          IPOINT I*4   Pointer (bytes)
c          PARAM  C*3   Name of parameter (e.g. "NPT")
c          IEND   I*4   Endian Flag (+-1)
c          MRS    I*2   Maximum value for RS
c
c  Outputs:
c          I2VAL(4)  I*2  Parameter value
c          RS        I*4  Length of parameter (bytes)
c
c On exit, rs=0 indicates failure to find PARAM
c          rs>0 indicates success, with the length of data field = 2*RS
c
      implicit none

      integer*4
     & luns,       ! Subroutine input argument (see above)
     & ipoint,     ! Subroutine input argument (see above)
     & iend,       ! Subroutine input argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & reclen,     ! Record length in bytes for file reads
     & irec,       ! Record number during file reads
     & is          ! Loop index during read of parameter value

      parameter (bigendian=1,reclen=2)
      integer*2 mrs,i2val(mrs),typ,rs
      character head*4, param*3

      head='....'
c  Search Parameter Block for PARAM label.
      irec=ipoint/reclen
      do while (head(:3).ne.param)
         read(luns,rec=1+irec) head(1:2)
         read(luns,rec=2+irec) head(3:4)
         read(luns,rec=3+irec) typ
         read(luns,rec=4+irec) rs
         if(iend.eq.bigendian) then
            call rbyte(rs,2,1)
            call rbyte(typ,2,1)
         endif
c         write(*,*)'irec,typ,rs=',irec,typ,rs,' ',head,' ',param
         if(rs.gt.mrs) then
            write(*,*) 'Increase parameter MRS to ',rs
            rs=0
         endif
         if(rs.le.0) return  ! failed to find PARAM
         irec=irec+4+rs
      end do  !  while (head(:3).ne.param)
c
c  Read the value of PARAM
      irec=irec-rs
      do is=1,rs
         irec=1+irec
         read(luns,rec=irec) i2val(is)
      end do
c
c  Byte-reverse Integers & Reals
      if(iend.eq.bigendian .and. typ.le.1) call rbyte(i2val,2*rs,1)
      return    ! success
      end

c  Bruker Parameter Types are:
c  0 = INT32
c  1 = REAL64
c  2 = STRING
c  3 = ENUM
c  4 = SENUM
