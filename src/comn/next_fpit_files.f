       subroutine next_fpit_files(jul, zpdtim, oblat, oblon,
     &     modname, vmrname)
c NEXT_FPIT_FILES(JUL, ZPDTIM, OBLAT, OBLON, MODNAME, VMRNAME)
c   Generate the base names of the .mod and .vmr files required
c   for the next observation
c 
c Inputs:
c   JUL - julian day number from the JULIAN subroutine
c   ZPDTIM - zero-path difference time, read from a runlog
c   OBLAT - observation latitude, read from a runlog
c   OBLON - observation longitude, read from a runlog
c
c Outputs:
c   MODNAME - the basename of the .mod file (does not include directory)
c   VMRNAME - the basename of the .vmr file (does not include directory)
c
c Josh Laughner, 29 Oct 2019

       implicit none
      
       integer*4 
     & jul,       ! julian day number (input)
     & iyyyy,     ! year number
     & imm,       ! month number
     & idd,       ! day number
     & ihh        ! hour number

       real*8 
     & oblat,     ! observation latitude (input)
     & oblon,     ! observation longitude (input)
     & zpdtim     ! time of ZPD (UT hours)

       character
     & modname*80,         ! new name of the model file (output)
     & vmrname*80,         ! new name of the vmr file (output)
     & ns*1,               ! north or south
     & ew*1                ! east or west

c----- Constants ------
       integer*4, parameter :: nhr = 3 ! number of hrs between FPIT files
c----------------------

       if(oblat.ge.0.0) then 
          ns='N'
       else
          ns='S'
       endif
       if(oblon.ge.0.0) then 
          ew='E'
       else
          ew='W'
       endif

       call caldat(jul,iyyyy,imm,idd)
       ihh = nhr*nint(zpdtim/nhr)
       if(ihh.lt.0) then 
          ihh=ihh+24
          call caldat(jul-1, iyyyy, imm, idd) 
       endif
       if(ihh.ge.24) then 
          ihh=ihh-24
          call caldat(jul+1, iyyyy, imm, idd) 
       endif

       write(modname,
     &  '(a5,i4.4,2i2.2,i2.2,a2,i2.2,2a1,i3.3,a1,a4)')
     &  'FPIT_',iyyyy,imm,idd,ihh,'Z_',nint(abs(oblat)),ns,
     &  '_',nint(abs(oblon)),ew,'.mod'
       write(vmrname,
     &  '(a4,i4.4,3i2.2,a2,i2.2,a1,a1,i3.3,a1,a4)')
     &  'JL1_',iyyyy,imm,idd,ihh,'Z_',nint(abs(oblat)),ns,
     &  '_',nint(abs(oblon)),ew,'.vmr'

       end subroutine
