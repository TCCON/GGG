c
c  These are constants used to create/read the header of OPUS files.
c  Except for 'mdb' and 'mhl', the names are taken from the Bruker
c  documentation or are selected to closely match the Bruker names.
c
      integer*4
     & mdb,        ! Maximum number of directory blocks
     & mhl,        ! Maximum Header Length, in bytes
     & magic,      ! Magic number for OPUS format
     & dbbampl,    ! Bruker part of block type for amplitude data
     & dbbsamp,    ! Bruker part of block type for sample data
     & dbbdstat,   ! Bruker part of block type for data status parameter
     & dbbinstr,   ! Bruker part of block type for instrument status parameter
     & dbbaqpar,   ! Bruker part of block type for standard acquisition param
     & dbbftpar,   ! Bruker part of block type for Fourier transform param
     & dbboptpar,  ! Made-up part of block type for optics parameter
     & dbborgpar,  ! Bruker part of block type for sample origin parameter
     & dbbspec,    ! Bruker part of block type for spectrum data
     & dbbigrm,    ! Bruker part of block type for interferogram data
     & dbbdir,     ! Bruker part of block type for directory
     & dbbslav,    ! Made-up part of block type for slave data
     & typ_i4,     ! Parameter type of I*4 (INT32)
     & typ_r8,     ! Parameter type of R*8 (REAL64)
     & typ_string, ! Parameter type of string
     & typ_enum,   ! Parameter type of string "enum"
     & typ_senum   ! Parameter type of string "senum"

      real*8
     & filevers    ! Expected file version date

      parameter (mdb=32, mhl=4096)
      parameter (magic=-16905718)
      parameter (filevers=920622.0d0) ! 1992-06-22
      parameter (dbbampl  = 3*(2**0))
      parameter (dbbsamp  = 1*(2**2))
      parameter (dbbdstat = 1*(2**4))
      parameter (dbbinstr = 2*(2**4))
      parameter (dbbaqpar = 3*(2**4))
      parameter (dbbftpar = 4*(2**4))
      parameter (dbboptpar= 6*(2**4)) ! Bruker doc calls this "dbbprcpar"
      parameter (dbborgpar=10*(2**4))
      parameter (dbbspec  = 1*(2**10))
      parameter (dbbigrm  = 2*(2**10))
      parameter (dbbdir   =13*(2**10))
      parameter (dbbslav  =32*(2**10))! Reverse-engineered, not a Bruker name
      parameter (typ_i4     = 0)
      parameter (typ_r8     = 1)
      parameter (typ_string = 2)
      parameter (typ_enum   = 3)
      parameter (typ_senum  = 4)
