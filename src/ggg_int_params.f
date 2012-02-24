
      integer*4
     & mvmode,
     & ncell,
     & mspeci,
     & mauxcol,
     & mcolvav,
     & mlev,
     & mrow_qc,
     & nchar,
     & mpath,
     & mfilepath

      parameter(mvmode=30)      !maximum number of vibrational modes
      parameter(ncell=2)        !number of gas cells in solar beam
      parameter(mspeci=160)     !maximum number of different species listed in ISOTOPOLOG.DAT
      parameter(mauxcol=25)     !maximum number of auxiliary parameters/columns
      parameter(mcolvav=150)    !maximum number of columns in .vav file
      parameter(mlev=250)       !maximum number of atmospheric levels
      parameter(mrow_qc=150)    !maximum number of rows in qc file
      parameter(nchar = 57)     !maximum length of character string for spectrum names
      parameter(mpath = 128)    !maximum length of character string for path variables
      parameter(mfilepath = mpath + 40)    !maximum length of absolute path of file names

